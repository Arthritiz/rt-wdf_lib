/*
 ==============================================================================

 This file is part of the RT-WDF library.
 Copyright (c) 2015,2016 - Maximilian Rest, Ross Dunkel, Kurt Werner.

 Permission is granted to use this software under the terms of either:
 a) the GPL v2 (or any later version)
 b) the Affero GPL v3

 Details of these licenses can be found at: www.gnu.org/licenses

 RT-WDF is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 -----------------------------------------------------------------------------
 To release a closed-source product which uses RT-WDF, commercial licenses are
 available: write to rt-wdf@e-rm.de for more information.

 ==============================================================================

 rt-wdf_nlSolvers.cpp
 Created: 2 Dec 2015 4:08:19pm
 Author:  mrest

 ==============================================================================
 */

#include "rt-wdf_nlSolvers.h"

#ifdef _WIN32
#include <Windows.h>
void debugOutput(const std::string& msg) {
    OutputDebugString(msg.c_str());
}
#endif

//==============================================================================
// Parent class for nlSolvers
//==============================================================================
nlSolver::nlSolver( ) :
#ifdef TRACKING
                pTracker("p"),
#endif
                numNLPorts( 0 ) {

}

nlSolver::~nlSolver( ) {

}

//----------------------------------------------------------------------
int nlSolver::getNumPorts( ) {
    return numNLPorts;
}


//==============================================================================
// Newton Solver
//==============================================================================
nlNewtonSolver::nlNewtonSolver( std::vector<nlModel*> nlList,
                        matData* myMatData ) :
                                myMatData ( myMatData )  {

    //std::cout << "newton tol: " << TOL << std::endl;

    nlModels = nlList;

    numNLPorts = 0;
    for ( nlModel* model : nlModels ) {
        numNLPorts += model->getNumPorts();
    }

    VEC_ALLOC(x0       , numNLPorts);
    VEC_ALLOC(F        , numNLPorts);
    MAT_ALLOC(J        , numNLPorts, numNLPorts);
    VEC_ALLOC(fNL      , numNLPorts);
    MAT_ALLOC(JNL      , numNLPorts, numNLPorts);
    VEC_ALLOC(Fmat_fNL , numNLPorts);
    VEC_ALLOC(Emat_in  , numNLPorts);
    MAT_ASSIGN_ID_BY_MAT(idJNL, *JNL);

#ifdef TRACKING
    pTracker.init(numNLPorts);
#endif
}

nlNewtonSolver::~nlNewtonSolver( ) {
    size_t modelCount = nlModels.size();
    for( size_t i = 0; i < modelCount; i++ ) {
        delete nlModels[i];
    }
    delete x0;
    delete F;
    delete J;
    delete fNL;
    delete JNL;
    delete Fmat_fNL;

#ifdef TRACKING
    //std::cout << "totalIter: " << totalIter << ", callCount: " << callCount << std::endl;
    //std::cout << "avgIter(vanilla): " << totalIter/(float)callCount << std::endl;
    avgIter = (avgIter * callCount/STEP + totalSubIter/(float)STEP)/(callCount/(float)STEP);
    //std::cout << "avgIter(my way): " << avgIter << std::endl;

    std::cout << totalIter/(float)callCount << "," << avgIter;

    //std::cout << "solver elapsed: " << totalElapsed << std::endl;
    //pTracker.print();
#endif
}

void nlNewtonSolver::iterWay()
{
    int iter = 0;            // # of iteration
    FloatType alpha = 1.0;

    MAT_SETZERO((*J));

    if ( firstRun ) {
        firstRun = false;
    }
    else {
#ifdef PREV_WAY
        (*x0) = (*Fmat_fNL) + myMatData->Emat*prevInWaves;
#else
        (*x0) = (*Fmat_fNL) + (*Emat_in);
#endif
    }

    evalNlModels( myMatData, x0 );

    FloatType normF = W_NORM(*F);
    //printf("iter alpha         ||F||_2\n");
    //printf(" %3g %9.2e %14.7e\n", iter, alpha, normF);

    Wvec xnew;
    FloatType normFnew;

#ifdef TRACKING
    auto start = std::chrono::high_resolution_clock::now();
#endif

    while ( (normF >= TOL) && (iter < ITMAX) )
    {
        Wvec p = - (*J).i() * (*F);
        xnew = (*x0) + alpha * p;
        evalNlModels( myMatData, &xnew);
        normFnew = W_NORM(*F);

#ifdef BTWAY
        if (normFnew < normF)
        {
            alpha = 1.0;
        } else
        {
            alpha /= 2;
        }
#endif

        (*x0) = xnew;
        normF = normFnew;

        iter++;

    //        printf(" %3g %9.2e %14.7e\n", iter, alpha, normF);
    }

#ifdef TRACKING
    totalElapsed += (double)std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count()/1e6;
#endif

    if (iter == ITMAX)
    {
        std::cout << "convergence failed" << std::endl;
        throw std::runtime_error("convergence failed");
    }

#ifdef TRACKING
    totalIter += iter;
    callCount += 1;
    totalSubIter += iter;

    if ( 0 == callCount % STEP )
    {
        avgIter = (avgIter * (callCount/STEP - 1) + totalSubIter/(float)STEP)/(callCount/STEP);
        totalSubIter = 0;
    }
#endif

//#ifdef _WIN32
//    std::ostringstream oss;
//    oss.precision(6);
//    oss << std::scientific << normF;
//    std::string resSN = oss.str();
//
//    debugOutput("iter: " + std::to_string(iter) + ", normF: " + resSN + "\n");
//#endif
}

FloatType interpolation(int startIndex, int endIndex, FloatType startVal, FloatType endVal, int currentIndex)
{
 return (endVal - startVal)*(currentIndex - startIndex)/(endIndex - startIndex) + startVal;
}

Wvec nlNewtonSolver::interpoCore(int p0Index, int p1Index, int p2Index, int fuzzIndex)
{
    int p1Count = std::get<2>(pRangeList[1]);
    int p2Count = std::get<2>(pRangeList[2]);

    int essenceStartIndex = fuzzIndex*(p1Count*p2Count*dimSize) + p1Index*(p2Count*dimSize) + p2Index*dimSize;


    Wvec iVec;
    iVec.set_size(fNL->n_elem);

    for (int relI = 0; relI< dimSize; relI++)
    {
        Essence& essence = essenceList[essenceStartIndex+relI];

        FloatType val;

        if (p0Index <= essence.firstStartKneeIndex) // compare with first start knee index
        {
           val = interpolation(0, essence.firstStartKneeIndex, essence.startVal, essence.firstStartKneeVal, p0Index);
        } else if (p0Index <= essence.firstEndKneeIndex)
        {
           val = interpolation(essence.firstStartKneeIndex, essence.firstEndKneeIndex, essence.firstStartKneeVal, essence.firstEndKneeVal, p0Index);
        } else
        {
           val = interpolation(essence.firstEndKneeIndex, std::get<2>(pRangeList[0]) - 1, essence.firstEndKneeVal, essence.endVal, p0Index);
        }

        iVec(relI) = val;
    }

    return iVec;
}

void nlNewtonSolver::interpoTabWay()
{
    if (0 == dimSize)
    {
        throw;
    }

    for (int i = 0; i < pRangeList.size(); i++)
    {
        if((*Emat_in)(i) > std::get<1>(pRangeList[i]))
        {
            (*Emat_in)(i) = std::get<1>(pRangeList[i]);
        }

        if((*Emat_in)(i) < std::get<0>(pRangeList[i]))
        {
            (*Emat_in)(i) = std::get<0>(pRangeList[i]);
        }
    }

    std::vector<int> pIndexList;

    for (int i = 0; i < pRangeList.size(); i++)
    {
        RangeInfo& rangeInfo = pRangeList[i];
        FloatType rawIndex = (std::get<2>(rangeInfo) - 1)*((*Emat_in)(i) - std::get<0>(rangeInfo))/(std::get<1>(rangeInfo) - std::get<0>(rangeInfo));

        if (rawIndex - std::round(rawIndex) > 0.5)
        {
            pIndexList.push_back(std::round(rawIndex) + 1);
        } else
        {
            pIndexList.push_back(std::round(rawIndex));
        }
    }

    FloatType rawFuzzIndex = (std::get<2>(fuzzRange) - 1)*(fuzz - std::get<0>(fuzzRange))/(std::get<1>(fuzzRange) - std::get<0>(fuzzRange));

    FloatType lowerFuzzIndex = std::round(rawFuzzIndex);
    Wvec lowerIVec = interpoCore(pIndexList[0], pIndexList[1], pIndexList[2], int(lowerFuzzIndex));

    //if (int(rawFuzzIndex) >= std::get<2>(fuzzRange) - 1)
    //{
    //    *fNL = lowerIVec;
    //    return;
    //}

    FloatType upperFuzzIndex = lowerFuzzIndex + 1;
    Wvec upperIVec = interpoCore(pIndexList[0], pIndexList[1], pIndexList[2], int(upperFuzzIndex));

    *fNL = (rawFuzzIndex - lowerFuzzIndex)*upperIVec + (upperFuzzIndex - rawFuzzIndex)*lowerIVec;
}

//----------------------------------------------------------------------
void nlNewtonSolver::nlSolve( Wvec* inWaves,
                          Wvec* outWaves ) {
#ifndef RECORD_TABLE
    *Emat_in = (myMatData->Emat)*(*inWaves);
#endif

#ifdef TRACKING
    pTracker.track(*Emat_in);
#endif

    iterWay();
    //interpoTabWay();

    (*outWaves) = (myMatData->Mmat) * (*inWaves) + (myMatData->Nmat) * (*fNL);

#ifdef RECORD_TABLE
    pVec.push_back(*Emat_in);
    iVec.push_back(*fNL);
#endif

#ifdef PREV_WAY
    prevInWaves = *inWaves;
#endif
}

//----------------------------------------------------------------------
void nlNewtonSolver::evalNlModels( matData* myMatData,
                                   Wvec* x ) {
    int currentPort = 0;
    (*JNL).zeros();


    for ( nlModel* model : nlModels ) {
       model->calculate( fNL, JNL, x, this->x0, &currentPort );
    }

    (*Fmat_fNL) = myMatData->Fmat*(*fNL);
    (*F) = (*Emat_in) + (*Fmat_fNL) - (*x);
    (*J) = (myMatData->Fmat)*(*JNL) - idJNL;
}

//nlTabSolver::nlTabSolver(std::vector<nlModel*> nlList, matData* myMatData, std::vector<std::tuple<double, double, int>> mapMeta): myMatData ( myMatData )
nlTabSolver::nlTabSolver(std::vector<nlModel*> nlList, matData* myMatData): myMatData ( myMatData )
{
    //std::vector<std::tuple<double, double, int>> mapMeta = {
    //       {0.0, 1.0, 2},
    //    {0.0, 1.0, 2},
    //    {0.0, 1.0, 3},
    //    {0.0, 1.0, 4}};

    std::vector<std::tuple<double, double, int>> mapMeta = {
        { -0.0594, 0.2334, 20 },
    { -1.3009, 0.2197, 50 },
    { -0.5825, 0.2763, 30 },
    { -8.9857, 0.2697, 100 }};

    std::vector<std::pair<double, double>> currentLimList = {
        {  7.7168e-09,  3.0183e-04 },
    { -2.7231e-04, -1.2915e-08 },
    { -1.0067e-09,  1.2291e-03 },
    { -9.7968e-04, -6.0400e-09 },
    };

    nlModels = nlList;

    numNLPorts = 0;
    for ( nlModel* model : nlModels ) {
        numNLPorts += model->getNumPorts();
    }

    int totalCount = 1;

    for (auto& tup: mapMeta)
    {
        totalCount *= std::get<2>(tup);
    }

    for (int i = 0; i < totalCount; i++)
    {
        double start, end;
        int count;
        int index;
        double factor = 1.0;

        Wvec vVec(numNLPorts);
        for (int j = 0; j < numNLPorts; j++)
        {
            count = std::get<2>(mapMeta[j]);
            index = (int)(i / factor) % count;

            start = std::get<0>(mapMeta[j]);
            end = std::get<1>(mapMeta[j]);

            vVec(j) = index*(end - start)/count + start;

            factor *= count;
        }

        int start_index = 0;
        Wvec iVec(numNLPorts);
        for ( nlModel* model : nlModels )
        {
            model->getCurrents(vVec, iVec, start_index);
            start_index += model->getNumPorts();
        }

        // sifting
        bool outOfRange = false;
        for (int j = 0; j < iVec.size(); j++)
        {
            if ( (iVec[j] > currentLimList[j].second) || (iVec[j] < currentLimList[j].first) )
            {
                outOfRange = true;
                break;
            }
        }

        if (outOfRange)
        {
            continue;
        }

        vsVec.push_back(vVec);
        isVec.push_back(iVec);
    }

    std::cout << "totalCount: " << totalCount << std::endl;
    std::cout << "after sifting: " << isVec.size() << std::endl;

    //std::cout << "-v-" << std::endl;
    //for (auto& vVec: vsVec)
    //{
    //    vVec.t().print();
    //    std::cout << "--" << std::endl;
    //}

    //std::cout << "-i-" << std::endl;
    //for (auto& iVec: isVec)
    //{
    //    iVec.t().print();
    //    std::cout << "--" << std::endl;
    //}

    pMat.resize( numNLPorts, isVec.size() );
}

nlTabSolver::~nlTabSolver()
{

}

void nlTabSolver::resetTab()
{
    Wmat vsMat(numNLPorts, vsVec.size());
    Wmat isMat(numNLPorts, isVec.size());

    int theCount = vsVec.size();
    for (int i = 0; i < theCount; i++)
    {
        vsMat.col(i) = vsVec[i];
        isMat.col(i) = isVec[i];
    }

    Wmat tmp = myMatData->Fmat*isMat;
    pMat = vsMat - tmp;
}

void nlTabSolver::nlSolve( Wvec* inWaves,
                          Wvec* outWaves )
{

    Wvec p = myMatData->Emat * (*inWaves);

    double minDis = arma::norm(p - pMat.col(0), 2);
    int minIndex = 0;

    for (int i = 1; i < pMat.n_cols; i++)
    {
        double tmp = arma::norm(p - pMat.col(i), 2);
        if (tmp < minDis)
        {
            minDis = tmp;
            minIndex = i;
        }
    }

    (*outWaves) = (myMatData->Mmat) * (*inWaves) + (myMatData->Nmat) * isVec[minIndex];
}
