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
nlSolver::nlSolver( ) : numNLPorts( 0 ) {

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
                        matData* myMatData ) : myMatData ( myMatData ) {

    std::cout << "newton tol: " << TOL << std::endl;

    nlModels = nlList;

    numNLPorts = 0;
    for ( nlModel* model : nlModels ) {
        numNLPorts += model->getNumPorts();
    }

    x0       = new Wvec( numNLPorts, arma::fill::zeros);
    F        = new Wvec( numNLPorts, arma::fill::zeros);
    J        = new Wmat( numNLPorts,numNLPorts, arma::fill::zeros);
    fNL      = new Wvec( numNLPorts, arma::fill::zeros);
    JNL      = new Wmat( numNLPorts,numNLPorts, arma::fill::zeros);
    Fmat_fNL = new Wvec( numNLPorts, arma::fill::zeros);
    Emat_in = new Wvec( numNLPorts, arma::fill::zeros);
    idJNL = arma::eye<Wmat>(size(*JNL));
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
    std::cout << "totalIter: " << totalIter << ", callCount: " << callCount << std::endl;
    std::cout << "avgIter(vanilla): " << totalIter/(float)callCount << std::endl;

    avgIter = (avgIter * callCount/STEP + totalSubIter/(float)STEP)/(callCount/(float)STEP);
    std::cout << "avgIter(my way): " << avgIter << std::endl;

    std::cout << "solver elapsed: " << totalElapsed << std::endl;
#endif
}

//----------------------------------------------------------------------
void nlNewtonSolver::nlSolve( Wvec* inWaves,
                          Wvec* outWaves ) {

    int iter = 0;            // # of iteration
    FloatType alpha = 1.0;

    *Emat_in = (myMatData->Emat)*(*inWaves);

    (*J).zeros();

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

    evalNlModels( inWaves, myMatData, x0 );

    FloatType normF = arma::norm(*F);
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
        evalNlModels(inWaves, myMatData, &xnew);
        normFnew = arma::norm(*F);

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

    (*outWaves) = (myMatData->Mmat) * (*inWaves) + (myMatData->Nmat) * (*fNL);

#ifdef PREV_WAY
    prevInWaves = *inWaves;
#endif
}

//----------------------------------------------------------------------
void nlNewtonSolver::evalNlModels( Wvec* inWaves,
                               matData* myMatData,
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
