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
}

//----------------------------------------------------------------------
void nlNewtonSolver::nlSolve( Wvec* inWaves,
                          Wvec* outWaves ) {

    int iter = 0;            // # of iteration
    FloatType alpha = 0;

    *Emat_in = (myMatData->Emat)*(*inWaves);

    (*J).zeros();

    if ( firstRun ) {
        firstRun = false;
    }
    else {
        (*x0) = (*Fmat_fNL) + (*Emat_in);
    }

    evalNlModels( inWaves, myMatData, x0 );

    FloatType normF = arma::norm(*F);
    //printf("iter alpha         ||F||_2\n");
    //printf(" %3g %9.2e %14.7e\n", iter, alpha, normF);

    Wvec xnew;
    FloatType normFnew;

    while ( (normF >= TOL) && (iter < ITMAX) )
    {
        Wvec p = - (*J).i() * (*F);
        alpha = 1;
        xnew = (*x0) + alpha * p;
        evalNlModels(inWaves, myMatData, &xnew);
        normFnew = arma::norm(*F);

        (*x0) = xnew;
        normF = normFnew;
        iter++;

    //        printf(" %3g %9.2e %14.7e\n", iter, alpha, normF);
    }

    if (iter == ITMAX)
    {
        std::cout << "convergence failed" << std::endl;
    }

//#ifdef _WIN32
//    std::ostringstream oss;
//    oss.precision(6);
//    oss << std::scientific << normF;
//    std::string resSN = oss.str();
//
//    debugOutput("iter: " + std::to_string(iter) + ", normF: " + resSN + "\n");
//#endif

    (*outWaves) = (myMatData->Mmat) * (*inWaves) + (myMatData->Nmat) * (*fNL);

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

