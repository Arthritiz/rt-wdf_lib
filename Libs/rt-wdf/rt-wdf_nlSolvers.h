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

 rt-wdf_nlSolvers.h
 Created: 2 Dec 2015 4:08:19pm
 Author:  mrest

 ==============================================================================
 */

#ifndef RTWDF_NLSOLVERS_H_INCLUDED
#define RTWDF_NLSOLVERS_H_INCLUDED

//==============================================================================
#include <float.h>

#include "rt-wdf_types.h"
#include "rt-wdf_nlModels.h"
#include <limits>


//#define PREV_WAY
//#define TRACKING
//#define BTWAY
#define RECORD_TABLE

//==============================================================================
// Define enums for solver identifiers


// Iterative:
// TODO: introduce enums!
/** Enum to specify a Newton Solver*/
#define NEWTON_SOLVER   1



//==============================================================================
// Newton Solver config parameters

/** tolerance for ||F||_2 */
#define TOL     1.0e-4                     // TODO: evaluate physically meaningful tolerance.
/** limit on function evaluations */
#define ITMAX   2000


//==============================================================================
// Forward declarations
class nlSolver;
class nlNewtonSolver;


//==============================================================================
class nlSolver {

public:
    RangeInfo fuzzRange;
    std::vector<RangeInfo> pRangeList;
    std::vector<Essence> essenceList;
    FloatType fuzz;
    int dimSize = 0;

#ifdef RECORD_TABLE
    std::vector<Wvec> pVec;
    std::vector<Wvec> iVec;
#endif

#ifdef TRACKING
    RangeTracker pTracker;
#endif

    /** variable to store Emat * inWaves */
    Wvec* Emat_in;

    //----------------------------------------------------------------------
    /**
     Parent class for all non-linear solvers.
    */
    nlSolver();

    //----------------------------------------------------------------------
    /**
     Deconstructor.
    */
    virtual ~nlSolver();

    //----------------------------------------------------------------------
    /**
     Function which returns the number of ports on an NL model for memory house-
     keeping

     @returns                   the total number of non-linear ports which are
                                present at the root.
    */
    int getNumPorts();

    //----------------------------------------------------------------------
    /**
     Virtual function that processes a vector of incoming waves and
     returns a vector of outgoing waves according to the specified
     nonlinearities.

     Must be implemented in an NL solver according to the desired behaviour.

     @param *inWaves            is a pointer to a vector of incoming waves
     @param *outWaves           is a pointer to a vector of outgoing waves
     @param *myMatData          is a pointer to the E,F,M,N (and S) matrices
    */
    virtual void nlSolve( Wvec* inWaves,
                          Wvec* outWaves ) = 0;

    //----------------------------------------------------------------------
    /**
     Vector of enums that specify the types on non-linearities in the solver
    */
    std::vector<nlModel*> nlModels;
    /**
     Total number of non-linear ports which are present at the root.
    */
    int numNLPorts;

};


//==============================================================================
class nlNewtonSolver : public nlSolver {

protected:
    //----------------------------------------------------------------------
    /** struct which holds all root NLSS matrices including variable conversion */
    matData* myMatData;
    /** latest guess to solve the NLSS system */
    Wvec* x0;
    /** result of the NL equations */
    Wvec* fNL;
    /** Jacobian of the NL equations */
    Wmat* JNL;
    /** F vector for newton method */
    Wvec* F;
    /** J matrix for the newton method */
    Wmat* J;
    /** variable to store Fmat * fNL for x0 next prediction */
    Wvec* Fmat_fNL;

    Wmat idJNL;
    Wvec prevInWaves;

    /** flag to detect first run of the solver for a clean first initial guess */
    bool firstRun = true;

#ifdef TRACKING
    // tracking iteration
    long totalIter = 0;
    long totalSubIter = 0;
    float avgIter = 0.0;
    int callCount = 0;

    const int STEP = 100;

    // tracking elapse
    double totalElapsed = 0.0;
#endif

public:
    //----------------------------------------------------------------------
    /**
     Multi-dimensional Newton Solver class.

     Newton solver as introduced by Michael Jorgen Olsen:
     ("Resolving grouped nonlinearities in wave digital filters using
     iterative techniques"), DAFx-16

     Creates a newton solver and it's nonlinearities

     @param nlList              is a vector of enums that specify the types of
                                nonlinearities
     @param *myMatData          is a pointer to the E,F,M,N (and S) matrices
    */
    nlNewtonSolver( std::vector<nlModel*> nlList,
                matData* myMatData );

    /**
     Deconstructor.
    */
    ~nlNewtonSolver( );

    //----------------------------------------------------------------------
    /**
     Solver function that processes a vector of incoming waves and
     returns a vector of outgoing waves according to the specified
     nonlinearities.

     The actual solver operates on J and F matrices which are calculated
     by evalNlModels based on fNL, JNL, Emat and Fmat.

     The result outWaves is calculated based on inWaves, x, Nmat, Mmat.

     @param inWaves             is a pointer to a vector of incoming waves
     @param outWaves            is a pointer to a vector of outgoing waves
    */
    void nlSolve( Wvec* inWaves,
                  Wvec* outWaves );

    //----------------------------------------------------------------------
    /**
     Evaluates all non-linear model members of a solver and sets J and F
     members properly.

     It calculates J and F matrices based on both fNL, JNL as returned
     from nlModel.calculate(x,..) and Emat, Fmat from myMatData.

     @param *inWaves            is a pointer to a vector of incoming waves
     @param *myMatData          is a pointer to the E,F,M,N matrices
     @param *x                  is a pointer to the input values x
    */
    void evalNlModels(matData* myMatData, Wvec* x);

    void iterWay();

    // interpolation tab way
    Wvec interpoCore(int p0Index, int p1Index, int p2Index, int fuzzIndex);
    void interpoTabWay();
};

// [tab solver]: unfinished and deprecated
class nlTabSolver: public nlSolver
{
    matData* myMatData;
    std::vector<Wvec> vsVec;
    std::vector<Wvec> isVec;
    Wmat pMat;


public:
    //nlTabSolver( std::vector<nlModel*> nlList,
    //            matData* myMatData,
    //            std::vector<std::tuple<double, double, int>> mapMeta );
    nlTabSolver( std::vector<nlModel*> nlList,
                matData* myMatData );

    ~nlTabSolver();

    void nlSolve( Wvec* inWaves,
                  Wvec* outWaves );

    void resetTab();
};

#endif  // RTWDF_NLSOLVERS_H_INCLUDED
