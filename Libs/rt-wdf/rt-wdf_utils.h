//
// Created by Sponge on 2024/11/14.
//

#ifndef RT_WDF_UTILS_H
#define RT_WDF_UTILS_H

using FloatType = double;

// armadillo lib
#include <armadillo>

using Wvec = arma::vec;
using Wmat = arma::mat;
using WcolVec = arma::colvec;

#define VEC_RESIZE(VEC, NELEM) (VEC).set_size(NELEM)
#define MAT_RESIZE(MAT, NROW, NCOL) (MAT).set_size(NROW, NCOL)
#define MAT_SETZERO(MAT) (MAT).zeros()

#define MAT_ALLOC(MAT, NROW, NCOL) MAT = new Wmat( NROW, NCOL, arma::fill::zeros)
#define VEC_ALLOC(VEC, NELEM) VEC = new Wvec( NELEM, arma::fill::zeros)

#define MAT_ASSIGN_ID_BY_MAT(MAT, SIZED_MAT) MAT = arma::eye<Wmat>(size(SIZED_MAT))
#define MAT_CREATE_ZERO_MAT(NROW, NCOL) arma::zeros<Wmat>(NROW, NCOL)

#define W_NORM(OBJ) arma::norm(OBJ)
#define W_INV(MAT) arma::inv(MAT)

#define MAT_DIAG_MAT(VEC) arma::diagmat(VEC)
#define MAT_JOIN_HORIZON(L, R) arma::join_horiz(L, R)
#define MAT_JOIN_VERTICAL(L, R) arma::join_vert(L, R)
#define MAT_CREATE_ID_BY_SIZE(NROW, NCOL) arma::eye<Wmat>(NROW, NCOL)
#define W_SUBVEC(VEC, START, END) (VEC).subvec(START, (END - START + 1))

// eigen lib
//#include "Eigen/Dense"
//
//using Wvec = Eigen::VectorXd;
//using Wmat = Eigen::MatrixXd;
//
//#define VEC_RESIZE(VEC, NELEM) (VEC).resize(NELEM)
//#define MAT_RESIZE(MAT, NROW, NCOL) (MAT).resize(NROW, NCOL)
//#define MAT_SETZERO(MAT) (MAT).setZero()
//
//#define MAT_ALLOC(MAT, NROW, NCOL) MAT = new Wvec( NROW, NCOL); MAT_SETZERO(MAT)
//#define VEC_ALLOC(VEC, NELEM) VEC = new Wvec( NELEM ); MAT_SETZERO(VEC)
//
//#define MAT_ASSIGN_ID_BY_MAT(MAT, SIZED_MAT) MAT = Wmat::Identity((SIZED_MAT).rows(), (SIZED_MAT).cols())
//#define MAT_CREATE_ZERO_MAT(NROW, NCOL) Eigen::MatrixXd::Zero(NROW, NCOL)
//
//#define W_NORM(OBJ) (OBJ).norm()
//#define W_INV(MAT) (MAT).inverse()
//
//#define MAT_DIAG_MAT(VEC) (VEC).asDiagonal()
//#define MAT_JOIN_HORIZON(L, R) Eigen::hstack(L, R)
//#define MAT_JOIN_VERTICAL(U, B) Eigen::hstack(U, B)
//#define MAT_CREATE_ID_BY_SIZE(NROW, NCOL) Wmat::Identity(NROW, NCOL)
//
//#define W_SUBVEC(VEC, START, END)

#endif //RT_WDF_UTILS_H
