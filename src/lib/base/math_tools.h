#ifndef _MATH_TOOLS_H_
#define _MATH_TOOLS_H_

//##################################################################################
//
// FEM Base
//
// Copyright (c) 2016 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   math_tools.h
 * @brief  FEMBase Header
 * @author T. Otani
 */

#include "allocation.h"

class mathTool {
 public:
    mathTool(){};
    ~mathTool(){};
    static double rnd();
    static double vectorNorm(const int &nump,ARRAY1D<double> &x);
    static double innerProduct(const int &nump,ARRAY1D<double> &x,ARRAY1D<double> &y);
    static void crossProduct(const double (&a)[3],const double (&b)[3],double (&c)[3],double &dc);
    static void calcInverseMatrix_3x3(double (&inv_a)[3][3],const double (&a)[3][3]);
    static double calcDeterminant_3x3(const double (&a)[3][3]);
    static void calcMatrix_x_matrix4(double (&ans)[4][4],const double (&a)[4][4],const double (&b)[4][4]);
    static int Jacobi3x3(const int &ct, const double &eps,const double (&A)[3][3], double (&A1)[3][3],double (&X1)[3][3]);
    static void skewSymmetricTensor(double (&M)[3][3],const double (&v)[3]);
    static void quaternion2rotation(double (&R)[3][3],const double (&q)[4]);
    static void calcRotationMatrix(double (&R)[3][3],const double (&rotAxis)[3],const double angle);
    static void calc_thetaFromRotationMatrix(double (&ql)[3],const double (&R)[3][3]);
};



#endif //_MATH_TOOLS_H_
