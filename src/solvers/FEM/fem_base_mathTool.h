#ifndef _FEM_BASE_MATH_TOOL_H_
#define _FEM_BASE_MATH_TOOL_H_


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
 * @file   fem_base_mathTooL.h
 * @brief  FEMBase Definition Header
 * @author T. Otani
 */

#include "allocation.h"

class FEM_MathTool{
  public:
  static void calc_dudX(double (&dudX)[3][3],DOUBLEARRAY2D &dNdX,DOUBLEARRAY2D &u,const int &numOfNodeInElm);
  static void calc_dXdr(double (&dXdr)[3][3],DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &x0,const int &numOfNodeInElm);
  static void calc_dxdr(double (&dxdr)[3][3],DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &x,const int &numOfNodeInElm);
  static void calc_dNdX(DOUBLEARRAY2D &dNdX,DOUBLEARRAY2D &dNdr,const double (&dXdr)[3][3],const int &numOfNodeInElm);
  static void calc_dNdx(DOUBLEARRAY2D &dNdx,DOUBLEARRAY2D &dNdr,const double (&dxdr)[3][3],const int &numOfNodeInElm);
  static void tensorPushForward_4order(double (&c4)[3][3][3][3],const double (&C4)[3][3][3][3],const double (&F)[3][3],const double J);

};

#endif