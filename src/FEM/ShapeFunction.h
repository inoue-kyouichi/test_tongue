#ifndef _SHAPE_FUNCTION_H_
#define _SHAPE_FUNCTION_H_


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
 * @file   shape_function.h
 * @brief  FEMBase Definition Header
 * @author T. Otani
 */

#include "fem_define.h"

class ShapeFunction{
 public:
  static void C2D4_N(DOUBLEARRAY1 &N,const double &g1,const double &g2);
  static void C2D6_N(DOUBLEARRAY1 &N,const double &L1,const double &L2,const double &L3);
  static void C2D8_N(DOUBLEARRAY1 &N,const double &g1,const double &g2);
  static void C2D9_N(DOUBLEARRAY1 &N,const double &g1,const double &g2);
  static void C3D4_N(DOUBLEARRAY1 &N,const double &L1,const double &L2,const double &L3,const double &L4);
  static void C3D8_N(DOUBLEARRAY1 &N,const double &g1,const double &g2,const double &g3);
  static void C3D20_N(DOUBLEARRAY1 &N,const double &g1,const double &g2,const double &g3);
  static void C3D27_N(DOUBLEARRAY1 &N,const double &g1,const double &g2,const double &g3);
  static void C3D10_N(DOUBLEARRAY1 &N,const double &L1,const double &L2,const double &L3,const double &L4);
  static void C2D3_N(DOUBLEARRAY1 &N,const double &L1,const double &L2,const double &L3);
  static void C3D6_N(DOUBLEARRAY1 &N,const double &L1,const double &L2,const double &L3,const double &g1);

  static void C2D4_dNdr(DOUBLEARRAY2 &dNdr,const double &g1,const double &g2);
  static void C2D6_dNdr(DOUBLEARRAY2 &dNdr,const double &L1,const double &L2,const double &L3);
  static void C2D8_dNdr(DOUBLEARRAY2 &dNdr,const double &g1,const double &g2);
  static void C2D9_dNdr(DOUBLEARRAY2 &dNdr,const double &g1,const double &g2);
  static void C3D4_dNdr(DOUBLEARRAY2 &dNdr,const double &L1,const double &L2,const double &L3,const double &L4);
  static void C3D8_dNdr(DOUBLEARRAY2 &dNdr,const double &g1,const double &g2,const double &g3);
  static void C3D20_dNdr(DOUBLEARRAY2 &dNdr,const double &g1,const double &g2,const double &g3);
  static void C3D27_dNdr(DOUBLEARRAY2 &dNdr,const double &g1,const double &g2,const double &g3);
  static void C3D10_dNdr(DOUBLEARRAY2 &dNdr,const double &L1,const double &L2,const double &L3,const double &L4);
  static void C2D3_dNdr(DOUBLEARRAY2 &dNdr,const double &L1,const double &L2,const double &L3);
  static void C3D6_dNdr(DOUBLEARRAY2 &dNdr,const double &L1,const double &L2,const double &L3,const double &g1);

  static void C2D9_d2Ndr2(DOUBLEARRAY2 &d2Ndr2,const double &g1,const double &g2);


 private:
  static double dH1(const double r);
  static double dH2(const double r);
  static double dH3(const double r);

  static double H1(const double r);
  static double H2(const double r);
  static double H3(const double r);

};

#endif //_SHAPE_FUNCTION_H_