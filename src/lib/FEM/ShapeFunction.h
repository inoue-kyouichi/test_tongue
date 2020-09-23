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
#include "allocation.h"

class ShapeFunction2D{
 public:

  static void C2D3_N(ARRAY1D<double> &N,const double &L1,const double &L2,const double &L3)
  {
    N(0) = L1;
    N(1) = L2;
    N(2) = L3;
  }
  static void C2D4_N(ARRAY1D<double> &N,const double &g1,const double &g2)
  {
    N(0) = 2.5e-1 * (1e+0-g1) * (1e+0-g2);
    N(1) = 2.5e-1 * (1e+0+g1) * (1e+0-g2);
    N(2) = 2.5e-1 * (1e+0+g1) * (1e+0+g2);
    N(3) = 2.5e-1 * (1e+0-g1) * (1e+0+g2);
  }
  static void C2D6_N(ARRAY1D<double> &N,const double &L1,const double &L2,const double &L3)
  {
    N(0) = L1 * (2e0*L1-1e0);
    N(1) = L2 * (2e0*L2-1e0);
    N(2) = L3 * (2e0*L3-1e0);
    N(3) = 4e0*L1*L2;
    N(4) = 4e0*L2*L3;
    N(5) = 4e0*L1*L3;
  }
  static void C2D8_N(ARRAY1D<double> &N,const double &g1,const double &g2)
  {
  N(0) = 2.5e-1 * (1e+0-g1) * (1e+0-g2) * (-1e+0-g1-g2);
  N(1) = 2.5e-1 * (1e+0+g1) * (1e+0-g2) * (-1e+0+g1-g2);
  N(2) = 2.5e-1 * (1e+0+g1) * (1e+0+g2) * (-1e+0+g1+g2);
  N(3) = 2.5e-1 * (1e+0-g1) * (1e+0+g2) * (-1e+0-g1+g2);
  N(4) = 5e-1   * (1e+0-g1*g1) * (1e+0-g2);
  N(5) = 5e-1   * (1e+0+g1)    * (1e+0-g2*g2);
  N(6) = 5e-1   * (1e+0-g1*g1) * (1e+0+g2);
  N(7) = 5e-1   * (1e+0-g1)    * (1e+0-g2*g2);
  }
  static void C2D9_N(ARRAY1D<double> &N,const double &g1,const double &g2)
  {
    N(0) =  2.5e-1 * (1e+0-g1) * (1e+0-g2) * g1 * g2;
    N(1) = -2.5e-1 * (1e+0+g1) * (1e+0-g2) * g1 * g2;
    N(2) =  2.5e-1 * (1e+0+g1) * (1e+0+g2) * g1 * g2;
    N(3) = -2.5e-1 * (1e+0-g1) * (1e+0+g2) * g1 * g2;
    N(4) = -5e-1   * (1e+0-g1*g1) * (1e+0-g2) * g2;
    N(5) =  5e-1   * (1e+0+g1) * (1e+0-g2*g2) * g1;
    N(6) =  5e-1   * (1e+0-g1*g1) * (1e+0+g2) * g2;
    N(7) = -5e-1   * (1e+0-g1) * (1e+0-g2*g2) * g1;
    N(8) =           (1e+0-g1*g1) * (1e+0-g2*g2);
  }

  static void C2D3_dNdr(ARRAY2D<double> &dNdr,const double &L1,const double &L2,const double &L3)
  {
    dNdr(0,0) = -1e0;
    dNdr(0,1) = -1e0;
    dNdr(1,0) = 1e0;
    dNdr(1,1) = 0e0;
    dNdr(2,0) = 0e0;
    dNdr(2,1) = 1e0;
  }

  static void C2D4_dNdr(ARRAY2D<double> &dNdr,const double &g1,const double &g2)
  {
    dNdr(0,0) = -2.5e-1 * (1e+0-g2);
    dNdr(0,1) = -2.5e-1 * (1e+0-g1);
    dNdr(1,0) =  2.5e-1 * (1e+0-g2);
    dNdr(1,1) = -2.5e-1 * (1e+0+g1);
    dNdr(2,0) =  2.5e-1 * (1e+0+g2);
    dNdr(2,1) =  2.5e-1 * (1e+0+g1);
    dNdr(3,0) = -2.5e-1 * (1e+0+g2);
    dNdr(3,1) =  2.5e-1 * (1e+0-g1);
  }
  static void C2D6_dNdr(ARRAY2D<double> &dNdr,const double &L1,const double &L2,const double &L3)
  {
    dNdr(0,0) = -4e0*L1+1e0;
    dNdr(0,1) = -4e0*L1+1e0;
    dNdr(1,0) = 4e0*L2-1e0;
    dNdr(1,1) = 0e0;
    dNdr(2,0) = 0e0;
    dNdr(2,1) = 4e0*L3-1e0;
    dNdr(3,0) = 4e0*(L1-L2);
    dNdr(3,1) = -4e0*L2;
    dNdr(4,0) = 4e0*L3;
    dNdr(4,1) = 4e0*L2;
    dNdr(5,0) = -4e0*L3;
    dNdr(5,1) = 4e0*(L1-L3);
  }

  static void C2D8_dNdr(ARRAY2D<double> &dNdr,const double &g1,const double &g2)
  {
    dNdr(0,0) = 2.5e-1 * (1e+0-g2) * (2e+0*g1 +      g2);
    dNdr(0,1) = 2.5e-1 * (1e+0-g1) * (     g1 + 2e+0*g2);
    dNdr(1,0) = 2.5e-1 * (1e+0-g2) * (2e+0*g1 -      g2);
    dNdr(1,1) = 2.5e-1 * (1e+0+g1) * (    -g1 + 2e+0*g2);
    dNdr(2,0) = 2.5e-1 * (1e+0+g2) * (2e+0*g1 +      g2);
    dNdr(2,1) = 2.5e-1 * (1e+0+g1) * (     g1 + 2e+0*g2);
    dNdr(3,0) = 2.5e-1 * (1e+0+g2) * (2e+0*g1 -      g2);
    dNdr(3,1) = 2.5e-1 * (1e+0-g1) * (    -g1 + 2e+0*g2);
    dNdr(4,0) = -g1 * (1e+0-g2);
    dNdr(4,1) = -5e-1 * (1e+0-g1*g1);
    dNdr(5,0) = 5e-1 * (1e+0-g2*g2);
    dNdr(5,1) = -g2 * (1e+0+g1);
    dNdr(6,0) = -g1 * (1e+0+g2);
    dNdr(6,1) = 5e-1 * (1e+0-g1*g1);
    dNdr(7,0) = -5e-1 * (1e+0-g2*g2);
    dNdr(7,1) = -g2 * (1e+0-g1);
  }

  static void C2D9_dNdr(ARRAY2D<double> &dNdr,const double &g1,const double &g2)
  {
  dNdr(0,0) = 2.5e-1  * (1e+0-g2) * g2 * (1e+0 - 2e+0*g1);
  dNdr(0,1) = 2.5e-1  * (1e+0-g1) * g1 * (1e+0 - 2e+0*g2);
  dNdr(1,0) = -2.5e-1 * (1e+0-g2) * g2 * (1e+0 + 2e+0*g1);
  dNdr(1,1) = -2.5e-1 * (1e+0+g1) * g1 * (1e+0 - 2e+0*g2);
  dNdr(2,0) = 2.5e-1  * (1e+0+g2) * g2 * (1e+0 + 2e+0*g1);
  dNdr(2,1) = 2.5e-1  * (1e+0+g1) * g1 * (1e+0 + 2e+0*g2);
  dNdr(3,0) = -2.5e-1 * (1e+0+g2) * g2 * (1e+0 - 2e+0*g1);
  dNdr(3,1) = -2.5e-1 * (1e+0-g1) * g1 * (1e+0 + 2e+0*g2);
  dNdr(4,0) = g1 * (1e+0-g2) * g2;
  dNdr(4,1) = -5e-1 * (1e+0-g1*g1) * (1e+0-2e+0*g2);
  dNdr(5,0) = 5e-1 * (1e+0-g2*g2) * (1e+0+2e+0*g1);
  dNdr(5,1) = -g2 * (1e+0+g1) * g1;
  dNdr(6,0) = -g1 * (1e+0+g2) * g2;
  dNdr(6,1) = 5e-1    * (1e+0-g1*g1) * (1e+0+2e+0*g2);
  dNdr(7,0) = -5e-1   * (1e+0-g2*g2) * (1e+0-2e+0*g1);
  dNdr(7,1) = g2 * (1e+0-g1) * g1;
  dNdr(8,0) = -2e+0 * g1 * (1e+0-g2*g2);
  dNdr(8,1) = -2e+0 * g2 * (1e+0-g1*g1);
  }
  static void C2D9_d2Ndr2(ARRAY2D<double> &d2Ndr2,const double &g1,const double &g2)
  {
  d2Ndr2(0,0) = -5e-1   * (1e0-g2) * g2;
  d2Ndr2(0,1) = -5e-1   * (1e0-g1) * g1;
  d2Ndr2(0,2) =  2.5e-1 * (1e0-2e0*g1) * (1e0-2e0*g2);
  d2Ndr2(1,0) = -5e-1   * (1e0-g2) * g2;
  d2Ndr2(1,1) =  5e-1   * (1e0+g1) * g1;
  d2Ndr2(1,2) = -2.5e-1 * (1e0+2e0*g1) * (1e0-2e0*g2);
  d2Ndr2(2,0) =  5e-1   * (1e0+g2) * g2;
  d2Ndr2(2,1) =  5e-1   * (1e0+g1) * g1;
  d2Ndr2(2,2) =  2.5e-1 * (1e0+2e0*g1) * (1e0+2e0*g2);
  d2Ndr2(3,0) =  5e-1   * (1e0+g2) * g2;
  d2Ndr2(3,1) = -5e-1   * (1e0-g1) * g1;
  d2Ndr2(3,2) = -2.5e-1 * (1e0-2e0*g1) * (1e0+2e0*g2);
  d2Ndr2(4,0) = (1e0-g2)    * g2;
  d2Ndr2(4,1) = 1e0-g1*g1;
  d2Ndr2(4,2) = g1 * (1e0-2e0*g2);
  d2Ndr2(5,0) =  (1e0-g2*g2);
  d2Ndr2(5,1) = -1e0 * (1e0+g1) * g1;
  d2Ndr2(5,2) = -g2  * (1e0+2e0*g1);
  d2Ndr2(6,0) = -1e0  *(1e0+g2) * g2;
  d2Ndr2(6,1) =  1e0-g1*g1;
  d2Ndr2(6,2) = -g1 * (1e0 + 2e0*g2);
  d2Ndr2(7,0) =  1e0-g2*g2;
  d2Ndr2(7,1) = (1e0-g1) * g1;
  d2Ndr2(7,2) = (1e0-2e0*g1) * g2;
  d2Ndr2(8,0) = -2e0 * (1e0-g2*g2);
  d2Ndr2(8,1) = -2e0 * (1e0-g1*g1);
  d2Ndr2(8,2) = 4e0 * g1 * g2;
  }
 private:
};

class ShapeFunction3D{
 public:
  static void C3D4_N(ARRAY1D<double> &N,const double &L1,const double &L2,const double &L3,const double &L4)
  {
  N(0)=L1;
  N(1)=L2;
  N(2)=L3;
  N(3)=L4;
  }
  static void C3D6_N(ARRAY1D<double> &N,const double &L1,const double &L2,const double &L3,const double &g1)
  {
  N(0) = 5e-1*L1*(1e0-g1);
  N(1) = 5e-1*L2*(1e0-g1);
  N(2) = 5e-1*L3*(1e0-g1);
  N(3) = 5e-1*L1*(1e0+g1);
  N(4) = 5e-1*L2*(1e0+g1);
  N(5) = 5e-1*L3*(1e0+g1);
  }
  static void C3D8_N(ARRAY1D<double> &N,const double &g1,const double &g2,const double &g3)
  {
  N(0)= 1.25e-1 * (1e0-g1) * (1e0-g2) * (1e0-g3);
  N(1) = 1.25e-1 * (1e0+g1) * (1e0-g2) * (1e0-g3);
  N(2) = 1.25e-1 * (1e0+g1) * (1e0+g2) * (1e0-g3);
  N(3) = 1.25e-1 * (1e0-g1) * (1e0+g2) * (1e0-g3);
  N(4) = 1.25e-1 * (1e0-g1) * (1e0-g2) * (1e0+g3);
  N(5) = 1.25e-1 * (1e0+g1) * (1e0-g2) * (1e0+g3);
  N(6) = 1.25e-1 * (1e0+g1) * (1e0+g2) * (1e0+g3);
  N(7) = 1.25e-1 * (1e0-g1) * (1e0+g2) * (1e0+g3);
  }
  static void C3D10_N(ARRAY1D<double> &N,const double &L1,const double &L2,const double &L3,const double &L4)
  {
  N(0) = L1*(2e0*L1-1e0);
  N(1) = L2*(2e0*L2-1e0);
  N(2) = L3*(2e0*L3-1e0);
  N(3) = L4*(2e0*L4-1e0);
  N(4) = 4e0*L1*L2;
  N(5) = 4e0*L1*L3;
  N(6) = 4e0*L1*L4;
  N(7) = 4e0*L2*L3;
  N(8) = 4e0*L3*L4;
  N(9) = 4e0*L4*L2;
  }
  static void C3D20_N(ARRAY1D<double> &N,const double &g1,const double &g2,const double &g3)
  {
  N(0)  = -1.25e-1 * (1e+0-g1) * (1e+0-g2) * (1e+0-g3) * (2e+0+g1+g2+g3);
  N(1)  = -1.25e-1 * (1e+0+g1) * (1e+0-g2) * (1e+0-g3) * (2e+0-g1+g2+g3);
  N(2)  = -1.25e-1 * (1e+0+g1) * (1e+0+g2) * (1e+0-g3) * (2e+0-g1-g2+g3);
  N(3)  = -1.25e-1 * (1e+0-g1) * (1e+0+g2) * (1e+0-g3) * (2e+0+g1-g2+g3);
  N(4)  = -1.25e-1 * (1e+0-g1) * (1e+0-g2) * (1e+0+g3) * (2e+0+g1+g2-g3);
  N(5)  = -1.25e-1 * (1e+0+g1) * (1e+0-g2) * (1e+0+g3) * (2e+0-g1+g2-g3);
  N(6)  = -1.25e-1 * (1e+0+g1) * (1e+0+g2) * (1e+0+g3) * (2e+0-g1-g2-g3);
  N(7)  = -1.25e-1 * (1e+0-g1) * (1e+0+g2) * (1e+0+g3) * (2e+0+g1-g2-g3);
  N(8)  =  2.5e-1  * (1e+0-g1*g1) * (1e+0-g2) * (1e+0-g3);
  N(9)  =  2.5e-1  * (1e+0+g1) * (1e+0-g2*g2) * (1e+0-g3);
  N(10) =  2.5e-1  * (1e+0-g1*g1) * (1e+0+g2) * (1e+0-g3);
  N(11) =  2.5e-1  * (1e+0-g1) * (1e+0-g2*g2) * (1e+0-g3);
  N(12) =  2.5e-1  * (1e+0-g1*g1) * (1e+0-g2) * (1e+0+g3);
  N(13) =  2.5e-1  * (1e+0+g1) * (1e+0-g2*g2) * (1e+0+g3);
  N(14) =  2.5e-1  * (1e+0-g1*g1) * (1e+0+g2) * (1e+0+g3);
  N(15) =  2.5e-1  * (1e+0-g1) * (1e+0-g2*g2) * (1e+0+g3);
  N(16) =  2.5e-1  * (1e+0-g1) * (1e+0-g2) * (1e+0-g3*g3);
  N(17) =  2.5e-1  * (1e+0+g1) * (1e+0-g2) * (1e+0-g3*g3);
  N(18) =  2.5e-1  * (1e+0+g1) * (1e+0+g2) * (1e+0-g3*g3);
  N(19) =  2.5e-1  * (1e+0-g1) * (1e+0-g2) * (1e+0-g3*g3);
  }

  static void C3D27_N(ARRAY1D<double> &N,const double &g1,const double &g2,const double &g3)
  {
  N(0)   = H1(g1) * H1(g2) * H1(g3);
  N(1)   = H2(g1) * H1(g2) * H1(g3);
  N(2)   = H2(g1) * H2(g2) * H1(g3);
  N(3)   = H1(g1) * H2(g2) * H1(g3);
  N(4)   = H1(g1) * H1(g2) * H2(g3);
  N(5)   = H2(g1) * H1(g2) * H2(g3);
  N(6)   = H2(g1) * H2(g2) * H2(g3);
  N(7)   = H1(g1) * H2(g2) * H2(g3);
  N(8)   = H3(g1) * H1(g2) * H1(g3);
  N(9)   = H2(g1) * H3(g2) * H1(g3);
  N(10)  = H3(g1) * H2(g2) * H1(g3);
  N(11)  = H1(g1) * H3(g2) * H1(g3);
  N(12)  = H3(g1) * H1(g2) * H2(g3);
  N(13)  = H2(g1) * H3(g2) * H2(g3);
  N(14)  = H3(g1) * H2(g2) * H2(g3);
  N(15)  = H1(g1) * H3(g2) * H2(g3);
  N(16)  = H1(g1) * H1(g2) * H3(g3);
  N(17)  = H2(g1) * H1(g2) * H3(g3);
  N(18)  = H2(g1) * H2(g2) * H3(g3);
  N(19)  = H1(g1) * H2(g2) * H3(g3);
  N(20)  = H1(g1) * H3(g2) * H3(g3);
  N(21)  = H2(g1) * H3(g2) * H3(g3);
  N(22)  = H3(g1) * H1(g2) * H3(g3);
  N(23)  = H3(g1) * H2(g2) * H3(g3);
  N(24)  = H3(g1) * H3(g2) * H1(g3);
  N(25)  = H3(g1) * H3(g2) * H2(g3);
  N(26)  = H3(g1) * H3(g2) * H3(g3);
  }

  static void C3D4_dNdr(ARRAY2D<double> &dNdr,const double &L1,const double &L2,const double &L3,const double &L4)
  {
  dNdr(0,0)=-1e0; dNdr(0,1)=-1e0; dNdr(0,2)=-1e0;
  dNdr(1,0)= 1e0; dNdr(1,1)= 0e0; dNdr(1,2)= 0e0;
  dNdr(2,0)= 0e0; dNdr(2,1)= 1e0; dNdr(2,2)= 0e0;
  dNdr(3,0)= 0e0; dNdr(3,1)= 0e0; dNdr(3,2)= 1e0;
  }
  static void C3D6_dNdr(ARRAY2D<double> &dNdr,const double &L1,const double &L2,const double &L3,const double &g1)
  {
  dNdr(0,0) = -5e-1*(1e0-g1);
  dNdr(0,1) = -5e-1*(1e0-g1);
  dNdr(0,2) = -5e-1*L1;

  dNdr(1,0) = 5e-1*(1e0-g1);
  dNdr(1,1) = 0e0;
  dNdr(1,2) = -5e-1*L2;

  dNdr(2,0) = 0e0;
  dNdr(2,1) = 5e-1*(1e0-g1);
  dNdr(2,2) = -5e-1*L3;

  dNdr(3,0) = -5e-1*(1e0+g1);
  dNdr(3,1) = -5e-1*(1e0+g1);
  dNdr(3,2) = 5e-1*L1;

  dNdr(4,0) = 5e-1*(1e0+g1);
  dNdr(4,1) = 0e0;
  dNdr(4,2) = 5e-1*L2;

  dNdr(5,0) = 0e0;
  dNdr(5,1) = 5e-1*(1e0+g1);
  dNdr(5,2) = 5e-1*L3;
  }

  static void C3D8_dNdr(ARRAY2D<double> &dNdr,const double &g1,const double &g2,const double &g3)
  {
  dNdr(0,0) = -1.25e-1 * (1e0-g2) * (1e0-g3);
  dNdr(0,1) = -1.25e-1 * (1e0-g1) * (1e0-g3);
  dNdr(0,2) = -1.25e-1 * (1e0-g1) * (1e0-g2);
  dNdr(1,0) =  1.25e-1 * (1e0-g2) * (1e0-g3);
  dNdr(1,1) = -1.25e-1 * (1e0+g1) * (1e0-g3);
  dNdr(1,2) = -1.25e-1 * (1e0+g1) * (1e0-g2);
  dNdr(2,0) =  1.25e-1 * (1e0+g2) * (1e0-g3);
  dNdr(2,1) =  1.25e-1 * (1e0+g1) * (1e0-g3);
  dNdr(2,2) = -1.25e-1 * (1e0+g1) * (1e0+g2);
  dNdr(3,0) = -1.25e-1 * (1e0+g2) * (1e0-g3);
  dNdr(3,1) =  1.25e-1 * (1e0-g1) * (1e0-g3);
  dNdr(3,2) = -1.25e-1 * (1e0-g1) * (1e0+g2);
  dNdr(4,0) = -1.25e-1 * (1e0-g2) * (1e0+g3);
  dNdr(4,1) = -1.25e-1 * (1e0-g1) * (1e0+g3);
  dNdr(4,2) =  1.25e-1 * (1e0-g1) * (1e0-g2);
  dNdr(5,0) =  1.25e-1 * (1e0-g2) * (1e0+g3);
  dNdr(5,1) = -1.25e-1 * (1e0+g1) * (1e0+g3);
  dNdr(5,2) =  1.25e-1 * (1e0+g1) * (1e0-g2);
  dNdr(6,0) =  1.25e-1 * (1e0+g2) * (1e0+g3);
  dNdr(6,1) =  1.25e-1 * (1e0+g1) * (1e0+g3);
  dNdr(6,2) =  1.25e-1 * (1e0+g1) * (1e0+g2);
  dNdr(7,0) = -1.25e-1 * (1e0+g2) * (1e0+g3);
  dNdr(7,1) =  1.25e-1 * (1e0-g1) * (1e0+g3);
  dNdr(7,2) =  1.25e-1 * (1e0-g1) * (1e0+g2);
  }
  static void C3D10_dNdr(ARRAY2D<double> &dNdr,const double &L1,const double &L2,const double &L3,const double &L4)
  {
  dNdr(0,0) = 4e0*L1-1e0;
  dNdr(0,1) = 0e0;
  dNdr(0,2) = 0e0;

  dNdr(2,0) = 0e0;
  dNdr(2,1) = 4e0*L2-1e0;
  dNdr(2,2) = 0e0;

  dNdr(1,0) = 0e0;
  dNdr(1,1) = 0e0;
  dNdr(1,2) = 4e0*L3-1e0;

  dNdr(3,0) = -4e0*L4+1e0;
  dNdr(3,1) = -4e0*L4+1e0;
  dNdr(3,2) = -4e0*L4+1e0;

  dNdr(6,0) = 4e0*L2;
  dNdr(6,1) = 4e0*L1;
  dNdr(6,2) = 0e0;

  dNdr(5,0) = 0e0;
  dNdr(5,1) = 4e0*L3;
  dNdr(5,2) = 4e0*L2;

  dNdr(4,0) = 4e0*L3;
  dNdr(4,1) = 0e0;
  dNdr(4,2) = 4e0*L1;

  dNdr(7,0) = 4e0*(L4-L1);
  dNdr(7,1) = -4e0*L1;
  dNdr(7,2) = -4e0*L1;

  dNdr(9,0) = -4e0*L2;
  dNdr(9,1) = 4e0*(L4-L2);
  dNdr(9,2) = -4e0*L2;

  dNdr(8,0) = -4e0*L3;
  dNdr(8,1) = -4e0*L3;
  dNdr(8,2) =  4e0*(L4-L3);

  }

  static void C3D20_dNdr(ARRAY2D<double> &dNdr,const double &g1,const double &g2,const double &g3)
  {
  dNdr(0,0) =  1.25e-1 * (1e+0-g2) * (1e+0-g3) * (1e+0 + 2e+0*g1 + g2      + g3);
  dNdr(0,1) =  1.25e-1 * (1e+0-g1) * (1e+0-g3) * (1e+0 + g1      + 2e+0*g2 + g3);
  dNdr(0,2) =  1.25e-1 * (1e+0-g1) * (1e+0-g2) * (1e+0 + g1      + g2      + 2e+0*g3);

  dNdr(1,0) = -1.25e-1 * (1e+0-g2) * (1e+0-g3) * (1e+0 - 2e+0*g1 + g2      + g3);
  dNdr(1,1) =  1.25e-1 * (1e+0+g1) * (1e+0-g3) * (1e+0 - g1      + 2e+0*g2 + g3);
  dNdr(1,2) =  1.25e-1 * (1e+0+g1) * (1e+0-g2) * (1e+0 - g1      + g2      + 2e+0*g3);

  dNdr(2,0) = -1.25e-1 * (1e+0+g2) * (1e+0-g3) * (1e+0 - 2e+0*g1 - g2      + g3);
  dNdr(2,1) = -1.25e-1 * (1e+0+g1) * (1e+0-g3) * (1e+0 - g1      - 2e+0*g2 + g3);
  dNdr(2,2) =  1.25e-1 * (1e+0+g1) * (1e+0+g2) * (1e+0 - g1      - g2      + 2e+0*g3);

  dNdr(3,0) =  1.25e-1 * (1e+0+g2) * (1e+0-g3) * (1e+0 + 2e+0*g1 - g2      + g3);
  dNdr(3,1) = -1.25e-1 * (1e+0-g1) * (1e+0-g3) * (1e+0 + g1      - 2e+0*g2 + g3);
  dNdr(3,2) =  1.25e-1 * (1e+0-g1) * (1e+0+g2) * (1e+0 + g1      - g2      + 2e+0*g3);

  dNdr(4,0) =  1.25e-1 * (1e+0-g2) * (1e+0+g3) * (1e+0 + 2e+0*g1 + g2      - g3);
  dNdr(4,1) =  1.25e-1 * (1e+0-g1) * (1e+0+g3) * (1e+0 + g1      + 2e+0*g2 - g3);
  dNdr(4,2) = -1.25e-1 * (1e+0-g1) * (1e+0-g2) * (1e+0 + g1      + g2      - 2e+0*g3);

  dNdr(5,0) = -1.25e-1 * (1e+0-g2) * (1e+0+g3) * (1e+0 - 2e+0*g1+g2-g3);
  dNdr(5,1) =  1.25e-1 * (1e+0+g1) * (1e+0+g3) * (1e+0 - g1+2e+0*g2-g3);
  dNdr(5,2) = -1.25e-1 * (1e+0+g1) * (1e+0-g2) * (1e+0 - g1+g2-2e+0*g3);

  dNdr(6,0) = -1.25e-1 * (1e+0+g2) * (1e+0+g3) * (1e+0 - 2e+0*g1-g2-g3);
  dNdr(6,1) = -1.25e-1 * (1e+0+g1) * (1e+0+g3) * (1e+0 - g1-2e+0*g2-g3);
  dNdr(6,2) = -1.25e-1 * (1e+0+g1) * (1e+0+g2) * (1e+0 - g1-g2-2e+0*g3);

  dNdr(7,0) =  1.25e-1 * (1e+0+g2) * (1e+0+g3) * (1e+0 + 2e+0*g1-g2-g3);
  dNdr(7,1) = -1.25e-1 * (1e+0-g1) * (1e+0+g3) * (1e+0 + g1-2e+0*g2-g3);
  dNdr(7,2) = -1.25e-1 * (1e+0-g1) * (1e+0+g2) * (1e+0 + g1-g2-2e+0*g3);

  dNdr(8,0) =    -5e-1 * g1 * (1e+0-g2) * (1e+0-g3);
  dNdr(8,1) =  -2.5e-1 * (1e+0-g1*g1)   * (1e+0-g3);
  dNdr(8,2) =  -2.5e-1 * (1e+0-g1*g1)   * (1e+0-g2);

  dNdr(9,0) =   2.5e-1 * (1e+0-g2*g2)   * (1e+0-g3);
  dNdr(9,1) =    -5e-1 * g2 * (1e+0+g1) * (1e+0-g3);
  dNdr(9,2) =  -2.5e-1 * (1e+0+g1)   * (1e+0-g2*g2);

  dNdr(10,0) =   -5e-1 * g1 * (1e+0+g2) * (1e+0-g3);
  dNdr(10,1) =  2.5e-1 * (1e+0-g1*g1)   * (1e+0-g3);
  dNdr(10,2) = -2.5e-1 * (1e+0-g1*g1)   * (1e+0+g2);

  dNdr(11,0) = -2.5e-1 * (1e+0-g2*g2)   * (1e+0-g3);
  dNdr(11,1) =   -5e-1 * g2 * (1e+0-g1) * (1e+0-g3);
  dNdr(11,2) = -2.5e-1 * (1e+0-g1)   * (1e+0-g2*g2);

  dNdr(12,0) =   -5e-1 * g1 * (1e+0-g2) * (1e+0+g3);
  dNdr(12,1) = -2.5e-1 * (1e+0-g1*g1)   * (1e+0+g3);
  dNdr(12,2) =  2.5e-1 * (1e+0-g1*g1)   * (1e+0-g2);

  dNdr(13,0) =  2.5e-1 * (1e+0-g2*g2)   * (1e+0+g3);
  dNdr(13,1) =   -5e-1 * g2 * (1e+0+g1) * (1e+0+g3);
  dNdr(13,2) =  2.5e-1 * (1e+0+g1)   * (1e+0-g2*g2);

  dNdr(14,0) =   -5e-1 * g1 * (1e+0+g2) * (1e+0+g3);
  dNdr(14,1) =  2.5e-1 * (1e+0-g1*g1)   * (1e+0+g3);
  dNdr(14,2) =  2.5e-1 * (1e+0-g1*g1)   * (1e+0+g2);

  dNdr(15,0) = -2.5e-1 * (1e+0-g2*g2)   * (1e+0+g3);
  dNdr(15,1) =   -5e-1 * g2 * (1e+0-g1) * (1e+0+g3);
  dNdr(15,2) =  2.5e-1 * (1e+0-g1)   * (1e+0-g2*g2);

  dNdr(16,0) = -2.5e-1 * (1e+0-g2)   * (1e+0-g3*g3);
  dNdr(16,1) = -2.5e-1 * (1e+0-g1)   * (1e+0-g3*g3);
  dNdr(16,2) =   -5e-1 * g3 * (1e+0-g1) * (1e+0-g2);

  dNdr(17,0) =  2.5e-1 * (1e+0-g2)   * (1e+0-g3*g3);
  dNdr(17,1) = -2.5e-1 * (1e+0+g1)   * (1e+0-g3*g3);
  dNdr(17,2) =   -5e-1 * g3 * (1e+0+g1) * (1e+0-g2);

  dNdr(18,0) =  2.5e-1 * (1e+0+g2)   * (1e+0-g3*g3);
  dNdr(18,1) =  2.5e-1 * (1e+0+g1)   * (1e+0-g3*g3);
  dNdr(18,2) =   -5e-1 * g3 * (1e+0+g1) * (1e+0+g2);

  dNdr(19,0) = -2.5e-1 * (1e+0+g2)   * (1e+0-g3*g3);
  dNdr(19,1) =  2.5e-1 * (1e+0-g1)   * (1e+0-g3*g3);
  dNdr(19,2) =   -5e-1 * g3 * (1e+0-g1) * (1e+0+g2);

  }
  static void C3D27_dNdr(ARRAY2D<double> &dNdr,const double &g1,const double &g2,const double &g3)
  {
  dNdr(0,0) = dH1(g1) * H1(g2)  * H1(g3);
  dNdr(0,1) = H1(g1)  * dH1(g2) * H1(g3);
  dNdr(0,2) = H1(g1)  * H1(g2)  * dH1(g3);

  dNdr(1,0) = dH2(g1) * H1(g2)  * H1(g3);
  dNdr(1,1) = H2(g1)  * dH1(g2) * H1(g3);
  dNdr(1,2) = H2(g1)  * H1(g2)  * dH1(g3);

  dNdr(2,0) = dH2(g1) * H2(g2)  * H1(g3);
  dNdr(2,1) = H2(g1)  * dH2(g2) * H1(g3);
  dNdr(2,2) = H2(g1)  * H2(g2)  * dH1(g3);

  dNdr(3,0) = dH1(g1) * H2(g2)  * H1(g3);
  dNdr(3,1) = H1(g1)  * dH2(g2) * H1(g3);
  dNdr(3,2) = H1(g1)  * H2(g2)  * dH1(g3);

  dNdr(4,0) = dH1(g1) * H1(g2)  * H2(g3);
  dNdr(4,1) = H1(g1)  * dH1(g2) * H2(g3);
  dNdr(4,2) = H1(g1)  * H1(g2)  * dH2(g3);

  dNdr(5,0) = dH2(g1) * H1(g2)  * H2(g3);
  dNdr(5,1) = H2(g1)  * dH1(g2) * H2(g3);
  dNdr(5,2) = H2(g1)  * H1(g2)  * dH2(g3);

  dNdr(6,0) = dH2(g1) * H2(g2)  * H2(g3);
  dNdr(6,1) = H2(g1)  * dH2(g2) * H2(g3);
  dNdr(6,2) = H2(g1)  * H2(g2)  * dH2(g3);

  dNdr(7,0) = dH1(g1) * H2(g2)  * H2(g3);
  dNdr(7,1) = H1(g1)  * dH2(g2) * H2(g3);
  dNdr(7,2) = H1(g1)  * H2(g2)  * dH2(g3);

  dNdr(8,0) = dH3(g1) * H1(g2)  * H1(g3);
  dNdr(8,1) = H3(g1)  * dH1(g2) * H1(g3);
  dNdr(8,2) = H3(g1)  * H1(g2)  * dH1(g3);

  dNdr(9,0) = dH2(g1) * H3(g2)  * H1(g3);
  dNdr(9,1) = H2(g1)  * dH3(g2) * H1(g3);
  dNdr(9,2) = H2(g1)  * H3(g2)  * dH1(g3);

  dNdr(10,0)  = dH3(g1) * H2(g2)  * H1(g3);
  dNdr(10,1)  = H3(g1)  * dH2(g2) * H1(g3);
  dNdr(10,2)  = H3(g1)  * H2(g2)  * dH1(g3);

  dNdr(11,0)  = dH1(g1) * H3(g2)  * H1(g3);
  dNdr(11,1)  = H1(g1)  * dH3(g2) * H1(g3);
  dNdr(11,2)  = H1(g1)  * H3(g2)  * dH1(g3);

  dNdr(12,0)  = dH3(g1) * H1(g2)  * H2(g3);
  dNdr(12,1)  = H3(g1)  * dH1(g2) * H2(g3);
  dNdr(12,2)  = H3(g1)  * H1(g2)  * dH2(g3);

  dNdr(13,0)  = dH2(g1) * H3(g2)  * H2(g3);
  dNdr(13,1)  = H2(g1)  * dH3(g2) * H2(g3);
  dNdr(13,2)  = H2(g1)  * H3(g2)  * dH2(g3);

  dNdr(14,0)  = dH3(g1) * H2(g2)  * H2(g3);
  dNdr(14,1)  = H3(g1)  * dH2(g2) * H2(g3);
  dNdr(14,2)  = H3(g1)  * H2(g2)  * dH2(g3);

  dNdr(15,0)  = dH1(g1) * H3(g2)  * H2(g3);
  dNdr(15,1)  = H1(g1)  * dH3(g2) * H2(g3);
  dNdr(15,2)  = H1(g1)  * H3(g2)  * dH2(g3);

  dNdr(16,0)  = dH1(g1) * H1(g2)  * H3(g3);
  dNdr(16,1)  = H1(g1)  * dH1(g2) * H3(g3);
  dNdr(16,2)  = H1(g1)  * H1(g2)  * dH3(g3);

  dNdr(17,0)  = dH2(g1) * H1(g2)  * H3(g3);
  dNdr(17,1)  = H2(g1)  * dH1(g2) * H3(g3);
  dNdr(17,2)  = H2(g1)  * H1(g2)  * dH3(g3);

  dNdr(18,0)  = dH2(g1) * H2(g2)  * H3(g3);
  dNdr(18,1)  = H2(g1)  * dH2(g2) * H3(g3);
  dNdr(18,2)  = H2(g1)  * H2(g2)  * dH3(g3);

  dNdr(19,0)  = dH1(g1) * H2(g2)  * H3(g3);
  dNdr(19,1)  = H1(g1)  * dH2(g2) * H3(g3);
  dNdr(19,2)  = H1(g1)  * H2(g2)  * dH3(g3);

  dNdr(20,0)  = dH1(g1) * H3(g2) * H3(g3);
  dNdr(20,1)  = H1(g1) * dH3(g2) * H3(g3);
  dNdr(20,2)  = H1(g1) * H3(g2) * dH3(g3);

  dNdr(21,0)  = dH2(g1) * H3(g2) * H3(g3);
  dNdr(21,1)  = H2(g1) * dH3(g2) * H3(g3);
  dNdr(21,2)  = H2(g1) * H3(g2) * dH3(g3);

  dNdr(22,0)  = dH3(g1) * H1(g2) * H3(g3);
  dNdr(22,1)  = H3(g1) * dH1(g2) * H3(g3);
  dNdr(22,2)  = H3(g1) * H1(g2) * dH3(g3);

  dNdr(23,0)  = dH3(g1) * H2(g2) * H3(g3);
  dNdr(23,1)  = H3(g1) * dH2(g2) * H3(g3);
  dNdr(23,2)  = H3(g1) * H2(g2) * dH3(g3);

  dNdr(24,0)  = dH3(g1) * H3(g2) * H1(g3);
  dNdr(24,1)  = H3(g1) * dH3(g2) * H1(g3);
  dNdr(24,2)  = H3(g1) * H3(g2) * dH1(g3);

  dNdr(25,0)  = dH3(g1) * H3(g2) * H2(g3);
  dNdr(25,1)  = H3(g1) * dH3(g2) * H2(g3);
  dNdr(25,2)  = H3(g1) * H3(g2) * dH2(g3);

  dNdr(26,0)  = dH3(g1) * H3(g2) * H3(g3);
  dNdr(26,1)  = H3(g1) * dH3(g2) * H3(g3);
  dNdr(26,2)  = H3(g1) * H3(g2) * dH3(g3);
  }

 private:
  static double dH1(const double r) { return -5e-1 * (1e+0-2e+0*r); }
  static double dH2(const double r) { return 5e-1 * (1e+0+2e+0*r); }
  static double dH3(const double r) { return -2e+0*r; }
  static double H1(const double r) { return -5e-1 * r * (1e+0-r); }
  static double H2(const double r) { return 5e-1 * r * (1e+0+r); }
  static double H3(const double r) { return (1e+0-r*r); }
};

#endif //_SHAPE_FUNCTION_H_