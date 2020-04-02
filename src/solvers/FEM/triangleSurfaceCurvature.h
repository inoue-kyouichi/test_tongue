#ifndef _TRIANGLE_SURFACE_CURVATURE_H_
#define _TRIANGLE_SURFACE_CURVATURE_H_

//##################################################################################
//
// surface curvature calculation from linear triangle elements
// Meyer M., Desbrun M., Schröder P., Barr A.H. (2003) Discrete Differential-Geometry Operators for Triangulated 2-Manifolds. In: Hege HC., Polthier K. (eds) Visualization and Mathematics III. Mathematics and Visualization. Springer, Berlin, Heidelberg
// Copyright (c) 2020 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   triangleSurfaceCurvature.h
 * @brief  FEMBase Definition Header
 * @author T. Otani
 */

#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <strings.h>
#include <sys/stat.h>
#include <omp.h>
#include <chrono>
#include <algorithm>
#include <vector>

#include "allocation.h"
#include "math_tools.h"

class triangleSurfaceCurvature{

  //fem.cpp
 public:
  triangleSurfaceCurvature(ARRAY2D<double> &x_input,ARRAY2D<int> &ie_input,const int numOfNode_input,const int numOfElm_input){
    numOfNode = numOfNode_input;
    numOfElm = numOfElm_input;

    x.allocate(numOfNode,3);
    ie.allocate(numOfElm,3);

    for(int i=0;i<numOfNode;i++){
      for(int j=0;j<3;j++) x(i,j) = x_input(i,j);
    }

    for(int i=0;i<numOfElm;i++){
      for(int j=0;j<3;j++) ie(i,j) = ie_input(i,j);
    }
    normal.allocate(numOfNode,3);
    meanCurvature.allocate(numOfNode);

    neighborElements.resize(numOfNode);
    neighborNodes.resize(numOfNode);
    calc_adjacent_elements();
    calc_adjacent_nodes();
  }

  int numOfNode,numOfElm;
  ARRAY2D<double> x;
  ARRAY2D<int> ie;
  ARRAY2D<double> normal;
  ARRAY1D<double> meanCurvature;

  void calcSurfaceMeanCurvature();
  void exportSurfaceNormal(ARRAY2D<double> &normal_export);
  void exportSurfaceMeanCurvature(ARRAY1D<double> &meanCurvature_export);
 private:

  void modifyNormalVectorDirection();
  bool acuteTriangle(double (&x1)[3],double (&x2)[3]);
  double cot(const double (&x1)[3],const double (&x2)[3]);
  std::vector<std::vector<int>> neighborElements,neighborNodes;
  void calc_adjacent_nodes();
  void calc_adjacent_elements();
};


#endif
