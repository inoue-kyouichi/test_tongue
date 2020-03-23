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

#include <vector>
#include "allocation.h"
#include "math_tools.h"

class triangleSurfaceCurvature{

  //fem.cpp
 public:
  triangleSurfaceCurvature(DOUBLEARRAY2D &x_input,INTARRAY2D &ie_input,const int numOfNode_input,const int numOfElm_input){
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
  }

  int numOfNode,numOfElm;
  DOUBLEARRAY2D x;
  INTARRAY2D ie;
  DOUBLEARRAY2D normal;
  DOUBLEARRAY1D meanCurvature;

  void calcSurfaceMeanCurvature();
  void exportSurfaceNormal(DOUBLEARRAY2D &normal_export);
  void exportSurfaceMeanCurvature(DOUBLEARRAY1D &meanCurvature_export);
 private:

  std::vector<std::vector<int>> neighborElements,neighborNodes;
  void calc_adjacent_nodes();
  void calc_adjacent_elements();
};


#endif
