/**
 * @file triangleSurfaceCurvature.cpp
 * @brief fem class
 * @author t. otani
 */

#include "triangleSurfaceCurvature.h"

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

void triangleSurfaceCurvature::calcSurfaceMeanCurvature()
{

}

void triangleSurfaceCurvature::exportSurfaceNormal(DOUBLEARRAY2D &normal_export)
{
  normal_export.allocate(numOfNode,3);

  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) normal_export(i,j) = normal(i,j);
  }
}

void triangleSurfaceCurvature::exportSurfaceMeanCurvature(DOUBLEARRAY1D &meanCurvature_export)
{
  meanCurvature_export.allocate(numOfNode);

  for(int i=0;i<numOfNode;i++){
    meanCurvature_export(i) = meanCurvature(i);
  }
}

// #################################################################
/**
 * @brief calc boundary conditions
 * @param [in] stress
 */
void triangleSurfaceCurvature::calc_adjacent_nodes()
{
  int tmp;
  bool b1;

  for(int ic=0;ic<numOfNode;ic++){
    for(int i=0,n=neighborElements[ic].size();i<n;i++){
      for(int j=0;j<3;j++){
        tmp = ie(neighborElements[ic][i],j);
        b1 = true;
        for(int k=0,n=neighborNodes[ic].size();k<n;k++){
          if(neighborNodes[ic][k]==tmp){
            b1 = false;
            break;
          }
        }
        if(b1==false) continue;
        neighborNodes[ic].push_back(tmp);
      }
    }
    sort(neighborNodes[ic].begin(),neighborNodes[ic].end());
  }
}

// #################################################################
/**
 * @brief calc boundary conditions
 * @param [in] stress
 */
void triangleSurfaceCurvature::calc_adjacent_elements()
{
  int tmp;
  for(int ic=0;ic<numOfElm;ic++){
    for(int j=0,n=3;j<n;j++){
      tmp = ie(ic,j);
      neighborElements[tmp].push_back(ic);
    }
  }
}
