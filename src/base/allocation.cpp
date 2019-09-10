/**
 * @file allocation.cpp
 * @brief BDIM class
 * @author T. Otani
 */
#include "allocation.h"
#include <iostream>
#include <cstdlib>
using namespace std;

#ifndef BOOST

// #################################################################
/**
 * @brief allocate 1d integr array
 * @param [in]     size          array size
 * @param [in]     guideCell     guide cell size
 */
  INTARRAY1 Allocation::allocate1dINT(const int &size,const int &guideCell)
  {
    INTARRAY1 var = (int *)malloc((size+2*guideCell) * sizeof(int));
    return var;
  }

// #################################################################
/**
 * @brief allocate 2d integr array [size1xsize2]
 * @param [in]     size1          array size
 * @param [in]     size2          array size
 * @param [in]     guideCell     guide cell size
 */
  INTARRAY2 Allocation::allocate2dINT(const int &size1,const int &size2,const int &guideCell)
  {
    int num1 = size1+2*guideCell;
    int num2 = size2+2*guideCell;

    INTARRAY2 var = (int**)malloc(num1*sizeof(int*));
    var[0] = (int*)malloc(num1*num2*sizeof(int));
    for(int i=1;i<num1;i++) var[i] = var[i-1] + num2;

    return var;
  }

// #################################################################
/**
 * @brief allocate 3d integr array [size1xsize2xsize3]
 * @param [in]     size1          array size
 * @param [in]     size2          array size
 * @param [in]     size3          array size
 * @param [in]     guideCell     guide cell size
 */
  INTARRAY3 Allocation::allocate3dINT(const int &size1,const int &size2,const int &size3,const int &guideCell)
  {
    int num1 = size1+2*guideCell;
    int num2 = size2+2*guideCell;
    int num3 = size3+2*guideCell;

    INTARRAY3 var = (int***)malloc(num1*sizeof(int**));
    var[0]    = (int**)malloc(num1*num2*sizeof(int*));
    var[0][0] = (int*)malloc(num1*num2*num3*sizeof(int));
    for(int i=0;i<num1;i++){
      var[i] = var[0] + i*num2;
      for(int j=0;j<num2;j++) var[i][j] = var[0][0] + i*num2*num3 + j*num3;
    }
    return var;
  }

// #################################################################
/**
 * @brief allocate 4d integr array [size1xsize2xsize3xsize4]
 * @param [in]     size1          array size
 * @param [in]     size2          array size
 * @param [in]     size3          array size
 * @param [in]     size4          array size
 * @param [in]     guideCell     guide cell size
 */
  INTARRAY4 Allocation::allocate4dINT(const int &size1,const int &size2,const int &size3,const int &size4,const int &guideCell)
  {
    int num1 = size1+2*guideCell;
    int num2 = size2+2*guideCell;
    int num3 = size3+2*guideCell;
    int num4 = size4+2*guideCell;

    INTARRAY4 var = (int****)malloc(num1*sizeof(int***));
    var[0]       = (int***)malloc(num1*num2*sizeof(int**));
    var[0][0]    = (int**)malloc(num1*num2*num3*sizeof(int*));
    var[0][0][0] = (int*)malloc(num1*num2*num3*num4*sizeof(int));

    for(int i=0;i<num1;i++){
      var[i] = var[0] + i*num2;
      for(int j=0;j<num2;j++){
        var[i][j] = var[0][0] + i*num2*num3 + j*num3;
        for(int k=0;k<num3;k++) var[i][j][k] = var[0][0][0] + i*num2*num3*num4 + j*num3*num4 + k*num4;
      }
    }
    return var;
  }

// #################################################################
/**
 * @brief allocate 5d integr array [size1xsize2xsize3xsize4xsize5]
 * @param [in]     size1          array size
 * @param [in]     size2          array size
 * @param [in]     size3          array size
 * @param [in]     size4          array size
 * @param [in]     size5          array size
 * @param [in]     guideCell     guide cell size
 */
  INTARRAY5 Allocation::allocate5dINT(const int &size1,const int &size2,const int &size3,const int &size4,const int &size5,const int &guideCell)
  {
    int num1 = size1+2*guideCell;
    int num2 = size2+2*guideCell;
    int num3 = size3+2*guideCell;
    int num4 = size4+2*guideCell;
    int num5 = size5+2*guideCell;

    INTARRAY5 var   = (int*****)malloc(num1*sizeof(int****));
    var[0]          = (int****)malloc(num1*num2*sizeof(int***));
    var[0][0]       = (int***)malloc(num1*num2*num3*sizeof(int**));
    var[0][0][0]    = (int**)malloc(num1*num2*num3*num4*sizeof(int*));
    var[0][0][0][0] = (int*)malloc(num1*num2*num3*num4*num5*sizeof(int));

    for(int i=0;i<num1;i++){
      var[i] = var[0] + i*num2;
      for(int j=0;j<num2;j++){
        var[i][j] = var[0][0] + i*num2*num3 + j*num3;
        for(int k=0;k<num3;k++){
          var[i][j][k] = var[0][0][0] + i*num2*num3*num4 + j*num3*num4 + k*num4;
          for(int l=0;l<num4;l++) var[i][j][k][l] = var[0][0][0][0] + i*num2*num3*num4*num5 + j*num3*num4*num5 + k*num4*num5 + l*num5;
        }
      }
    }
    return var;
  }

// #################################################################
/**
 * @brief allocate 1d double array [size1]
 * @param [in]     size          array size
 * @param [in]     guideCell     guide cell size
 */
  DOUBLEARRAY1 Allocation::allocate1dDOUBLE(const int &size,const int &guideCell)
  {
    DOUBLEARRAY1 var = (double *)malloc((size+2*guideCell) * sizeof(double));
    return var;
  }

// #################################################################
/**
 * @brief allocate 2d double array [size1xsize2]
 * @param [in]     size1          array size
 * @param [in]     size2          array size
 * @param [in]     guideCell     guide cell size
 */
  DOUBLEARRAY2 Allocation::allocate2dDOUBLE(const int &size1,const int &size2,const int &guideCell)
  {
    int num1 = size1+2*guideCell;
    int num2 = size2+2*guideCell;

    DOUBLEARRAY2 var  = (double**)malloc(num1*sizeof(double*));
    var[0] = (double*)malloc(num1*num2*sizeof(double));
    for(int i=1;i<num1;i++) var[i] = var[i-1] + num2;
    return var;
  }

// #################################################################
/**
 * @brief allocate 3d double array [size1xsize2xsize3]
 * @param [in]     size1          array size
 * @param [in]     size2          array size
 * @param [in]     size3          array size
 * @param [in]     guideCell     guide cell size
 */
  DOUBLEARRAY3 Allocation::allocate3dDOUBLE(const int &size1,const int &size2,const int &size3,const int &guideCell)
  {
    int num1 = size1+2*guideCell;
    int num2 = size2+2*guideCell;
    int num3 = size3+2*guideCell;

    DOUBLEARRAY3 var = (double***)malloc(num1*sizeof(double**));
    var[0]    = (double**)malloc(num1*num2*sizeof(double*));
    var[0][0] = (double*)malloc(num1*num2*num3*sizeof(double));
    for(int i=0;i<num1;i++){
      var[i] = var[0] + i*num2;
      for(int j=0;j<num2;j++) var[i][j] = var[0][0] + i*num2*num3 + j*num3;
    }
    return var;
  }

// #################################################################
/**
 * @brief allocate 4d double array [size1xsize2xsize3xsize4]
 * @param [in]     size1          array size
 * @param [in]     size2          array size
 * @param [in]     size3          array size
 * @param [in]     size4          array size
 * @param [in]     guideCell     guide cell size
 */
  DOUBLEARRAY4 Allocation::allocate4dDOUBLE(const int &size1,const int &size2,const int &size3,const int &size4,const int &guideCell)
  {
    int num1 = size1+2*guideCell;
    int num2 = size2+2*guideCell;
    int num3 = size3+2*guideCell;
    int num4 = size4+2*guideCell;

    DOUBLEARRAY4 var = (double****)malloc(num1*sizeof(double***));
    var[0]       = (double***)malloc(num1*num2*sizeof(double**));
    var[0][0]    = (double**)malloc(num1*num2*num3*sizeof(double*));
    var[0][0][0] = (double*)malloc(num1*num2*num3*num4*sizeof(double));

    for(int i=0;i<num1;i++){
      var[i] = var[0] + i*num2;
      for(int j=0;j<num2;j++){
        var[i][j] = var[0][0] + i*num2*num3 + j*num3;
        for(int k=0;k<num3;k++) var[i][j][k] = var[0][0][0] + i*num2*num3*num4 + j*num3*num4 + k*num4;
      }
    }
    return var;
  }

// #################################################################
/**
 * @brief allocate 4d double array [size1xsize2xsize3xsize4xsize5]
 * @param [in]     size1          array size
 * @param [in]     size2          array size
 * @param [in]     size3          array size
 * @param [in]     size4          array size
 * @param [in]     size5          array size
 * @param [in]     guideCell     guide cell size
 */
  DOUBLEARRAY5 Allocation::allocate5dDOUBLE(const int &size1,const int &size2,const int &size3,const int &size4,const int &size5,const int &guideCell)
  {
    int num1 = size1+2*guideCell;
    int num2 = size2+2*guideCell;
    int num3 = size3+2*guideCell;
    int num4 = size4+2*guideCell;
    int num5 = size5+2*guideCell;

    DOUBLEARRAY5 var = (double*****)malloc(num1*sizeof(double****));
    var[0]           = (double****)malloc(num1*num2*sizeof(double***));
    var[0][0]        = (double***)malloc(num1*num2*num3*sizeof(double**));
    var[0][0][0]     = (double**)malloc(num1*num2*num3*num4*sizeof(double*));
    var[0][0][0][0]  = (double*)malloc(num1*num2*num3*num4*num5*sizeof(double));

    for(int i=0;i<num1;i++){
      var[i] = var[0] + i*num2;
      for(int j=0;j<num2;j++){
        var[i][j] = var[0][0] + i*num2*num3 + j*num3;
        for(int k=0;k<num3;k++){
          var[i][j][k] = var[0][0][0] + i*num2*num3*num4 + j*num3*num4 + k*num4;
          for(int l=0;l<num4;l++) var[i][j][k][l] = var[0][0][0][0] + i*num2*num3*num4*num5 + j*num3*num4*num5 + k*num4*num5 + l*num5;
        }
      }
    }
    return var;
  }
#endif