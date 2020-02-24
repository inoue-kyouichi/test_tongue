#ifndef _PDL_H_
#define _PDL_H_


//##################################################################################
//
// PDL
//
// Copyright (c) 2020 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   PDL.h
 * @brief  PDL Header
 * @author T. Otani
 */

#include "fem.h"

class PeriodontalLigament:public Fem{
 public:
  double bulkModulusRatio;
  DOUBLEARRAY3D fiberDirection;
  DOUBLEARRAY2D lambda_ave;

  DOUBLEARRAY2D fiberDirection_elm;

  void calcStressTensor();
  void allocatePDLvariables();

  void postProcess_PDL_element_2018(const int &ic,DOUBLEARRAY2D &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint);

  void inputFiberInfo(TextParser &tp);

  void setFiberDirection();
  void setFiberDirection_KogaModel(TextParser &tp);

  // void postProcess_PDL_element_spatialForm_hexa_SRI(const int &ic,DOUBLEARRAY2D &U_tmp,
  // const int &numOfNodeInElm,const int &numOfGaussPoint);

 private:
  void calc_normal_quad(double (&normal)[3],double (&X)[4][3]);
  double arcLength(const double xMin,const double xMax);

  // void calcStressTensor_PDL_element_spatialForm_hexa_Fbar(const int &ic,DOUBLEARRAY2D &U_tmp,
  // const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  // void calcStressTensor_PDL_element_spatialForm_hexa_SRI(const int &ic,DOUBLEARRAY2D &U_tmp,
  // const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  // void calcStressTensor_PDL_element_spatialForm_hexa_2018(const int &ic,DOUBLEARRAY2D &U_tmp,
  // const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  // void calcStressTensor_PDL_element_spatialForm_hexa_SRI_2018(const int &ic,DOUBLEARRAY2D &U_tmp,
  // const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  // int calcStressTensor_PDL_element_fibreStretch(const int &ic,DOUBLEARRAY2D &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint);

  void calcStressTensor_PDL_element_2018(const int &ic,DOUBLEARRAY2D &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint);
  void calcStressTensor_PDL_element_spatialForm_hexa_2018_inGaussIntegral(const int &ic,DOUBLEARRAY2D &U_tmp,
      const int &numOfNodeInElm,const Gauss &gauss,DOUBLEARRAY2D &x_current,DOUBLEARRAY2D &x_ref,DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &dNdx,const int i1,const int i2,const int i3,double (&stress)[3][3],const bool mainLoop);
  void calcStressTensor_hyperFoam_element_spatialForm_hexa_inGaussIntegral(const int &ic,DOUBLEARRAY2D &U_tmp,
  const int &numOfNodeInElm,const Gauss &gauss,DOUBLEARRAY2D &x_current,DOUBLEARRAY2D &x_ref,DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &dNdx,const int i1,const int i2,const int i3,double (&stress)[3][3],const bool mainLoop);
  //fem_hyperFoam_spatialForm.cpp
  // void calcStressTensor_hyperFoam_element_spatialForm_hexa(const int &ic,DOUBLEARRAY2D &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  void calcLambda(double (&stretch)[3],double (&stretchDirection)[3][3],const double (&C)[3][3]);
};

#endif