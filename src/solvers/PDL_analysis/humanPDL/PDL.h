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
  ARRAY3D<double> fiberDirection;
  ARRAY2D<double> lambda_ave;

  ARRAY2D<double> fiberDirection_elm;

  void calcStressTensor();
  void allocatePDLvariables();

  void postProcess_PDL_element_2018(const int &ic,ARRAY2D<double> &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint);

  void inputFiberInfo(TextParser &tp);

  void setFiberDirection();
  void setFiberDirection_KogaModel(TextParser &tp);

  // void postProcess_PDL_element_spatialForm_hexa_SRI(const int &ic,ARRAY2D<double> &U_tmp,
  // const int &numOfNodeInElm,const int &numOfGaussPoint);

 private:
  void calc_normal_quad(double (&normal)[3],double (&X)[4][3]);
  double arcLength(const double xMin,const double xMax);

  // void calcStressTensor_PDL_element_spatialForm_hexa_Fbar(const int &ic,ARRAY2D<double> &U_tmp,
  // const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  // void calcStressTensor_PDL_element_spatialForm_hexa_SRI(const int &ic,ARRAY2D<double> &U_tmp,
  // const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  // void calcStressTensor_PDL_element_spatialForm_hexa_2018(const int &ic,ARRAY2D<double> &U_tmp,
  // const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  // void calcStressTensor_PDL_element_spatialForm_hexa_SRI_2018(const int &ic,ARRAY2D<double> &U_tmp,
  // const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  // int calcStressTensor_PDL_element_fibreStretch(const int &ic,ARRAY2D<double> &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint);

  double calcI4bar(double (&F)[3][3],ARRAY2D<double> &dNdr,ARRAY2D<double> &dNdx,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,const int numOfNodeInElm,const int ic);

  void calcStressTensor_PDL_element_2018(const int &ic,ARRAY2D<double> &U_tmp);
  void calcStressTensor_PDL_element_spatialForm_2018_inGaussIntegral(const int &ic,ARRAY2D<double> &U_tmp,
      const int &numOfNodeInElm,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,ARRAY2D<double> &dNdr,ARRAY2D<double> &dNdx,const double weight,double (&stress)[3][3],const bool mainLoop);
  void calcStressTensor_hyperFoam_element_spatialForm_inGaussIntegral(const int &ic,ARRAY2D<double> &U_tmp,
  const int &numOfNodeInElm,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,ARRAY2D<double> &dNdr,ARRAY2D<double> &dNdx,const double weight,double (&stress)[3][3],const bool mainLoop);
  //fem_hyperFoam_spatialForm.cpp
  // void calcStressTensor_hyperFoam_element_spatialForm_hexa(const int &ic,ARRAY2D<double> &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  void calcLambda(double (&stretch)[3],double (&stretchDirection)[3][3],const double (&C)[3][3]);
  void summation_postProcess(double &averageLambda,const double (&stress)[3][3],const double (&F)[3][3],const double Ic4bar,const int ic);
};

#endif