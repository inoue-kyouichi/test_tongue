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
  void postProcess_LinearElastic_element_spatialForm(const int ic,const bool option);
  void calcStressTensor_LinearElastic_element_spatialForm(const int ic,const bool option);

  void inputFiberInfo(TextParser &tp);

  void setFiberDirection();
  void setFiberDirection_KogaModel(TextParser &tp);

  // void postProcess_PDL_element_spatialForm_hexa_SRI(const int &ic,DOUBLEARRAY2D &U_tmp,
  // const int &numOfNodeInElm,const int &numOfGaussPoint);
  void calcEigen(const double (&A)[3][3],double (&AEigen)[3],double (&AEigenVector)[3][3]);
  void normalize(DOUBLEARRAY2D &AEigen,DOUBLEARRAY3D &AEigenVector_Ave,const int ic);

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

  double LinearElastic_inGaussIntegral(DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &x_ref,const int numOfNodeInElm,const double weight,const int ic,const double lambda,const double mu,const bool option);
  double postProcess_LinearElastic_inGaussIntegral(double (&sigmaEigen)[3],double (&sigmaEigenVector)[3][3],DOUBLEARRAY2D &u,DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &x_ref,DOUBLEARRAY2D &dNdx,const int numOfNodeInElm,const double weight,const int ic,const double lambda,const double mu,const bool option);

};

#endif