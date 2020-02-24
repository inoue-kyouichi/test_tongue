#ifndef _FEM_H_
#define _FEM_H_


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
 * @file   fem.h
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

#include "allocation.h"
#include "domain.h"
#include "gauss.h"
#include "math_tools.h"
#include "fem_define.h"
#include "fem_base_mathTool.h"
#include "ShapeFunction.h"
#include "fileIO.h"
#include "pardiso_solver.h"
#include "rigidBody.h"


class Fem : public Domain{

  //fem.cpp
 public:
  Fem(){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) I2[i][j]=0e0;
        I2[i][i]=1e0;
    }

    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        for(int k=0;k<3;k++){
          for(int l=0;l<3;l++){
            I4[i][j][k][l]=5e-1*(I2[i][k]*I2[j][l]+I2[i][l]*I2[j][k]);
          }
        }
      }
    }
  }

  double I2[3][3],I4[3][3][3][3];

  double totalVolume;

  DOUBLEARRAY2D U, innerForce, externalForce, RHS;
  DOUBLEARRAY5D Ku,K;
  DOUBLEARRAY1D volume,volume0,volumeChangeRatio;

  void rotationalDirichlet(const int loop);
  void calcStressTensor();
  void set_rhs_statics();
  void calc_MassMatrix();
  void corrector_statistics(const double *u,const double relaxation);
  void calcVolume_hexa(const int &ic,DOUBLEARRAY1D &elementVolume,const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);

  void stress_tensor_initialize();

 private:
  double rho;
  int ConstitutiveLawName;
  DOUBLEARRAY3D BFe;
  DOUBLEARRAY3D Qu;
  DOUBLEARRAY3D Mass;
  INTARRAY1D boundaryNode_femur,boundaryNode_tibia;
  void exportRestartData(const int loop);

  //line search
  double line_search(const double *u);
  double line_search_innerProduct(DOUBLEARRAY2D &Q,const double *u);
  void calc_Q(DOUBLEARRAY2D &innerForce_tmp,DOUBLEARRAY2D &U_tmp);

  //fem_preprocessing.cpp
 public:
  void initialize(TextParser &tp);
  void allocate();
 private:
  void inputDomainInfo(TextParser &tp);
  void inputMaterialInfo(TextParser &tp);
  void inputDirichletBoundaryInfo(TextParser &tp);
  void inputNeumannBoundaryInfo(TextParser &tp);
  void inputFiberInfo(TextParser &tp);
  void restart_setting(const int dataNumber,const bool Restart,TextParser &tp);
  void calc_normal_quad(double (&normal)[3],double (&X)[4][3]);
  void setFiberDirection();
  void setFiberDirection_KogaModel(TextParser &tp);
  double arcLength(const double xMin,const double xMax);

  //fem_postprocessing.cpp
  public:
  DOUBLEARRAY1D Mises;
  DOUBLEARRAY2D AEigen_Ave,sigmaEigen_Ave;
  DOUBLEARRAY3D AEigenVector_Ave,sigmaEigenVector_Ave;
  void postProcess_PDL_element_spatialForm_hexa_SRI(const int &ic,DOUBLEARRAY2D &U_tmp,
  const int &numOfNodeInElm,const int &numOfGaussPoint);
  private:
  void calcEigen(const double (&A)[3][3],double (&AEigen)[3],double (&AEigenVector)[3][3]);
  void normalize(DOUBLEARRAY2D &AEigen,DOUBLEARRAY3D &AEigenVector_Ave,const int ic);


  public:
  double bulkModulusRatio;
  DOUBLEARRAY3D fiberDirection;
  DOUBLEARRAY2D lambda_ave;

  //fem_PDL_spatialForm2018.cpp
public:
  DOUBLEARRAY2D fiberDirection_elm;
  void postProcess_PDL_element_2018(const int &ic,DOUBLEARRAY2D &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint);
  void postProcess_LinearElastic_element_spatialForm(const int ic,const bool option);
  void calcStressTensor_LinearElastic_element_spatialForm(const int ic,const bool option);
private:
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

  //fem_boundary.cpp
//  public:
//   double boundaryPressure;
//   void calc_surfaceBoundaryForce();
//  private:
//   void calc_TractionByPressure_element(const int &ic);
//   void calc_Traction_element_quad(const int &ic,const int numOfNodeInBdElm,const int numOfGaussPoint);
//   void calc_TractionByPressure_element(const int &ic,const Element &boundaryElement);
//   void calcBFe_inGaussIntegral(DOUBLEARRAY3D &BFe,DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &X,const int numOfNodeInBdElm,const double weight,const int ic);

};

#endif //_FEM_H_
