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
  TextParser tp;
  std::string outputDir,fileName;

  int dataNumber;
  int Restart;
  int OMPnumThreads;
  int maxIteration,NRiteration;
  double NRtolerance;
  double totalVolume;
  double relaxation;

  int numOfGaussPoint;
  DOUBLEARRAY2D U, innerForce, externalForce, RHS;
  DOUBLEARRAY5D Ku,K;
  DOUBLEARRAY1D volume,volume0,volumeChangeRatio;
  DOUBLEARRAY1D bundle;
  INTARRAY1D bundleElement;
  double tibiaRotation,femurRotation;

  void rotationalDirichlet(const int loop);
  // int NRscheme(PARDISO_solver &PARDISO,RigidBody &RBdy);
  void calcStressTensor();
  void set_rhs_statics();
  void calc_MassMatrix();
  void corrector_statistics(const double *u,const double relaxation);
  // void femSolidAnalysis(PARDISO_solver &PARDISO,RigidBody &RBdy);
  void calcVolume_hexa(const int &ic,DOUBLEARRAY1D &elementVolume,const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);

  void calc_B_matrix(DOUBLEARRAY3D &B,DOUBLEARRAY2D &dNdr,const double (&dXdr)[3][3],DOUBLEARRAY2D &u,const int &numOfNodeInElm);
  void stress_tensor_initialize();
  // void calcTemporalFw(RigidBody &RBdy);


 private:
  double rho;
  int ConstitutiveLawName;
  DOUBLEARRAY3D BFe;
  DOUBLEARRAY3D Qu;
  DOUBLEARRAY3D Mass;
  INTARRAY1D boundaryNode_femur,boundaryNode_tibia;
  void exportRestartData(const int loop);
  // void calc_thetaFromRotationMatrix(double (&ql)[3],const double (&R)[3][3]);

  //line search
  double line_search(const double *u);
  double line_search_innerProduct(DOUBLEARRAY2D &Q,const double *u);
  void calc_Q(DOUBLEARRAY2D &innerForce_tmp,DOUBLEARRAY2D &U_tmp);

  //fem_preprocessing.cpp
 public:
  void initialize();
  void allocate();
 private:
  void inputDomainInfo();
  void inputDirichletBoundaryInfo();
  void inputNeumannBoundaryInfo();
  void inputFiberInfo();
  void inputSolverInfo();
  void inputOutputInfo();
  void restart_setting();
  void calc_normal_quad(double (&normal)[3],double (&X)[4][3]);
  void setFiberDirection();
  void setFiberDirection_KogaModel();
  void calcRotationMatrix(double (&R)[3][3],const double (&rotAxis)[3],const double angle);
  double arcLength(const double xMin,const double xMax);

  //fem_postprocessing.cpp
  public:
  DOUBLEARRAY2D AEigen_Ave,sigmaEigen_Ave;
  DOUBLEARRAY3D AEigenVector_Ave,sigmaEigenVector_Ave;
  void postProcess_PDL_element_spatialForm_hexa_SRI(const int &ic,DOUBLEARRAY2D &U_tmp,
  const int &numOfNodeInElm,const int &numOfGaussPoint);
  void postProcess_ACL_element_spatialForm_hexa_SRI(const int &ic,DOUBLEARRAY2D &U_tmp,
  const int &numOfNodeInElm,const int &numOfGaussPoint);
  void postProcess_ACL_element_spatialForm_hexa_Fbar(const int &ic,DOUBLEARRAY2D &U_tmp,
    const int &numOfNodeInElm,const int &numOfGaussPoint);
    private:
  void calcEigen(const double (&A)[3][3],double (&AEigen)[3],double (&AEigenVector)[3][3]);
  void normalize(DOUBLEARRAY2D &AEigen,DOUBLEARRAY3D &AEigenVector_Ave,const int ic);

  void tensorPushForward_4order(double (&c4)[3][3][3][3],const double (&C4)[3][3][3][3],const double (&F)[3][3],const double J);

  //fem_neoHookean_spatialForm.cpp
  void calcStressTensor_NeoHookean_element_spatialForm_hexa(const int &ic,DOUBLEARRAY2D &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  void calcStressTensor_NeoHookean_element_spatialForm_hexa_Fbar(const int &ic,DOUBLEARRAY2D &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);

  //fem_ACL_spatialForm.cpp
  public:
  double bulkModulusRatio;
  DOUBLEARRAY3D fiberDirection;
  DOUBLEARRAY2D lambda_ave;

  //fem_PDL_spatialForm.cpp
public:
  DOUBLEARRAY2D fiberDirection_elm;
private:
  void calcStressTensor_PDL_element_spatialForm_hexa_Fbar(const int &ic,DOUBLEARRAY2D &U_tmp,
  const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  void calcStressTensor_PDL_element_spatialForm_hexa_SRI(const int &ic,DOUBLEARRAY2D &U_tmp,
  const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  void calcStressTensor_PDL_element_spatialForm_hexa_2018(const int &ic,DOUBLEARRAY2D &U_tmp,
  const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  void calcStressTensor_PDL_element_spatialForm_hexa_SRI_2018(const int &ic,DOUBLEARRAY2D &U_tmp,
  const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  int calcStressTensor_PDL_element_fibreStretch(const int &ic,DOUBLEARRAY2D &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint);

  void calcStressTensor_PDL_element_2018(const int &ic,DOUBLEARRAY2D &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint);
  void calcStressTensor_PDL_element_spatialForm_hexa_2018_inGaussIntegral(const int &ic,DOUBLEARRAY2D &U_tmp,
      const int &numOfNodeInElm,const Gauss &gauss,DOUBLEARRAY2D &x_current,DOUBLEARRAY2D &x_ref,DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &dNdx,const int i1,const int i2,const int i3);
  void calcStressTensor_hyperFoam_element_spatialForm_hexa_inGaussIntegral(const int &ic,DOUBLEARRAY2D &U_tmp,
  const int &numOfNodeInElm,const Gauss &gauss,DOUBLEARRAY2D &x_current,DOUBLEARRAY2D &x_ref,DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &dNdx,const int i1,const int i2,const int i3);

  //fem_hyperFoam_spatialForm.cpp
  void calcStressTensor_hyperFoam_element_spatialForm_hexa(const int &ic,DOUBLEARRAY2D &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  void calcLambda(double (&stretch)[3],double (&stretchDirection)[3][3],const double (&C)[3][3]);

  //fem_boundary.cpp
 public:
  double boundaryPressure;
  void calc_surfaceBoundaryForce();
 private:
  void calc_TractionByPressure_element(const int &ic);
  void calc_Traction_element_quad(const int &ic,const int numOfNodeInBdElm,const int numOfGaussPoint);
  void calc_TractionByPressure_element(const int &ic,const Element &boundaryElement);
  void calcBFe_inGaussIntegral(DOUBLEARRAY3D &BFe,DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &X,const int numOfNodeInBdElm,const double weight,const int ic);

  //fem_base.cpp
 private:
  void calc_dudX(double (&dudX)[3][3],DOUBLEARRAY2D &dNdX,DOUBLEARRAY2D &u,const int &numOfNodeInElm);
  void calc_dXdr(double (&dXdr)[3][3],DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &x0,const int &numOfNodeInElm);
  void calc_dxdr(double (&dxdr)[3][3],DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &x,const int &numOfNodeInElm);
  void calc_dNdX(DOUBLEARRAY2D &dNdX,DOUBLEARRAY2D &dNdr,const double (&dXdr)[3][3],const int &numOfNodeInElm);
  void calc_dNdx(DOUBLEARRAY2D &dNdx,DOUBLEARRAY2D &dNdr,const double (&dxdr)[3][3],const int &numOfNodeInElm);

  //rigidBodyInteraction
//  public:
//   int numOfCP;
//   double FU[3],Fw[3],FU_input[3],FUpoint[3],initialMomentArm[3];
//   double Kqq[3][3],QU[3],Qw[3];
//   INTARRAY1D CP,iCP;
//   DOUBLEARRAY2D b,b0,Qlambda;
//   DOUBLEARRAY2D LAMBDA;
//   DOUBLEARRAY3D Rb;
//   void preprocess_rigidBodyInteraction(const RigidBody &RBdy);

// private:
//   void inputRigidBodyInterface();
//   void rigidBodyInteraction(const RigidBody &RBdy);
//   void corrector_statics(const double *u,const double relaxation,RigidBody &RBdy);

//   void updateb(const RigidBody &RBdy);
//   void tildeRB(const RigidBody &RBdy);
//   void calcKqq(const RigidBody &RBdy);
//   void updateRotationMatrix_spatialForm(double (&R)[3][3],const double (&w)[3]);
//   void calc_Qlambda(const RigidBody &RBdy);
//   void calc_Q_rigid(const RigidBody &RBdy);

};

#endif //_FEM_H_
