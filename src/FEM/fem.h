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

#include "fem.h"
#include "domain.h"
#include "gauss.h"
#include "math_tools.h"
#include "fem_define.h"
#include "ShapeFunction.h"
#include "fileIO.h"
#include "pardiso_solver.h"
#include "rigidBody.h"

// class SurfaceSet{
//  public:
//   int numOfBoundaryNode,numOfNode,numOfElm;
//   INTARRAY2 elm;
//   INTARRAY1 boundaryNode;
//   DOUBLEARRAY2 x;

//   void setBoundaryNode(const int numberOfTotalNode);
//   void readDAT(const std::string &file,const int elementNumber);
//   void readVTU(const std::string &file,const int nodeNumber,const int elementNumber);
//   void translation(const double (&center)[3]);
//   void rotation(const double (&q)[4]);
//   void exportVTU(const std::string &file);
// };

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

  // SurfaceSet internalSurface,externalSurface;
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
  DOUBLEARRAY2 U, innerForce, externalForce, RHS;
  DOUBLEARRAY5 Ku,K;
  DOUBLEARRAY1 volume,volume0,volumeChangeRatio;
  DOUBLEARRAY1 bundle;
  INTARRAY1 bundleElement;
  double AMBinitialStretch,PLBinitialStretch;
  double tibiaRotation,femurRotation;

  void rotationalDirichlet(const int loop);
  int NRscheme(PARDISO_solver &PARDISO,RigidBody &RBdy);
  void calcStressTensor();
  void set_rhs_statics();
  void calc_MassMatrix();
  void corrector_statistics(const double *u,const double relaxation);
  void femSolidAnalysis(PARDISO_solver &PARDISO,RigidBody &RBdy);
  void calcVolume_hexa(const int &ic,DOUBLEARRAY1 &elementVolume,const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);

  void calc_B_matrix(DOUBLEARRAY3 &B,const DOUBLEARRAY2 &dNdr,const double (&dXdr)[3][3],
            const DOUBLEARRAY2 &u,const int &numOfNodeInElm);
  void stress_tensor_initialize();
  void calcTemporalFw(RigidBody &RBdy);


 private:
  double rho;
  int ConstitutiveLawName;
  DOUBLEARRAY3 BFe;
  DOUBLEARRAY3 Qu;
  DOUBLEARRAY3 Mass;
  INTARRAY1 boundaryNode_femur,boundaryNode_tibia;
  void exportRestartData(const int loop);
  void calc_thetaFromRotationMatrix(double (&ql)[3],const double (&R)[3][3]);

  //line search
  double line_search(const double *u);
  double line_search_innerProduct(const DOUBLEARRAY2 &Q,const double *u);
  void calc_Q(DOUBLEARRAY2 &innerForce_tmp,const DOUBLEARRAY2 &U_tmp);

  //fem_preprocessing.cpp
 public:
  void initialize();
  void allocate();
 private:
  void inputDomainInfo();
  // void planeDirichletCondition(SurfaceSet &plane,const INTARRAY1 &boundaryNode);
  // void planeDirichletCondition(SurfaceSet &plane,const int num);
  void inputDirichletBoundaryInfo();
  void inputNeumannBoundaryInfo();
  void inputFiberInfo();
  void inputSolverInfo();
  void inputOutputInfo();
  void restart_setting();
  // void read_boundaryNode(const std::string &file,INTARRAY1 &boundaryNode);
  void calc_normal_quad(double (&normal)[3],double (&X)[4][3]);
  void setFiberDirection();
  void setFiberDirection_KogaModel();
  void calcRotationMatrix(double (&R)[3][3],const double (&rotAxis)[3],const double angle);
  double arcLength(const double xMin,const double xMax);

  //fem_postprocessing.cpp
  DOUBLEARRAY2 AEigen_Ave,sigmaEigen_Ave;
  DOUBLEARRAY3 AEigenVector_Ave,sigmaEigenVector_Ave;

  void postProcess_PDL_element_spatialForm_hexa_SRI(const int &ic,const DOUBLEARRAY2 &U_tmp,
  const int &numOfNodeInElm,const int &numOfGaussPoint);
  void postProcess_ACL_element_spatialForm_hexa_SRI(const int &ic,const DOUBLEARRAY2 &U_tmp,
  const int &numOfNodeInElm,const int &numOfGaussPoint);
  void postProcess_ACL_element_spatialForm_hexa_Fbar(const int &ic,const DOUBLEARRAY2 &U_tmp,
    const int &numOfNodeInElm,const int &numOfGaussPoint);
  void calcEigen(const double (&A)[3][3],double (&AEigen)[3],double (&AEigenVector)[3][3]);
  void normalize(const DOUBLEARRAY2 &AEigen,DOUBLEARRAY3 &AEigenVector_Ave,const int ic);

  //fem_SantVenant_spatialForm.cpp
  void calcStressTensor_SantVenant_element_spatialForm(const int ic,const DOUBLEARRAY2 &U_tmp,const bool option);
  double SantVenant_inGaussIntegral(const DOUBLEARRAY2 &dNdr,const DOUBLEARRAY2 &x_current,const DOUBLEARRAY2 &x_ref,
        DOUBLEARRAY2 &dNdx,const int numOfNodeInElm,const double weight,const int ic,const bool option);
  void tensorPushForward_4order(double (&c4)[3][3][3][3],const double (&C4)[3][3][3][3],const double (&F)[3][3],const double J);

  //fem_neoHookean_spatialForm.cpp
  void calcStressTensor_NeoHookean_element_spatialForm_hexa(const int &ic,const DOUBLEARRAY2 &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  void calcStressTensor_NeoHookean_element_spatialForm_hexa_Fbar(const int &ic,const DOUBLEARRAY2 &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);

  //fem_ACL_spatialForm.cpp
  double bulkModulusRatio;
  DOUBLEARRAY3 fiberDirection;
  DOUBLEARRAY2 lambda_ave;
  void calcStressTensor_ACL_element_spatialForm_hexa_SRI(const int &ic,const DOUBLEARRAY2 &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  void calcStressTensor_ACL_element_spatialForm_hexa_Fbar(const int &ic,const DOUBLEARRAY2 &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  void calcStressTensor_ACL_element_spatialForm(const int ic,const DOUBLEARRAY2 &U_tmp,const bool option);
  void ACL_inGaussIntegral(const DOUBLEARRAY2 &dNdr,const DOUBLEARRAY2 &x_current,const DOUBLEARRAY2 &x_ref,DOUBLEARRAY2 &dNdx,const int numOfNodeInElm,const double weight,const int ic,const bool option);
  void calc_F_initial(double (&F_initial)[3][3],const double (&a0)[3],const double lambda);

  //fem_PDL_spatialForm.cpp
public:
  DOUBLEARRAY2 fiberDirection_elm;
private:
  void calcStressTensor_PDL_element_spatialForm_hexa_Fbar(const int &ic,const DOUBLEARRAY2 &U_tmp,
  const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  void calcStressTensor_PDL_element_spatialForm_hexa_SRI(const int &ic,const DOUBLEARRAY2 &U_tmp,
  const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  void calcStressTensor_PDL_element_spatialForm_hexa_2018(const int &ic,const DOUBLEARRAY2 &U_tmp,
  const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  void calcStressTensor_PDL_element_spatialForm_hexa_SRI_2018(const int &ic,const DOUBLEARRAY2 &U_tmp,
  const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  int calcStressTensor_PDL_element_fibreStretch(const int &ic,const DOUBLEARRAY2 &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint);

  void calcStressTensor_PDL_element_2018(const int &ic,const DOUBLEARRAY2 &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint);
  void calcStressTensor_PDL_element_spatialForm_hexa_2018_inGaussIntegral(const int &ic,const DOUBLEARRAY2 &U_tmp,
      const int &numOfNodeInElm,const Gauss &gauss,const DOUBLEARRAY2 &x_current,const DOUBLEARRAY2 &x_ref,const DOUBLEARRAY2 &dNdr,const DOUBLEARRAY2 &dNdx,const int i1,const int i2,const int i3);
  void calcStressTensor_hyperFoam_element_spatialForm_hexa_inGaussIntegral(const int &ic,const DOUBLEARRAY2 &U_tmp,
  const int &numOfNodeInElm,const Gauss &gauss,const DOUBLEARRAY2 &x_current,const DOUBLEARRAY2 &x_ref,const DOUBLEARRAY2 &dNdr,const DOUBLEARRAY2 &dNdx,const int i1,const int i2,const int i3);

  //fem_hyperFoam_spatialForm.cpp
  void calcStressTensor_hyperFoam_element_spatialForm_hexa(const int &ic,const DOUBLEARRAY2 &U_tmp,
  const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);
  void calcLambda(double (&stretch)[3],double (&stretchDirection)[3][3],const double (&C)[3][3]);

  //fem_boundary.cpp
 public:
  double boundaryPressure;
  void calc_surfaceBoundaryForce();
 private:
  void calc_TractionByPressure_element(const int &ic);
  void calc_Traction_element_quad(const int &ic,const int numOfNodeInBdElm,const int numOfGaussPoint);
  void calc_TractionByPressure_element(const int &ic,const Element &boundaryElement);
  void calcBFe_inGaussIntegral(DOUBLEARRAY3 &BFe,const DOUBLEARRAY2 &dNdr,
      const DOUBLEARRAY2 &X,const int numOfNodeInBdElm,const double weight,const int ic);

  //fem_base.cpp
 private:
  void calc_dudX(double (&dudX)[3][3],const DOUBLEARRAY2 &dNdX,const DOUBLEARRAY2 &u,const int &numOfNodeInElm);
  void calc_dXdr(double (&dXdr)[3][3],const DOUBLEARRAY2 &dNdr,const DOUBLEARRAY2 &x0,const int &numOfNodeInElm);
  void calc_dxdr(double (&dxdr)[3][3],const DOUBLEARRAY2 &dNdr,const DOUBLEARRAY2 &x,const int &numOfNodeInElm);
  void calc_dNdX(DOUBLEARRAY2 &dNdX,const DOUBLEARRAY2 &dNdr,const double (&dXdr)[3][3],const int &numOfNodeInElm);
  void calc_dNdx(DOUBLEARRAY2 &dNdx,const DOUBLEARRAY2 &dNdr,const double (&dxdr)[3][3],const int &numOfNodeInElm);

  //rigidBodyInteraction
 public:
  int numOfCP;
  double FU[3],Fw[3],FU_input[3],FUpoint[3],initialMomentArm[3];
  double Kqq[3][3],QU[3],Qw[3];
  INTARRAY1 CP,iCP;
  DOUBLEARRAY2 b,b0,Qlambda;
  DOUBLEARRAY2 LAMBDA;
  DOUBLEARRAY3 Rb;
  void preprocess_rigidBodyInteraction(const RigidBody &RBdy);

private:
  void inputRigidBodyInterface();
  void rigidBodyInteraction(const RigidBody &RBdy);
  void corrector_statics(const double *u,const double relaxation,RigidBody &RBdy);

  void updateb(const RigidBody &RBdy);
  void tildeRB(const RigidBody &RBdy);
  void calcKqq(const RigidBody &RBdy);
  void updateRotationMatrix_spatialForm(double (&R)[3][3],const double (&w)[3]);
  void calc_Qlambda(const RigidBody &RBdy);
  void calc_Q_rigid(const RigidBody &RBdy);

};

#endif //_FEM_H_
