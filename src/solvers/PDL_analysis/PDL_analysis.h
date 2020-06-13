#ifndef _PDL_ANALYSIS_H_
#define _PDL_ANALYSIS_H_

//##################################################################################
//
// FEM Base
//
// Copyright (c) 2019 Mechanical and Bioengineering Systems Lab.,
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

// #include "PDL.h"
// #include "Rat_PDL.h"
#include "fem.h"
#include "RigidBody.h"

class RigidElasticInteraction_base {
 public:

  int numOfCP;
  double FU[3],Fw[3],FU_input[3],FUpoint[3],initialMomentArm[3];
  double Kqq[3][3],QU[3],Qw[3];
  ARRAY1D<int> CP,iCP;
  ARRAY2D<double> b,b0,Qlambda;
  ARRAY2D<double> LAMBDA;
  ARRAY3D<double> Rb;

  void calcTemporalFw(const RigidBody &RBdy);
  void updateRotationMatrix_spatialForm(double (&R)[3][3],const double (&w)[3]);
  void corrector_statics(ARRAY2D<double> &U,const double *u,RigidBody &RBdy,const int numOfNode,const double relaxation);
  void calc_thetaFromRotationMatrix(double (&ql)[3],const double (&R)[3][3]);
  void calcRigidBodyInteractionTerm(ARRAY2D<double> &U,const RigidBody &RBdy);
 private:
  void updateb(const RigidBody &RBdy);
  void tildeRB(const RigidBody &RBdy);
  void calcKqq(const RigidBody &RBdy);
  void calc_Qlambda(ARRAY2D<double> &U,const RigidBody &RBdy);
  void calc_Q_rigid(const RigidBody &RBdy);

};

namespace humanPDL{
class PeriodontalLigament:public Fem{
 public:
  double bulkModulusRatio;
  ARRAY3D<double> fiberDirection;
  ARRAY2D<double> lambda_ave;
  ARRAY1D<double> fibreStress;
  ARRAY1D<double> angleVariation;
  ARRAY1D<double> fibreAngle;

  ARRAY2D<double> fiberDirection_elm;

  void calcStressTensor();
  void allocatePDLvariables();

  void postProcess_PDL_element_2018(const int &ic,ARRAY2D<double> &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint);

  void inputFiberInfo(TextParser &tp);

  void setFiberDirection();
  void setFiberDirection_KogaModel(TextParser &tp);

  void export_vtu(const std::string &file);
  void export_vtu_Mises(const std::string &file);
  void export_vtu_boundary(const std::string &file);
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

class RigidElasticInteraction :public RigidElasticInteraction_base{

  //rigidBodyInteraction
 public:
  PeriodontalLigament ElasticBody;
  void initialize_rigidBodyInteraction();
  void mainLoop();

  TextParser tp;
  std::string outputDir,fileName;

 private:

  int dataNumber;
  int Restart;
  int OMPnumThreads;
  int maxIteration,NRiteration;
  double NRtolerance;
  double relaxation;

  RigidBody RBdy;
  PARDISO_solver PARDISO;

  void inputRigidBodyInterface();
  int NRscheme();
  void set_rhs_statics();
  void inputSolverInfo(TextParser &tp);
  void inputOutputInfo(TextParser &tp);
};
}

namespace ratPDL{

typedef enum {
  PULP       = 0,
  PDL        = 1,
  NUMBER_OF_MATERIALTYPE
} MATERIALType;

class PeriodontalLigament : public Fem{
 public:
  void postProcess_LinearElastic_element_spatialForm(const int ic,const bool option);
  void calcStressTensor_LinearElastic_element_spatialForm(const int ic,const bool option);
  void inputMaterialInfo(TextParser &tp);

 private:
  double LinearElastic_inGaussIntegral(ARRAY2D<double> &dNdr,ARRAY2D<double> &x_ref,const int numOfNodeInElm,const double weight,const int ic,const double lambda,const double mu,const bool option);
  double postProcess_LinearElastic_inGaussIntegral(double (&sigmaEigen)[3],double (&sigmaEigenVector)[3][3],ARRAY2D<double> &u,ARRAY2D<double> &dNdr,ARRAY2D<double> &x_ref,ARRAY2D<double> &dNdx,const int numOfNodeInElm,const double weight,const int ic,const double lambda,const double mu,const bool option);

};


class RigidElasticInteraction :public RigidElasticInteraction_base{

  //rigidBodyInteraction
 public:
  PeriodontalLigament ElasticBody;
  void initialize_rigidBodyInteraction();
  void mainLoop();

  TextParser tp;
  std::string outputDir,fileName;

 private:

  int dataNumber;
  int Restart;
  int OMPnumThreads;
  int maxIteration,NRiteration;
  double NRtolerance;
  double relaxation;

  RigidBody RBdy;
  PARDISO_solver PARDISO;

  void inputRigidBodyInterface();
  int NRscheme();
  void set_rhs_statics();
  void inputSolverInfo(TextParser &tp);
  void inputOutputInfo(TextParser &tp);
};
}


#endif //_FEM_H_
