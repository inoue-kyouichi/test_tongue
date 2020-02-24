#ifndef _RIGID_ELASTIC_INTERACTION_H_
#define _RIGID_ELASTIC_INTERACTION_H_

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

#include "PDL.h"
#include "RigidBody.h"

class RigidElasticInteraction_base {
 public:


  int numOfCP;
  double FU[3],Fw[3],FU_input[3],FUpoint[3],initialMomentArm[3];
  double Kqq[3][3],QU[3],Qw[3];
  INTARRAY1D CP,iCP;
  DOUBLEARRAY2D b,b0,Qlambda;
  DOUBLEARRAY2D LAMBDA;
  DOUBLEARRAY3D Rb;

  void calcTemporalFw(const RigidBody &RBdy);
  void updateRotationMatrix_spatialForm(double (&R)[3][3],const double (&w)[3]);
  void corrector_statics(DOUBLEARRAY2D &U,const double *u,RigidBody &RBdy,const int numOfNode,const double relaxation);
  void calc_thetaFromRotationMatrix(double (&ql)[3],const double (&R)[3][3]);
  void calcRigidBodyInteractionTerm(DOUBLEARRAY2D &U,const RigidBody &RBdy);
 private:
  void updateb(const RigidBody &RBdy);
  void tildeRB(const RigidBody &RBdy);
  void calcKqq(const RigidBody &RBdy);
  void calc_Qlambda(DOUBLEARRAY2D &U,const RigidBody &RBdy);
  void calc_Q_rigid(const RigidBody &RBdy);

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
  void inputSolverInfo(TextParser &tp);
  void inputOutputInfo(TextParser &tp);

};

#endif //_FEM_H_
