#ifndef _RBD_EBD_INTERACTION_H_
#define _RBD_EBD_INTERACTION_H_

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

#include "fem.h"


class RigidElasticInteraction : public Fem {

  //rigidBodyInteraction
 public:
  int numOfCP;
  double FU[3],Fw[3],FU_input[3],FUpoint[3],initialMomentArm[3];
  double Kqq[3][3],QU[3],Qw[3];
  INTARRAY1D CP,iCP;
  DOUBLEARRAY2D b,b0,Qlambda;
  DOUBLEARRAY2D LAMBDA;
  DOUBLEARRAY3D Rb;
  void preprocess_rigidBodyInteraction(const RigidBody &RBdy);
  void inputRigidBodyInterface();

  void femSolidAnalysis(PARDISO_solver &PARDISO,RigidBody &RBdy);
  int NRscheme(PARDISO_solver &PARDISO,RigidBody &RBdy);
  void calcTemporalFw(RigidBody &RBdy);
private:
  void rigidBodyInteraction(const RigidBody &RBdy);
  void corrector_statics(const double *u,const double relaxation,RigidBody &RBdy);
  void calc_thetaFromRotationMatrix(double (&ql)[3],const double (&R)[3][3]);
  void updateb(const RigidBody &RBdy);
  void tildeRB(const RigidBody &RBdy);
  void calcKqq(const RigidBody &RBdy);
  void updateRotationMatrix_spatialForm(double (&R)[3][3],const double (&w)[3]);
  void calc_Qlambda(const RigidBody &RBdy);
  void calc_Q_rigid(const RigidBody &RBdy);

};

#endif //_FEM_H_
