/**
 * @file fem_rigidBodyInteraction.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "fem.h"

using namespace std;

// #################################################################
/**
 * @brief corrector scheme.
 * @param [in] u           displacement vector
 * @param [in] relaxation  relaxation parameters
 * @param [in] RBdy        Rigid body class
 */
void Fem::corrector_statics(const double *u,const double relaxation,RigidBody &RBdy)
{
  double w[3];

  #pragma omp parallel for
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) U[i][j] += u[i+j*numOfNode]*relaxation;
  }

  #pragma omp parallel for
  for(int i=0;i<numOfCP;i++){
    for(int j=0;j<3;j++) LAMBDA[i][j] += u[3*numOfNode+i+j*numOfCP]*relaxation;
  }

  for(int j=0;j<3;j++) RBdy.U[j] += u[3*numOfNode+3*numOfCP+j]*relaxation;
  for(int j=0;j<3;j++) w[j]       = u[3*numOfNode+3*numOfCP+3+j]*relaxation;

  RBdy.updateRotationMatrix_spatialForm(w);
}

// #################################################################
/**
 * @brief calc b0
 * @param [in] RBdy          rigid body class
 */
void Fem::preprocess_rigidBodyInteraction(const RigidBody &RBdy)
{
  for(int ic=0;ic<numOfCP;ic++){
    for(int j=0;j<3;j++) b0[ic][j]=x[CP[ic]][j]-RBdy.xg[j];
  }
}

// #################################################################
/**
 * @brief rigid body term
 * @param [in] RBdy          rigid body class
 */
void Fem::rigidBodyInteraction(const RigidBody &RBdy)
{
  updateb(RBdy);
  tildeRB(RBdy);
  calcKqq(RBdy);
  calc_Qlambda(RBdy);
  calc_Q_rigid(RBdy);
}

// #################################################################
/**
 * @brief update b
 * @param [in] RBdy          rigid body class
 */
void Fem::updateb(const RigidBody &RBdy)
{
  for(int ic=0;ic<numOfCP;ic++){
    for(int i=0;i<3;i++){
      b[ic][i]=0e0;
      for(int j=0;j<3;j++) b[ic][i]+=RBdy.R[i][j]*b0[ic][j];
    }
  }
}

// #################################################################
/**
 * @brief skew-symmetric matrix determined by b
 * @param [in] RBdy          rigid body class
 */
void Fem::tildeRB(const RigidBody &RBdy)
{
  double tmp[3],Rbtmp[3][3];
  for(int ic=0;ic<numOfCP;ic++){
    for(int i=0;i<3;i++) tmp[i]=b[ic][i];
    mathTool::skewSymmetricTensor(Rbtmp,tmp);
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) Rb[ic][i][j]=Rbtmp[i][j];
    }
  }
}

// #################################################################
/**
 * @brief calc Kqq of Kww
 * @param [in] RBdy          rigid body class
 */
void Fem::calcKqq(const RigidBody &RBdy)
{
  double lambda_tmp[3],matrix[3][3],TildeLambda[3][3];

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) Kqq[i][j]=0e0;
  }

  for(int ic=0;ic<numOfCP;ic++){

   for(int i=0;i<3;i++) lambda_tmp[i]=LAMBDA[ic][i];
    mathTool::skewSymmetricTensor(TildeLambda,lambda_tmp);

    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        matrix[i][j]=0e0;
        for(int k=0;k<3;k++){
          matrix[i][j]+=TildeLambda[i][k]*Rb[ic][k][j];
        }
      }
    }

    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) Kqq[i][j]-=matrix[i][j];
    }

  }
}

// #################################################################
/**
 * @brief calc Q_lambda vector
 * @param [in] RBdy          rigid body class
 */
void Fem::calc_Qlambda(const RigidBody &RBdy)
{
  for(int ic=0;ic<numOfCP;ic++){
    for(int j=0;j<3;j++) Qlambda[ic][j]=U[CP[ic]][j]-RBdy.U[j]+b0[ic][j]-b[ic][j];
  }
}

// #################################################################
/**
 * @brief calc Q_rigid vector
 * @param [in] RBdy          rigid body class
 */
void Fem::calc_Q_rigid(const RigidBody &RBdy)
{
  for(int i=0;i<3;i++) QU[i]=0e0;
  for(int ic=0;ic<numOfCP;ic++){
    for(int i=0;i<3;i++) QU[i]-=LAMBDA[ic][i];
  }

  for(int i=0;i<3;i++) Qw[i]=0e0;
  for(int ic=0;ic<numOfCP;ic++){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) Qw[i]-=Rb[ic][i][j]*LAMBDA[ic][j];
    }
  }
}
