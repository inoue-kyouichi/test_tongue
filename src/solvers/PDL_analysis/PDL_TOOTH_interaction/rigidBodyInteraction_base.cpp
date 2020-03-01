/**
 * @file NRscheme_rigidBodyInteraction.cpp
 * @brief RigidElasticInteraction class
 * @author T. Otani
 */

#include "rigidElasticInteraction.h"

using namespace std;


// #################################################################
/**
 * @brief rigid body term
 * @param [in] RBdy          rigid body class
 */
void RigidElasticInteraction_base::calcRigidBodyInteractionTerm(DOUBLEARRAY2D &U,const RigidBody &RBdy)
{
  updateb(RBdy);
  tildeRB(RBdy);
  calcKqq(RBdy);
  calc_Qlambda(U,RBdy);
  calc_Q_rigid(RBdy);
}

// #################################################################
/**
 * @brief update b
 * @param [in] RBdy          rigid body class
 */
void RigidElasticInteraction_base::updateb(const RigidBody &RBdy)
{
  for(int ic=0;ic<numOfCP;ic++){
    for(int i=0;i<3;i++){
      b(ic,i)=0e0;
      for(int j=0;j<3;j++) b(ic,i)+=RBdy.R[i][j]*b0(ic,j);
    }
  }
}

// #################################################################
/**
 * @brief skew-symmetric matrix determined by b
 * @param [in] RBdy          rigid body class
 */
void RigidElasticInteraction_base::tildeRB(const RigidBody &RBdy)
{
  double tmp[3],Rbtmp[3][3];
  for(int ic=0;ic<numOfCP;ic++){
    for(int i=0;i<3;i++) tmp[i]=b(ic,i);
    mathTool::skewSymmetricTensor(Rbtmp,tmp);
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) Rb(ic,i,j)=Rbtmp[i][j];
    }
  }
}

// #################################################################
/**
 * @brief calc Kqq of Kww
 * @param [in] RBdy          rigid body class
 */
void RigidElasticInteraction_base::calcKqq(const RigidBody &RBdy)
{
  double lambda_tmp[3],matrix[3][3],TildeLambda[3][3];

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) Kqq[i][j]=0e0;
  }

  for(int ic=0;ic<numOfCP;ic++){

   for(int i=0;i<3;i++) lambda_tmp[i]=LAMBDA(ic,i);
    mathTool::skewSymmetricTensor(TildeLambda,lambda_tmp);

    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        matrix[i][j]=0e0;
        for(int k=0;k<3;k++){
          matrix[i][j]+=TildeLambda[i][k]*Rb(ic,k,j);
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
void RigidElasticInteraction_base::calc_Qlambda(DOUBLEARRAY2D &U,const RigidBody &RBdy)
{
  for(int ic=0;ic<numOfCP;ic++){
    for(int j=0;j<3;j++) Qlambda(ic,j)=U(CP(ic),j)-RBdy.U[j]+b0(ic,j)-b(ic,j);
  }
}

// #################################################################
/**
 * @brief calc Q_rigid vector
 * @param [in] RBdy          rigid body class
 */
void RigidElasticInteraction_base::calc_Q_rigid(const RigidBody &RBdy)
{
  for(int i=0;i<3;i++) QU[i]=0e0;
  for(int ic=0;ic<numOfCP;ic++){
    for(int i=0;i<3;i++) QU[i]-=LAMBDA(ic,i);
  }

  for(int i=0;i<3;i++) Qw[i]=0e0;
  for(int ic=0;ic<numOfCP;ic++){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) Qw[i]-=Rb(ic,i,j)*LAMBDA(ic,j);
    }
  }
}

// #################################################################
/**
 * @brief temporal fw.
 */
void RigidElasticInteraction_base::calcTemporalFw(const RigidBody &RBdy)
{
  double momentArm[3];
  for(int i=0;i<3;i++){
    momentArm[i]=0e0;
    for(int j=0;j<3;j++) momentArm[i] += RBdy.R[i][j] * initialMomentArm[j];
  }
  //printf("momentArm=%e %e %e\n",momentArm[0],momentArm[1],momentArm[2]);

  double tmp;
  mathTool::crossProduct(momentArm,FU,Fw,tmp);
}

// #################################################################
/**
 * @brief corrector scheme.
 * @param [in] u           displacement vector
 * @param [in] relaxation  relaxation parameters
 */
void RigidElasticInteraction_base::corrector_statics(DOUBLEARRAY2D &U,const double *u,RigidBody &RBdy,const int numOfNode,const double relaxation)
{
  double w[3];

  #pragma omp parallel for
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) U(i,j) += u[i+j*numOfNode]*relaxation;
  }

  #pragma omp parallel for
  for(int i=0;i<numOfCP;i++){
    for(int j=0;j<3;j++) LAMBDA(i,j) += u[3*numOfNode+i+j*numOfCP]*relaxation;
  }

  for(int j=0;j<3;j++) RBdy.U[j] += u[3*numOfNode+3*numOfCP+j]*relaxation;
  for(int j=0;j<3;j++) w[j]       = u[3*numOfNode+3*numOfCP+3+j]*relaxation;

  RBdy.updateRotationMatrix_spatialForm(w);
}
