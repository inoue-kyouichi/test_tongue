/**
 * @file line_search_innerProduct.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "fem.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

using namespace std;

// #################################################################
/**
 * @brief line search algorithm
 * @param [in] u    incremental displacement vector
 * @detail detail of the line seach algorithm is shown in Crisfield, Nonlinear finite element analysis of solid and structure, Vol.1,1991.
 * In Japanese, Hisada and Noguchi (1996) summarizes this technique.
 */
double Fem::line_search(const double *u)
{
  double s0,sn,sn1,beta_tmp=1e0,beta;

  DOUBLEARRAY2D U_tmp(numOfNode,3);
  DOUBLEARRAY2D innerForce_tmp(numOfNode,3);

  #pragma omp parallel for
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) U_tmp(i,j) = U(i,j)+u[i+j*numOfNode]*beta_tmp;
  }
  calc_Q(innerForce_tmp,U_tmp);
  s0=line_search_innerProduct(innerForce,u);
  sn=line_search_innerProduct(innerForce_tmp,u);
  double alpha=s0/sn;

  beta=beta_tmp;
  for(int loop=0;loop<5;loop++){
    printf("\nbeta=%e alpha=%e s0=%e sn=%e\n",beta_tmp,alpha,s0,sn);

    if(abs(alpha)>2e0){
      if(fabs(beta_tmp)<1e-3){
        beta_tmp=1e-3;
      }else if(beta_tmp>1e0){
        beta_tmp=1e0;
      }else{
        beta=beta_tmp;
      }
      beta=beta_tmp;
      break;
    }
    #pragma omp parallel for
    for(int i=0;i<numOfNode;i++){
      for(int j=0;j<3;j++) U_tmp(i,j) = U(i,j)+u[i+j*numOfNode]*beta_tmp;
    }

    calc_Q(innerForce_tmp,U_tmp);
    sn=line_search_innerProduct(innerForce_tmp,u);
    alpha=s0/sn;

    beta_tmp=s0*beta_tmp/(s0-sn);
  }

  #pragma omp parallel for
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) U(i,j) += beta * u[i+j*numOfNode];
  }

  double residual=0e0;
  #pragma omp parallel for reduction(+:residual)
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) residual += ibd(i,j)*innerForce_tmp(i,j)*innerForce_tmp(i,j);
  }

  return sqrt(residual);
}

// #################################################################
/**
 * @brief inner product (Q dot u)
 * @param [in] Q    inner force vector
 * @param [in] u    displacement vector
 */
double Fem::line_search_innerProduct(DOUBLEARRAY2D &Q,const double *u)
{
  double s=0e0;
  #pragma omp parallel for reduction(+:s)
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++){
      s += -ibd(i,j)*Q(i,j)*u[i+j*numOfNode];
    }
  }
  return s;
}

// #################################################################
/**
 * @brief calc innerForce
 * @param [out] innerForce_tmp   inner force vector
 * @param [in] U_tmp             displacement vector
 */
void Fem::calc_Q(DOUBLEARRAY2D &innerForce_tmp,DOUBLEARRAY2D &U_tmp)
{
  double elementVolume,volume;

  #pragma omp parallel for
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) innerForce_tmp(i,j) = 0e0;
  }

  #pragma omp parallel for
  for(int ic=0;ic<numOfElm;ic++){
    for(int p=0;p<element[ic].node.size();p++){
      for(int i=0;i<3;i++) Qu(ic,p,i) = 0e0;
    }
  }

  #pragma omp parallel for
  for(int ic=0;ic<numOfElm;ic++){
    // calcStressTensor_ACL_element_spatialForm_hexa_Fbar(ic,U_tmp,8,2,false);
    calcStressTensor_PDL_element_spatialForm_hexa_SRI(ic,U_tmp,8,2,false);
  }

  for(int ic=0;ic<numOfElm;ic++){
    for(int p=0;p<element[ic].node.size();p++){
      for(int i=0;i<3;i++){
        innerForce_tmp(element[ic].node[p],i) += Qu(ic,p,i);
      }
    }
  }
}

