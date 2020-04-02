/**
 * @file math_tools.cpp
 * @brief mathTool class
 * @author T. Otani
 */

#include "math_tools.h"
#include <omp.h>
#include <string>
#include <iostream>
#include <cmath>
#include <random>

using namespace std;

// #################################################################
/**
 * @brief calc random function
 */
double mathTool::rnd()
{
  random_device rd;
  mt19937 mt(rd());
  uniform_real_distribution<double> score(0.0e0,1.0e0);
  return score(mt);
}

// #################################################################
/**
 * @brief calc norm of vector
 */
double mathTool::vectorNorm(const int &nump,ARRAY1D<double> &x)
{
  double norm=0e0;

  #pragma omp parallel for reduction(+:norm)
  for(int i=0;i<nump;i++){
    norm += fabs(x(i));
  }
  return norm;
}

// #################################################################
/**
 * @brief calc inner product
 */
double mathTool::innerProduct(const int &nump,ARRAY1D<double> &x,ARRAY1D<double> &y)
{
  double dot_p=0e0;

  #pragma omp parallel for reduction(+:dot_p)
  for(int i=0;i<nump;i++){
    dot_p += x(i) * y(i);
  }
  return dot_p;
}

// #################################################################
/**
 * @brief calc cross product
 */
void mathTool::crossProduct(const double (&a)[3],const double (&b)[3],double (&c)[3],double &dc)
{
  dc = 0e0;

  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];

  for(int j=0;j<3;j++) dc += (c[j]*c[j]);

  dc = sqrt(dc);
}
// #################################################################
/**
 * @brief calc determinant
 */
double mathTool::calcDeterminant_3x3(const double (&a)[3][3])
{
  double det  = a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2] + a[2][0] * a[0][1] * a[1][2]
              - a[2][0] * a[1][1] * a[0][2] - a[1][0] * a[0][1] * a[2][2] - a[0][0] * a[2][1] * a[1][2];
  return det;
}
// #################################################################
/**
 * @brief calc inverse matrix
 */
void mathTool::calcInverseMatrix_3x3(double (&inv_a)[3][3],const double (&a)[3][3])
{
  double det;

  det = calcDeterminant_3x3(a);

  inv_a[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
  inv_a[0][1] = a[0][2]*a[2][1] - a[0][1]*a[2][2];
  inv_a[0][2] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
  inv_a[1][0] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
  inv_a[1][1] = a[0][0]*a[2][2] - a[0][2]*a[2][0];
  inv_a[1][2] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
  inv_a[2][0] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
  inv_a[2][1] = a[0][1]*a[2][0] - a[0][0]*a[2][1];
  inv_a[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) inv_a[i][j] = inv_a[i][j] / det;
  }
}


void mathTool::calcMatrix_x_matrix4(double (&ans)[4][4],const double (&a)[4][4],const double (&b)[4][4])
{
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      ans[i][j] = 0e0;
      for(int k=0;k<4;k++) ans[i][j] += a[i][k] * b[k][j];
    }
  }

}

/*************************************************************/
/* 実対称行列の固有値・固有ベクトル（ヤコビ法）              */
/*      n : 次数                                             */
/*      ct : 最大繰り返し回数                                */
/*      eps : 収束判定条件                                   */
/*      A : 対象とする行列                                   */
/*      A1, A2 : 作業域（nxnの行列），A1の対角要素が固有値   */
/*      X1, X2 : 作業域（nxnの行列），X1の各列が固有ベクトル */
/*      return : =0 : 正常                                   */
/*               =1 : 収束せず                               */
/*      coded by Y.Suganuma                                  */
/*************************************************************/
int mathTool::Jacobi3x3(const int &ct, const double &eps,const double (&A)[3][3], double (&A1)[3][3],double (&X1)[3][3])
{
  int n=3;
  double A2[3][3],X2[3][3];
  double max, s, t, v, sn, cs;
  int i1, i2, k = 0, ind = 1, p = 0, q = 0;
          // 初期設定
  for (i1 = 0; i1 < n; i1++) {
    for (i2 = 0; i2 < n; i2++) {
      A1[i1][i2] = A[i1][i2];
      X1[i1][i2] = 0.0;
    }
    X1[i1][i1] = 1.0;
  }
          // 計算
  while (ind > 0 && k < ct) {
            // 最大要素の探索
    max = 0.0;
    for (i1 = 0; i1 < n; i1++) {
      for (i2 = 0; i2 < n; i2++) {
        if (i2 != i1) {
          if (fabs(A1[i1][i2]) > max) {
            max = fabs(A1[i1][i2]);
            p   = i1;
            q   = i2;
          }
        }
      }
    }
            // 収束判定
              // 収束した
    if (max < eps)
      ind = 0;
              // 収束しない
    else {
                // 準備
      s  = -A1[p][q];
      t  = 0.5 * (A1[p][p] - A1[q][q]);
      v  = fabs(t) / sqrt(s * s + t * t);
      sn = sqrt(0.5 * (1.0 - v));
      if (s*t < 0.0)
        sn = -sn;
      cs = sqrt(1.0 - sn * sn);
                // Akの計算
      for (i1 = 0; i1 < n; i1++) {
        if (i1 == p) {
          for (i2 = 0; i2 < n; i2++) {
            if (i2 == p)
              A2[p][p] = A1[p][p] * cs * cs + A1[q][q] * sn * sn -
                                       2.0 * A1[p][q] * sn * cs;
            else if (i2 == q)
              A2[p][q] = 0.0;
            else
              A2[p][i2] = A1[p][i2] * cs - A1[q][i2] * sn;
          }
        }
        else if (i1 == q) {
          for (i2 = 0; i2 < n; i2++) {
            if (i2 == q)
              A2[q][q] = A1[p][p] * sn * sn + A1[q][q] * cs * cs +
                                       2.0 * A1[p][q] * sn * cs;
            else if (i2 == p)
              A2[q][p] = 0.0;
            else
              A2[q][i2] = A1[q][i2] * cs + A1[p][i2] * sn;
          }
        }
        else {
          for (i2 = 0; i2 < n; i2++) {
            if (i2 == p)
              A2[i1][p] = A1[i1][p] * cs - A1[i1][q] * sn;
            else if (i2 == q)
              A2[i1][q] = A1[i1][q] * cs + A1[i1][p] * sn;
            else
              A2[i1][i2] = A1[i1][i2];
          }
        }
      }
                // Xkの計算
      for (i1 = 0; i1 < n; i1++) {
        for (i2 = 0; i2 < n; i2++) {
          if (i2 == p)
            X2[i1][p] = X1[i1][p] * cs - X1[i1][q] * sn;
          else if (i2 == q)
            X2[i1][q] = X1[i1][q] * cs + X1[i1][p] * sn;
          else
            X2[i1][i2] = X1[i1][i2];
        }
      }
                // 次のステップへ
      k++;
      for (i1 = 0; i1 < n; i1++) {
        for (i2 = 0; i2 < n; i2++) {
          A1[i1][i2] = A2[i1][i2];
          X1[i1][i2] = X2[i1][i2];
        }
      }
    }
  }

  return ind;
}
// #########################################################
/**
 * @brief skew symmetric tensor
 */
void mathTool::skewSymmetricTensor(double (&M)[3][3],const double (&v)[3])
{
  M[0][0]=0e0;   M[0][1]=-v[2]; M[0][2]=v[1];
  M[1][0]=v[2];  M[1][1]=0e0;   M[1][2]=-v[0];
  M[2][0]=-v[1]; M[2][1]=v[0];  M[2][2]=0e0;
}

// #########################################################
/**
 * @brief rotation tensor
 */
void mathTool::quaternion2rotation(double (&R)[3][3],const double (&q)[4])
{
  R[0][0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
  R[0][1]=2e0*(q[1]*q[2]-q[0]*q[3]);
  R[0][2]=2e0*(q[1]*q[3]+q[0]*q[2]);
  R[1][0]=2e0*(q[1]*q[2]+q[0]*q[3]);
  R[1][1]=q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
  R[1][2]=2e0*(q[2]*q[3]-q[0]*q[1]);
  R[2][0]=2e0*(q[1]*q[3]-q[0]*q[2]);
  R[2][1]=2e0*(q[2]*q[3]+q[0]*q[1]);
  R[2][2]=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
}

// #################################################################
/**
 * @brief calcRotationMatrix
 */
void mathTool::calcRotationMatrix(double (&R)[3][3],const double (&rotAxis)[3],const double angle)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) R[i][j]=0e0;
    R[i][i]=1e0;
  }

  double S[3][3],SS[3][3];
  mathTool::skewSymmetricTensor(S,rotAxis);
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      SS[i][j]=0e0;
      for(int k=0;k<3;k++) SS[i][j] += S[i][k] * S[k][j];
    }
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) R[i][j] += sin(angle) * S[i][j] + (1e0-cos(angle)) * SS[i][j];
  }
}

// #################################################################
/**
 * @brief calc theta (Nour-Omid and Rankin, Compt. Methods Appl. Mech. Eng., 1991.)
 * @param [out] ql      rotation angle in local coordinates
 * @param [in] Rbar     rotation matrix [reference to current coorindate (node level)]
 * @param [in] ic element number
 */
void mathTool::calc_thetaFromRotationMatrix(double (&ql)[3],const double (&R)[3][3])
{
  double trR,Ra[3][3],tau;

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      Ra[i][j] = R[i][j]-R[j][i];
    }
  }

  ql[0] = Ra[2][1];
  ql[1] = Ra[0][2];
  ql[2] = Ra[1][0];
  tau=5e-1*sqrt(ql[0]*ql[0]+ql[1]*ql[1]+ql[2]*ql[2]);
  if(tau<1e-15){
    for(int i=0;i<3;i++) ql[i]=0e0;
  }else{
    for(int i=0;i<3;i++) ql[i]=5e-1*asin(tau)/tau*ql[i];
  }
}