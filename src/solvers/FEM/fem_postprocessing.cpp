
/**
 * @file fem_PDL_spatialForm2018.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "fem.h"

#include <Eigen/Core>
#include <Eigen/Eigen>

using namespace Eigen;
using namespace std;

// #################################################################
/**
 * @brief normalize vectors
 */
void Fem::normalize(DOUBLEARRAY2D &AEigen,DOUBLEARRAY3D &AEigenVector_Ave,const int ic)
{
  double tmp=0e0;
  for(int i=0;i<3;i++){
    tmp=     AEigenVector_Ave(ic,i,0)*AEigenVector_Ave(ic,i,0)
            +AEigenVector_Ave(ic,i,1)*AEigenVector_Ave(ic,i,1)
            +AEigenVector_Ave(ic,i,2)*AEigenVector_Ave(ic,i,2);
    if(tmp<1e-15) tmp=1e0;
    tmp=sqrt(tmp);
    for(int j=0;j<3;j++) AEigenVector_Ave(ic,i,j)=AEigenVector_Ave(ic,i,j)/tmp*AEigen(ic,i);
  }
}
// #################################################################
/**
 * @brief calc eigenvalues and eigenVectors using Eigen library
 */
void Fem::calcEigen(const double (&A)[3][3],double (&AEigen)[3],double (&AEigenVector)[3][3])
{
  int order[3];
  double E[3][3],Vec[3][3];

  Matrix3d M;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) M(i,j)=A[i][j];
  }
  SelfAdjointEigenSolver<Matrix3d> ES(M);
  // cout << "The eigenvalues of M=\n" << ES.eigenvalues() << endl;
  // cout << "The corresponding eigenvectors of M=\n"<< ES.eigenvectors() << endl;
  // Vector3d min_eigen_vector2 = ES.eigenvectors().col(2);
  // cout << "The␣eigenvector␣corresponding␣the␣minimum␣eigenvalue␣=␣"
  //  << min_eigen_vector2.transpose() << endl;

  AEigen[0] = ES.eigenvalues()(2);  //max
  AEigen[1] = ES.eigenvalues()(1);  //med
  AEigen[2] = ES.eigenvalues()(0);  //min
  // cout << "The␣minimum␣eigenvalue␣=␣" << min_eigen << endl;
  Vector3d max_eigen_vector = ES.eigenvectors().col(2); //max
  Vector3d med_eigen_vector = ES.eigenvectors().col(1);  //med
  Vector3d min_eigen_vector = ES.eigenvectors().col(0);  //min
  for(int j=0;j<3;j++){
    AEigenVector[0][j]=max_eigen_vector(j);
    AEigenVector[1][j]=med_eigen_vector(j);
    AEigenVector[2][j]=min_eigen_vector(j);
  }
}