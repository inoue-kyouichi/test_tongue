/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */
#include "PDL.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

// #################################################################
/**
 * @brief calc stress tensor of SantVenant material
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] option           true or faluse: calculate tangential stiffness matrix or not.
 */
void PeriodontalLigament::postProcess_LinearElastic_element_spatialForm(const int ic,const bool option)
{
  double young = 66.7e-3; //(MPa)
  double poisson = 0.49e0;
  double PDL_lambda = young * poisson / ((1e0+poisson) * (1e0-2e0*poisson));
  double PDL_mu = 5e-1 * young / (1e0+poisson);

  young = 2.1e-3; //(MPa)
  poisson = 0.45e0;
  double PULP_lambda = young * poisson / ((1e0+poisson) * (1e0-2e0*poisson));
  double PULP_mu = 5e-1 * young / (1e0+poisson);
  // double young = 66.7e0; //(MPa)
  // double poisson = 0.49e0;
  // double lambda = young * poisson / ((1e0+poisson) * (1e0-2e0*poisson));
  // double mu = 5e-1 * young / (1e0+poisson);

  int numOfNodeInElm=element[ic].node.size();
  DOUBLEARRAY2D x_ref(numOfNodeInElm,3);
  DOUBLEARRAY2D u(numOfNodeInElm,3);
  DOUBLEARRAY2D dNdr(numOfNodeInElm,3);
  DOUBLEARRAY2D dNdX(numOfNodeInElm,3);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      x_ref(p,i) = x0(element[ic].node[p],i);
      u(p,i)     = U(element[ic].node[p],i);
    }
  }

  Gauss g(1),g2(2);
  GaussTetra gTet(1),gTet2(2);
  GaussTriangle gTri(1),gTri2(2);
  double sigmaEigen[3],sigmaEigenVector[3][3];

  switch(element[ic].meshType){
    // case VTK_TETRA:
    //   ShapeFunction::C3D4_dNdr(dNdr,gTet.point[0][0],gTet.point[0][1],gTet.point[0][2],gTet.point[0][3]);
    //   SantVenant_inGaussIntegral(dNdr,x_current,x_ref,dNdx,numOfNodeInElm,gTet.weight[0]*1e0/6e0,ic,option);
    //   break;
    // case VTK_HEXAHEDRON:
    //   for(int i1=0;i1<2;i1++){
    //     for(int i2=0;i2<2;i2++){
    //       for(int i3=0;i3<2;i3++){
    //         ShapeFunction::C3D8_dNdr(dNdr,g.point[i1],g.point[i2],g.point[i3]);
    //         SantVenant_inGaussIntegral(dNdr,x_current,x_ref,dNdx,numOfNodeInElm,g.weight[i1]*g.weight[i2]*g.weight[i3],ic,option);
    //       }
    //     }
    //   }
    //   break;
    case VTK_QUADRATIC_TETRA:
      for(int i=0;i<3;i++){
        sigmaEigen_Ave(ic,i)=0e0;
        for(int j=0;j<3;j++) sigmaEigenVector_Ave(ic,i,j)=0e0;
      }
      for(int i1=0;i1<4;i1++){
        ShapeFunction3D::C3D10_dNdr(dNdr,gTet2.point[i1][0],gTet2.point[i1][1],gTet2.point[i1][2],gTet2.point[i1][3]);
        if(element[ic].materialType==PDL){
          postProcess_LinearElastic_inGaussIntegral(sigmaEigen,sigmaEigenVector,u,dNdr,x_ref,dNdX,numOfNodeInElm,gTet2.weight[i1]*1e0/6e0,ic,PDL_lambda,PDL_mu,option);
        }else if(element[ic].materialType==PULP){
          postProcess_LinearElastic_inGaussIntegral(sigmaEigen,sigmaEigenVector,u,dNdr,x_ref,dNdX,numOfNodeInElm,gTet2.weight[i1]*1e0/6e0,ic,PULP_lambda,PULP_mu,option);
        }
        for(int i=0;i<3;i++){
          sigmaEigen_Ave(ic,i)+=sigmaEigen[i];
          for(int j=0;j<3;j++){
            sigmaEigenVector_Ave(ic,i,j)+=sigmaEigenVector[i][j];
          }
        }
      }
      for(int i=0;i<3;i++) sigmaEigen_Ave(ic,i)/=4e0;
      Mises(ic)=sqrt(5e-1*(pow(sigmaEigen_Ave(ic,0)-sigmaEigen_Ave(ic,1),2e0)+pow(sigmaEigen_Ave(ic,1)-sigmaEigen_Ave(ic,2),2e0)+pow(sigmaEigen_Ave(ic,2)-sigmaEigen_Ave(ic,0),2e0)));
      normalize(sigmaEigen_Ave,sigmaEigenVector_Ave,ic);
      break;
    case VTK_QUADRATIC_HEXAHEDRON:
      break;
    // case VTK_WEDGE:
    //   for(int i1=0;i1<2;i1++){
    //     ShapeFunction::C3D6_dNdr(dNdr,gTri.point[0][0],gTri.point[0][1],gTri.point[0][2],g.point[i1]);
    //     SantVenant_inGaussIntegral(dNdr,x_current,x_ref,dNdx,numOfNodeInElm,5e-1*gTri.weight[0]*g.weight[i1],ic,option);
    //   }
    //   break;
    default:
      cout << "undefined mesh type" << endl;
      exit(1);
  }
}

// #################################################################
/**
 * @brief calc stress tensor of SantVenant material
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] numOfNodeInElm   number of node in each element
 * @param [in] numOfGaussPoint  number of Gauss point set in each element
 * @param [in] option           true or faluse: calculate tangential stiffness matrix or not.
 */
double PeriodontalLigament::postProcess_LinearElastic_inGaussIntegral(double (&sigmaEigen)[3],double (&sigmaEigenVector)[3][3],DOUBLEARRAY2D &u,DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &x_ref,DOUBLEARRAY2D &dNdx,const int numOfNodeInElm,const double weight,const int ic,const double lambda,const double mu,const bool option)
{
  double detJ,volume,J;
  double dXdr[3][3],C4[3][3][3][3];

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) C4[i][j][k][l] = lambda*I2[i][j]*I2[k][l]+mu*2e0*I4[i][j][k][l];
      }
    }
  }

  FEM_MathTool::calc_dXdr(dXdr,dNdr,x_ref,numOfNodeInElm);
  FEM_MathTool::calc_dNdx(dNdx,dNdr,dXdr,numOfNodeInElm);
  detJ = mathTool::calcDeterminant_3x3(dXdr);
  volume = detJ * weight;

  double dudx[3][3];
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      dudx[i][j]=0e0;
      for(int p=0;p<numOfNodeInElm;p++) dudx[i][j]+=dNdx(p,j)*u(p,i);
    }
  }

  double strain[3][3];
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) strain[i][j]=5e-1*(dudx[i][j]+dudx[j][i]);
  }

  double sigma[3][3];
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      sigma[i][j]=0e0;
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) sigma[i][j]+=C4[i][j][k][l]*strain[k][l];
      }
    }
  }

  calcEigen(sigma,sigmaEigen,sigmaEigenVector);

  if(option==false) return volume;
  return volume;
}

// #################################################################
/**
 * @brief normalize vectors
 */
// void Fem::normalize(DOUBLEARRAY2D &AEigen,DOUBLEARRAY3D &AEigenVector_Ave,const int ic)
// {
//   double tmp=0e0;
//   for(int i=0;i<3;i++){
//     tmp=     AEigenVector_Ave(ic,i,0)*AEigenVector_Ave(ic,i,0)
//             +AEigenVector_Ave(ic,i,1)*AEigenVector_Ave(ic,i,1)
//             +AEigenVector_Ave(ic,i,2)*AEigenVector_Ave(ic,i,2);
//     if(tmp<1e-15) tmp=1e0;
//     tmp=sqrt(tmp);
//     for(int j=0;j<3;j++) AEigenVector_Ave(ic,i,j)=AEigenVector_Ave(ic,i,j)/tmp*AEigen(ic,i);
//   }
// }

// #################################################################
/**
 * @brief calc eigenvalues and eigenVectors using Eigen library
 */
// void Fem::calcEigen(const double (&A)[3][3],double (&AEigen)[3],double (&AEigenVector)[3][3])
// {
//   int order[3];
//   double E[3][3],Vec[3][3];

//   Matrix3d M;
//   for(int i=0;i<3;i++){
//     for(int j=0;j<3;j++) M(i,j)=A[i][j];
//   }
//   SelfAdjointEigenSolver<Matrix3d> ES(M);
//   // cout << "The eigenvalues of M=\n" << ES.eigenvalues() << endl;
//   // cout << "The corresponding eigenvectors of M=\n"<< ES.eigenvectors() << endl;
//   // Vector3d min_eigen_vector2 = ES.eigenvectors().col(2);
//   // cout << "The␣eigenvector␣corresponding␣the␣minimum␣eigenvalue␣=␣"
//   //  << min_eigen_vector2.transpose() << endl;

//   AEigen[0] = ES.eigenvalues()(2);  //max
//   AEigen[1] = ES.eigenvalues()(1);  //med
//   AEigen[2] = ES.eigenvalues()(0);  //min
//   // cout << "The␣minimum␣eigenvalue␣=␣" << min_eigen << endl;
//   Vector3d max_eigen_vector = ES.eigenvectors().col(2); //max
//   Vector3d med_eigen_vector = ES.eigenvectors().col(1);  //med
//   Vector3d min_eigen_vector = ES.eigenvectors().col(0);  //min
//   for(int j=0;j<3;j++){
//     AEigenVector[0][j]=max_eigen_vector(j);
//     AEigenVector[1][j]=med_eigen_vector(j);
//     AEigenVector[2][j]=min_eigen_vector(j);
//   }
// }