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
 * @brief calc stress tensor of PDL with using selective reduced integration
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] numOfNodeInElm   number of node in each element
 * @param [in] numOfGaussPoint  number of Gauss point set in each element
 * @detail
   PDL model and parameters: Ortun-Terrazas et al., J. Mech. Behavior Biomed. Mat., 2018
 */
void Fem::calcStressTensor_PDL_element_2018(const int &ic,DOUBLEARRAY2D &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint)
{
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  double C[3][3],F[3][3];
  double stress[3][3];

  Gauss gauss(numOfGaussPoint);
  DOUBLEARRAY2D x_current(numOfNodeInElm,3);
  DOUBLEARRAY2D x_ref(numOfNodeInElm,3);
  DOUBLEARRAY2D dNdr(numOfNodeInElm,3);
  DOUBLEARRAY2D dNdx(numOfNodeInElm,3);

  double Ic4bar;
  double term4,term4_2,a0[3],a[3];

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      x_current(p,i) = x0(element[ic].node[p],i)+U_tmp(element[ic].node[p],i);
      x_ref(p,i)     = x0(element[ic].node[p],i);
    }
  }

  //--- Selective reduced integration ---
  for(int i1=0;i1<numOfGaussPoint;i1++){
    for(int i2=0;i2<numOfGaussPoint;i2++){
      for(int i3=0;i3<numOfGaussPoint;i3++){

        ShapeFunction3D::C3D8_dNdr(dNdr,gauss.point[i1],gauss.point[i2],gauss.point[i3]);

        FEM_MathTool::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
        detJ = mathTool::calcDeterminant_3x3(dxdr);
        FEM_MathTool::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
        FEM_MathTool::calc_dXdr(dXdr,dNdr,x_ref,numOfNodeInElm);
        mathTool::calcInverseMatrix_3x3(drdX,dXdr);
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            F[i][j]=0e0;
            for(int k=0;k<3;k++) F[i][j] += dxdr[i][k]*drdX[k][j];
          }
        }
        J = mathTool::calcDeterminant_3x3(F);

        if(J<0e0){
          printf("Jacobian is nagtive (element number=%d). Exit...\n",ic);
          for(int p=0;p<numOfNodeInElm;p++){
            printf("%e %e %e\n",x_current(p,0),x_current(p,1),x_current(p,2));
          }
          exit(1);
        }

        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            C[i][j]=0e0;
            for(int k=0;k<3;k++) C[i][j]+=F[k][i]*F[k][j];
          }
        }

        //sigma an-isotropic term
        for(int i=0;i<3;i++) a0[i]=fiberDirection_elm(ic,i);

        Ic4bar=0e0;
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) Ic4bar += a0[i]*C[i][j]*a0[j]*pow(J,-2e0/3e0);
        }

        if(Ic4bar<1e0){
          calcStressTensor_hyperFoam_element_spatialForm_hexa_inGaussIntegral(ic,U,8,gauss,x_current,x_ref,dNdr,dNdx,i1,i2,i3,stress,true);
        }else{
          calcStressTensor_PDL_element_spatialForm_hexa_2018_inGaussIntegral(ic,U,8,gauss,x_current,x_ref,dNdr,dNdx,i1,i2,i3,stress,true);
        }

      }
    }
  }

}



// #################################################################
/**
 * @brief calc stress tensor of PDL with using selective reduced integration
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] numOfNodeInElm   number of node in each element
 * @param [in] numOfGaussPoint  number of Gauss point set in each element
 * @param [in] option           true or faluse: calculate tangential stiffness matrix or not.
 * @detail
   PDL model and parameters: Ortun-Terrazas et al., J. Mech. Behavior Biomed. Mat., 2018
 */
void Fem::postProcess_PDL_element_2018(const int &ic,DOUBLEARRAY2D &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint)
{
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  double C[3][3],F[3][3];
  double stress[3][3];

  Gauss gauss(numOfGaussPoint);
  DOUBLEARRAY2D x_current(numOfNodeInElm,3);
  DOUBLEARRAY2D x_ref(numOfNodeInElm,3);
  DOUBLEARRAY2D dNdr(numOfNodeInElm,3);
  DOUBLEARRAY2D dNdx(numOfNodeInElm,3);

  //PDL
  double Ic4bar;
  double lambda,averageLambda=0e0,tmp;
  // int fiberNum=0;
  double term4,term4_2,a0[3],a[3];

  for(int j=0;j<3;j++) lambda_ave(ic,j)=0e0;

  //post process
  double sigmaEigen[3],sigmaEigenVector[3][3];
  for(int i=0;i<3;i++){
    sigmaEigen_Ave(ic,i)=0e0;
    for(int j=0;j<3;j++){
      sigmaEigenVector_Ave(ic,i,j)=0e0;
    }
  }

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      x_current(p,i) = x0(element[ic].node[p],i)+U_tmp(element[ic].node[p],i);
      x_ref(p,i)     = x0(element[ic].node[p],i);
    }
  }

  //--- Selective reduced integration ---
  for(int i1=0;i1<numOfGaussPoint;i1++){
    for(int i2=0;i2<numOfGaussPoint;i2++){
      for(int i3=0;i3<numOfGaussPoint;i3++){

        ShapeFunction3D::C3D8_dNdr(dNdr,gauss.point[i1],gauss.point[i2],gauss.point[i3]);

        FEM_MathTool::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
        detJ = mathTool::calcDeterminant_3x3(dxdr);
        FEM_MathTool::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
        FEM_MathTool::calc_dXdr(dXdr,dNdr,x_ref,numOfNodeInElm);
        mathTool::calcInverseMatrix_3x3(drdX,dXdr);
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            F[i][j]=0e0;
            for(int k=0;k<3;k++) F[i][j] += dxdr[i][k]*drdX[k][j];
          }
        }
        J = mathTool::calcDeterminant_3x3(F);

        if(J<0e0){
          printf("Jacobian is nagtive (element number=%d). Exit...\n",ic);
          for(int p=0;p<numOfNodeInElm;p++){
            printf("%e %e %e\n",x_current(p,0),x_current(p,1),x_current(p,2));
          }
          exit(1);
        }

        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            C[i][j]=0e0;
            for(int k=0;k<3;k++) C[i][j]+=F[k][i]*F[k][j];
          }
        }

        //sigma an-isotropic term
        for(int i=0;i<3;i++) a0[i]=fiberDirection_elm(ic,i);

        Ic4bar=0e0;
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) Ic4bar += a0[i]*C[i][j]*a0[j]*pow(J,-2e0/3e0);
        }
        lambda=sqrt(Ic4bar);
        averageLambda+=lambda;

      if(Ic4bar<1e0){
        calcStressTensor_hyperFoam_element_spatialForm_hexa_inGaussIntegral(ic,U,8,gauss,x_current,x_ref,dNdr,dNdx,i1,i2,i3,stress,false);
      }else{
        calcStressTensor_PDL_element_spatialForm_hexa_2018_inGaussIntegral(ic,U,8,gauss,x_current,x_ref,dNdr,dNdx,i1,i2,i3,stress,false);
      }

        calcEigen(stress,sigmaEigen,sigmaEigenVector);
        for(int i=0;i<3;i++){
          sigmaEigen_Ave(ic,i)+=sigmaEigen[i];
          for(int j=0;j<3;j++){
            sigmaEigenVector_Ave(ic,i,j)+=sigmaEigenVector[i][j];
          }
        }

        for(int i=0;i<3;i++){
          a[i]=0e0;
          for(int j=0;j<3;j++) a[i] += F[i][j] * a0[j];
        }
        tmp=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
        for(int j=0;j<3;j++) lambda_ave(ic,j)+=sqrt(Ic4bar)*a[j]/tmp;
      }
    }
  }

  tmp=sqrt(lambda_ave(ic,0)*lambda_ave(ic,0)+lambda_ave(ic,1)*lambda_ave(ic,1)+lambda_ave(ic,2)*lambda_ave(ic,2));
  for(int j=0;j<3;j++) lambda_ave(ic,j)/=tmp;
  for(int j=0;j<3;j++) lambda_ave(ic,j)*=averageLambda/8e0;

  for(int j=0;j<3;j++) sigmaEigen_Ave(ic,j)/=8e0;
  normalize(sigmaEigen_Ave,sigmaEigenVector_Ave,ic);
}

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


// #################################################################
/**
 * @brief calc stress tensor of PDL with using selective reduced integration
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] numOfNodeInElm   number of node in each element
 * @param [in] numOfGaussPoint  number of Gauss point set in each element
 * @param [in] option           true or faluse: calculate tangential stiffness matrix or not.
 * @detail
   PDL model and parameters: Ortun-Terrazas et al., J. Mech. Behavior Biomed. Mat., 2018
 */
void Fem::calcStressTensor_PDL_element_spatialForm_hexa_2018_inGaussIntegral(const int &ic,DOUBLEARRAY2D &U_tmp,
const int &numOfNodeInElm,const Gauss &gauss,DOUBLEARRAY2D &x_current,DOUBLEARRAY2D &x_ref,DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &dNdx,const int i1,const int i2,const int i3,double (&stress)[3][3],const bool mainLoop)
{
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  double sigma[3][3],S[3][3],C[3][3],invC[3][3],F[3][3],B[3][3];
  double elasticityTensor_ref[3][3][3][3],elasticityTensor_current[3][3][3][3],tangentCoefficient[3][3][3][3];

  //nearly-incompressible material
  double pressure,dpressure,Siso[3][3],S_aniso[3][3],Sbar[3][3],P4[3][3][3][3],P4bar[3][3][3][3];
  double invC_odot[3][3][3][3],C4iso[3][3][3][3],C4vol[3][3][3][3],C4bar[3][3][3][3],term2;

  //PDL
  double Ic4bar;
  int fiberNum=0;
  double lambda;

  const double c1=1e-2*1e6;
  const double bulkModulus=1e6/9.078e0;
  const double k1=0.298e0*1e6;
  const double k2=1.525e0;
  double term4,term4_2,a0[3],a[3];
  double F_initial[3][3],Ftmp[3][3];

  FEM_MathTool::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
  detJ = mathTool::calcDeterminant_3x3(dxdr);
  volume += detJ * gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];

  FEM_MathTool::calc_dXdr(dXdr,dNdr,x_ref,numOfNodeInElm);
  mathTool::calcInverseMatrix_3x3(drdX,dXdr);

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      F[i][j]=0e0;
      for(int k=0;k<3;k++) F[i][j] += dxdr[i][k]*drdX[k][j];
    }
  }

  J = mathTool::calcDeterminant_3x3(F);

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      C[i][j]=0e0;
      for(int k=0;k<3;k++){
        C[i][j]+=F[k][i]*F[k][j];
      }
    }
  }

  mathTool::calcInverseMatrix_3x3(invC,C);

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) P4[i][j][k][l]=I4[i][j][k][l]-1e0/3e0*invC[i][j]*C[k][l];
      }
    }
  }

  //------------------------------specific routine-----------------------------
  //S_bar
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) Sbar[i][j]=2e0*c1*I2[i][j];
  }

  //sigma an-isotropic term
  for(int i=0;i<3;i++) a0[i]=fiberDirection_elm(ic,i);

  Ic4bar=0e0;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) Ic4bar += a0[i]*C[i][j]*a0[j]*pow(J,-2e0/3e0);
    // for(int j=0;j<3;j++) Ic4bar += a0[i]*C[i][j]*a0[j];
  }

  lambda=sqrt(Ic4bar);
  for(int i=0;i<3;i++){
    a[i]=0e0;
    for(int j=0;j<3;j++) a[i] += F[i][j] * a0[j];
  }
  for(int i=0;i<3;i++) a[i] = a[i]/lambda;

  if(Ic4bar<1e0){
    term4=0e0;
    term4_2=0e0;
  }else{
    term4=k1*(Ic4bar-1e0)*exp(k2*(Ic4bar-1e0)*(Ic4bar-1e0));
    term4_2=k1*(1e0+2e0*k2*(Ic4bar-1e0)*(Ic4bar-1e0))*exp(k2*(Ic4bar-1e0)*(Ic4bar-1e0));
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) Sbar[i][j]+=2e0*term4*a0[i]*a0[j];
  }

  //S_iso
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      Siso[i][j]=0e0;
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) Siso[i][j] += pow(J,-2e0/3e0)*P4[i][j][k][l]*Sbar[k][l];
      }
    }
  }

  //S_vol or hydrostatic pressure
  pressure=bulkModulus*(J-1e0/J);
  dpressure=bulkModulus*(1e0+1e0/(J*J));
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      sigma[i][j] = pressure*I2[i][j];
    }
  }

  //transformation
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) sigma[i][j] += F[i][k]*Siso[k][l]*F[j][l]/J;
      }
    }
  }

  if(mainLoop==false){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) stress[i][j]=sigma[i][j];
    }
    return;
  }

  //calc_internal force vector
  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) Qu(ic,p,i) += sigma[i][j] * dNdx(p,j) * detJ * gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
    }
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) invC_odot[i][j][k][l]=5e-1*(invC[i][k]*invC[j][l]+invC[i][l]*invC[j][k]);
      }
    }
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) P4bar[i][j][k][l]=invC_odot[i][j][k][l]-1e0/3e0*(invC[i][j]*invC[k][l]);
      }
    }
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) C4vol[i][j][k][l]=J*(pressure+J*dpressure)*invC[i][j]*invC[k][l]-2e0*J*pressure*invC_odot[i][j][k][l];
      }
    }
  }

  //------------------------------specific routine-----------------------------
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) C4bar[i][j][k][l]=term4_2*a0[i]*a0[j]*a0[k]*a0[l];
                  // for(int l=0;l<3;l++) C4bar[i][j][k][l]=0e0;
      }
    }
  }

  term2=0e0;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) term2+=pow(J,-2e0/3e0)*(Sbar[i][j]*C[i][j]);
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++){
          C4iso[i][j][k][l]=0e0;
          for(int m=0;m<3;m++){
            for(int n=0;n<3;n++){
              for(int p=0;p<3;p++){
                for(int q=0;q<3;q++) C4iso[i][j][k][l]+=P4[i][j][m][n]*C4bar[m][n][p][q]*P4[k][l][p][q];
              }
            }
          }
          C4iso[i][j][k][l] += 2e0/3e0*term2*P4bar[i][j][k][l]-2e0/3e0*(invC[i][j]*Siso[k][l]+Siso[i][j]*invC[k][l]);
          C4iso[i][j][k][l] += C4vol[i][j][k][l];
                // C4iso[i][j][k][l]+=term4_2*a0[i]*a0[j]*a0[k]*a0[l];
        }
      }
    }
  }

        //------------------------------end specific routine-----------------------------

  FEM_MathTool::tensorPushForward_4order(elasticityTensor_current,C4iso,F,J);

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++){
          tangentCoefficient[i][j][k][l] = elasticityTensor_current[i][k][j][l]+sigma[k][l]*I2[i][j];
        }
      }
    }
  }

        //-------------------------------------------------------

        //calc_tangential_stiffness_matrix
  for(int p=0;p<numOfNodeInElm;p++){
    for(int q=0;q<numOfNodeInElm;q++){
      for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
          for(int k=0;k<3;k++){
            for(int l=0;l<3;l++){
              K(ic,p,q,i,j) += dNdx(p,k)*(tangentCoefficient[i][j][k][l]*dNdx(q,l)) * detJ * gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
            }
          }
        }
      }
    }
  }
}


// #################################################################
/**
 * @brief calc stress tensor of hyperFoam
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] numOfNodeInElm   number of node in each element
 * @param [in] numOfGaussPoint  number of Gauss point set in each element
 * @param [in] option           true or faluse: calculate tangential stiffness matrix or not.
 * @detail
   PDL model and parameters: Bergomi et al., J. Biomech., 2011
   W=2mu/alpha*[lambda_1^alpha+lambda_2^alpha+lambda_3^alpha-3+1/beta*(J^(alpha*beta)-1)]
 */
void Fem::calcStressTensor_hyperFoam_element_spatialForm_hexa_inGaussIntegral(const int &ic,DOUBLEARRAY2D &U_tmp,
const int &numOfNodeInElm,const Gauss &gauss,DOUBLEARRAY2D &x_current,DOUBLEARRAY2D &x_ref,DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &dNdx,const int i1,const int i2,const int i3,double (&stress)[3][3],const bool mainLoop)
{
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  double sigma[3][3],S[3][3],C[3][3],invC[3][3],F[3][3];
  double elasticityTensor_ref[3][3][3][3],elasticityTensor_current[3][3][3][3],tangentCoefficient[3][3][3][3];

  //nearly-incompressible material
  double pressure,dpressure,Siso[3][3],Sp[3];
  double invC_odot[3][3][3][3],C4iso[3][3][3][3],C4vol[3][3][3][3];

  //hyperFoam
  double term;
  double stretch[3],stretchDirection[3][3],stretchDiff;
  const double alpha=20.9e0;
  const double mu=3e-2*1e6;
  const double poisson=0.257e0;
  const double beta=poisson/(1e0-2e0*poisson);

  FEM_MathTool::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
  detJ = mathTool::calcDeterminant_3x3(dxdr);
  // volume += detJ * gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];

  FEM_MathTool::calc_dXdr(dXdr,dNdr,x_ref,numOfNodeInElm);
  mathTool::calcInverseMatrix_3x3(drdX,dXdr);

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      F[i][j]=0e0;
      for(int k=0;k<3;k++) F[i][j] += dxdr[i][k]*drdX[k][j];
    }
  }

  J = mathTool::calcDeterminant_3x3(F);

  // if(J<0e0){
  //   cout << "Jacobian is nagtive. Exit..." << endl;
  //   exit(1);
  // }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      C[i][j]=0e0;
      for(int k=0;k<3;k++){
        C[i][j]+=F[k][i]*F[k][j];
      }
    }
  }

  mathTool::calcInverseMatrix_3x3(invC,C);

  //------------------------------specific routine-----------------------------
  pressure= -2e0*mu/alpha*pow(J,-1e0*alpha*beta-1e0);
  dpressure=-2e0*mu/alpha*(-1e0*alpha*beta-1e0)*pow(J,-1e0*alpha*beta-2e0);
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      sigma[i][j] = pressure*I2[i][j];
    }
  }

  calcLambda(stretch,stretchDirection,C);
  for(int a=0;a<3;a++) Sp[a]=(2e0*mu/alpha)*pow(stretch[a],alpha-1e0)/stretch[a];

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      Siso[i][j]=0e0;
      for(int a=0;a<3;a++){
        Siso[i][j]+=Sp[a]*stretchDirection[a][i]*stretchDirection[a][j];
      }
    }
  }
  //transformation
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) sigma[i][j] += F[i][k]*Siso[k][l]*F[j][l]/J;
      }
    }
  }
  //---------------------------------------------------------------------------

  if(mainLoop==false){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) stress[i][j]=sigma[i][j];
    }
    return;
  }

  //calc_internal force vector
  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) Qu(ic,p,i) += sigma[i][j] * dNdx(p,j) * detJ * gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
    }
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) invC_odot[i][j][k][l]=5e-1*(invC[i][k]*invC[j][l]+invC[i][l]*invC[j][k]);
      }
    }
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) C4iso[i][j][k][l]=J*(pressure+J*dpressure)*invC[i][j]*invC[k][l]-2e0*J*pressure*invC_odot[i][j][k][l];
      }
    }
  }


  //------------------------------specific routine-----------------------------
  //isotropic term
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++){
          for(int a=0;a<3;a++){
              C4iso[i][j][k][l]+=(2e0*mu/alpha)*(alpha-2e0)*pow(stretch[a],alpha-3e0)/stretch[a]
                                *stretchDirection[a][i]*stretchDirection[a][j]*stretchDirection[a][k]*stretchDirection[a][l];
          }
          for(int a=0;a<3;a++){
            for(int b=0;b<3;b++){
              if(a==b) continue;
              stretchDiff=stretch[b]*stretch[b]-stretch[a]*stretch[a];
              if(fabs(stretchDiff)<1e-15){
                term=(2e0*mu/alpha)*(alpha-2e0)*pow(stretch[b],alpha-4e0);
              }else{
                term=(Sp[b]-Sp[a])/stretchDiff;
              }
              C4iso[i][j][k][l]+=term*(stretchDirection[a][i]*stretchDirection[b][j]*stretchDirection[a][k]*stretchDirection[b][l]
                                      +stretchDirection[a][i]*stretchDirection[b][j]*stretchDirection[b][k]*stretchDirection[a][l]);
            }
          }
        }
      }
    }
  }
  //------------------------------end specific routine-----------------------------

  FEM_MathTool::tensorPushForward_4order(elasticityTensor_current,C4iso,F,J);

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++){
          tangentCoefficient[i][j][k][l] = elasticityTensor_current[i][k][j][l]+sigma[k][l]*I2[i][j];
        }
      }
    }
  }

  //-------------------------------------------------------

  //calc_tangential_stiffness_matrix
  for(int p=0;p<numOfNodeInElm;p++){
    for(int q=0;q<numOfNodeInElm;q++){
      for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
          for(int k=0;k<3;k++){
            for(int l=0;l<3;l++){
              K(ic,p,q,i,j) += dNdx(p,k)*(tangentCoefficient[i][j][k][l]*dNdx(q,l)) * detJ * gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
            }
          }
        }
      }
    }
  }
}

// #################################################################
/**
 * @brief calc stretching direction
 * @param [out] stretch    principal stretches
 * @param [out] stretch direction     associating direction vector
 * @param [in]  C     right Cauchy-Green deformation tensor
 */
void Fem::calcLambda(double (&stretch)[3],double (&stretchDirection)[3][3],const double (&C)[3][3])
{
  Matrix3d M;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) M(i,j)=C[i][j];
  }
  SelfAdjointEigenSolver<Matrix3d> ES(M);
  // cout << "The eigenvalues of M=\n" << ES.eigenvalues() << endl;
  // cout << "The corresponding eigenvectors of M=\n"<< ES.eigenvectors() << endl;
  // Vector3d min_eigen_vector2 = ES.eigenvectors().col(2);
  // cout << "The␣eigenvector␣corresponding␣the␣minimum␣eigenvalue␣=␣"
  //  << min_eigen_vector2.transpose() << endl;

  stretch[0] = sqrt(ES.eigenvalues()(2));  //max
  stretch[1] = sqrt(ES.eigenvalues()(1));  //med
  stretch[2] = sqrt(ES.eigenvalues()(0));  //min
  // cout << "The␣minimum␣eigenvalue␣=␣" << min_eigen << endl;
  Vector3d max_eigen_vector = ES.eigenvectors().col(2); //max
  Vector3d med_eigen_vector = ES.eigenvectors().col(1);  //med
  Vector3d min_eigen_vector = ES.eigenvectors().col(0);  //min
  for(int j=0;j<3;j++){
    stretchDirection[0][j]=max_eigen_vector(j);
    stretchDirection[1][j]=med_eigen_vector(j);
    stretchDirection[2][j]=min_eigen_vector(j);
  }
}
