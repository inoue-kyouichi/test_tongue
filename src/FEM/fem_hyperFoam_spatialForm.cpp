
/**
 * @file fem_PDL_spatialForm.cpp
 * @brief Fem class
 * @author T. Otani
 */
#include <Eigen/Core>
#include <Eigen/Eigen>

using namespace Eigen;

#include "fem.h"
using namespace std;

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
 */
void Fem::calcStressTensor_hyperFoam_element_spatialForm_hexa(const int &ic,const DOUBLEARRAY2 &U_tmp,
const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option)
{
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  double sigma[3][3],S[3][3],C[3][3],invC[3][3],F[3][3];
  double elasticityTensor_ref[3][3][3][3],elasticityTensor_current[3][3][3][3],tangentCoefficient[3][3][3][3];

  DOUBLEARRAY2 x_current=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 x_ref=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 dNdr=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 dNdx=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);

  //nearly-incompressible material
  double pressure,dpressure,Siso[3][3],Sp[3];
  double invC_odot[3][3][3][3],C4iso[3][3][3][3],C4vol[3][3][3][3];

  //hyperFoam
  double term;
  double stretch[3],stretchDirection[3][3],stretchDiff;
  double alpha=20.9e0;
  double mu=3e-2*1e6;
  double poisson=0.257e0;
  double beta=poisson/(1e0-2e0*poisson);

  //------two point---------
  Gauss gauss(numOfGaussPoint);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      x_current[p][i] = x0[element[ic].node[p]][i]+U_tmp[element[ic].node[p]][i];
      x_ref[p][i]     = x0[element[ic].node[p]][i];
    }
  }

  for(int i1=0;i1<numOfGaussPoint;i1++){
    for(int i2=0;i2<numOfGaussPoint;i2++){
      for(int i3=0;i3<numOfGaussPoint;i3++){

        switch(element[ic].meshType){
          case VTK_HEXAHEDRON:
            ShapeFunction::C3D8_dNdr(dNdr,gauss.point[i1],gauss.point[i2],gauss.point[i3]);
            break;
          case VTK_TRIQUADRATIC_HEXAHEDRON:
            ShapeFunction::C3D27_dNdr(dNdr,gauss.point[i1],gauss.point[i2],gauss.point[i3]);
            break;
          default:
            cout << "error " << endl;
            cout << element[ic].meshType << endl;
            break;
        }

        calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
        calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
        detJ = mathTool::calcDeterminant_3x3(dxdr);
        volume += detJ * gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];

        calc_dXdr(dXdr,dNdr,x_ref,numOfNodeInElm);
        mathTool::calcInverseMatrix_3x3(drdX,dXdr);

        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            F[i][j]=0e0;
            for(int k=0;k<3;k++) F[i][j] += dxdr[i][k]*drdX[k][j];
          }
        }

        J = mathTool::calcDeterminant_3x3(F);

        if(J<0e0){
          cout << "Jacobian is nagtive. Exit..." << endl;
          exit(1);
        }

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

        //calc_internal force vector
        for(int p=0;p<numOfNodeInElm;p++){
          for(int i=0;i<3;i++){
            for(int j=0;j<3;j++) Qu[ic][p][i] += sigma[i][j] * dNdx[p][j] * detJ * gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
          }
        }

        if(option==false) continue;

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

        tensorPushForward_4order(elasticityTensor_current,C4iso,F,J);

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
                    K[ic][p][q][i][j] += dNdx[p][k]*(tangentCoefficient[i][j][k][l]*dNdx[q][l]) * detJ * gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
                  }
                }
              }
            }
          }
        }

      }
    }
  }
  // if(averageLambda/8e0>1.1e0) cout << "test" << endl;
  // printf("volume=%e\n",volume);
  Allocation::free2d(x_current);
  Allocation::free2d(x_ref);
  Allocation::free2d(dNdr);
  Allocation::free2d(dNdx);
}

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
