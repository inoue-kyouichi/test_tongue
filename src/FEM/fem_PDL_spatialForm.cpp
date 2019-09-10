
/**
 * @file fem_PDL_spatialForm.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "fem.h"
using namespace std;


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
int Fem::calcStressTensor_PDL_element_fibreStretch(const int &ic,const DOUBLEARRAY2 &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint)
{
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  double C[3][3],F[3][3];

  DOUBLEARRAY2 x_current=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 x_ref=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 dNdr=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 dNdx=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);

  //PDL
  double Ic4bar,lambda,tmp;
  int fiberNum=0;
  double term4,term4_2,a0[3],a[3];

  //------two point---------
  Gauss gauss0(numOfGaussPoint-1);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      x_current[p][i] = x0[element[ic].node[p]][i]+U_tmp[element[ic].node[p]][i];
      x_ref[p][i]     = x0[element[ic].node[p]][i];
    }
  }

  //--- Selective reduced integration ---
  for(int i1=0;i1<numOfGaussPoint-1;i1++){
    for(int i2=0;i2<numOfGaussPoint-1;i2++){
      for(int i3=0;i3<numOfGaussPoint-1;i3++){

        ShapeFunction::C3D8_dNdr(dNdr,gauss0.point[i1],gauss0.point[i2],gauss0.point[i3]);
        // ShapeFunction::C3D8_dNdr(dNdr,0e0,0e0,0e0);
        calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
        detJ = mathTool::calcDeterminant_3x3(dxdr);
        calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
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
          cout << "centroid Jacobian is nagtive. Exit..." << endl;
          for(int p=0;p<numOfNodeInElm;p++){
            printf("%e %e %e\n",x_current[p][0],x_current[p][1],x_current[p][2]);
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
        for(int i=0;i<3;i++) a0[i]=fiberDirection_elm[ic][i];

        Ic4bar=0e0;
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) Ic4bar += a0[i]*C[i][j]*a0[j]*pow(J,-2e0/3e0);
        }
        for(int i=0;i<3;i++){
          a[i]=0e0;
          for(int j=0;j<3;j++) a[i] += F[i][j] * a0[j];
        }
        tmp=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
        for(int j=0;j<3;j++) lambda_ave[ic][j]=sqrt(Ic4bar)*a[j]/tmp;
      }
    }
  }
  Allocation::free2d(x_current);
  Allocation::free2d(x_ref);
  Allocation::free2d(dNdr);
  Allocation::free2d(dNdx);

  if(Ic4bar-1e0<-1e-12){
    return 0;
  }else{
    return 1;
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
void Fem::calcStressTensor_PDL_element_spatialForm_hexa_2018(const int &ic,const DOUBLEARRAY2 &U_tmp,
const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option)
{
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  double sigma[3][3],S[3][3],C[3][3],invC[3][3],F[3][3],B[3][3];
  double elasticityTensor_ref[3][3][3][3],elasticityTensor_current[3][3][3][3],tangentCoefficient[3][3][3][3];

  DOUBLEARRAY2 x_current=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 x_ref=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 dNdr=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 dNdx=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);

  double c1=1e-2*1e6;
  double bulkModulus=1e6/9.078e0;

  //nearly-incompressible material
  double pressure,dpressure,Siso[3][3],S_aniso[3][3],Sbar[3][3],P4[3][3][3][3],P4bar[3][3][3][3];
  double invC_odot[3][3][3][3],C4iso[3][3][3][3],C4vol[3][3][3][3],C4bar[3][3][3][3],term2;

  //PDL
  double Ic4bar;
  int fiberNum=0;
  double lambda;
  double k1=0.298e0*1e6;
  double k2=1.525e0;
  double term4,term4_2,a0[3],a[3];
  // double averageLambda=0e0;
  // for(int i=0;i<3;i++) lambda_ave[ic][i]=0e0;
  double F_initial[3][3],Ftmp[3][3];

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
            B[i][j]=0e0;
            C[i][j]=0e0;
            for(int k=0;k<3;k++){
              B[i][j]+=F[i][k]*F[j][k];
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
        for(int i=0;i<3;i++) a0[i]=fiberDirection_elm[ic][i];

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
        // averageLambda+=lambda;
        // for(int i=0;i<3;i++) lambda_ave[ic][i]+=a[i];

        if(Ic4bar<1e0){
          term4=0e0;
          term4_2=0e0;
        }else{
          term4=k1*(Ic4bar-1e0)*exp(k2*(Ic4bar-1e0)*(Ic4bar-1e0));
          term4_2=k1*(1e0+2e0*k2*(Ic4bar-1e0)*(Ic4bar-1e0))*exp(k2*(Ic4bar-1e0)*(Ic4bar-1e0));
        }

        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) Sbar[i][j]+=2e0*term4*a0[i]*a0[j];
          // for(int j=0;j<3;j++) S_aniso[i][j]=term4*a0[i]*a0[j];
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
              // for(int l=0;l<3;l++) sigma[i][j] += F[i][k]*(Siso[k][l]+S_aniso[k][l])*F[j][l]/J;
              for(int l=0;l<3;l++) sigma[i][j] += F[i][k]*Siso[k][l]*F[j][l]/J;
            }
          }
        }

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
  // double tmp=sqrt(lambda_ave[ic][0]*lambda_ave[ic][0]+lambda_ave[ic][1]*lambda_ave[ic][1]+lambda_ave[ic][2]*lambda_ave[ic][2]);
  // for(int j=0;j<3;j++) lambda_ave[ic][j]/=tmp;
  // for(int j=0;j<3;j++) lambda_ave[ic][j]*=averageLambda/8e0;
  // if(averageLambda/8e0>1.1e0) cout << "test" << endl;
  // printf("volume=%e\n",volume);
  Allocation::free2d(x_current);
  Allocation::free2d(x_ref);
  Allocation::free2d(dNdr);
  Allocation::free2d(dNdx);
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
void Fem::calcStressTensor_PDL_element_spatialForm_hexa_SRI_2018(const int &ic,const DOUBLEARRAY2 &U_tmp,
const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option)
{
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  double sigma[3][3],S[3][3],C[3][3],invC[3][3],F[3][3],B[3][3];
  double elasticityTensor_ref[3][3][3][3],elasticityTensor_current[3][3][3][3],tangentCoefficient[3][3][3][3];

  DOUBLEARRAY2 x_current=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 x_ref=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 dNdr=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 dNdx=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);

  double c1=1e-2*1e6;
  double bulkModulus=1e6/9.078e0;

  //nearly-incompressible material
  double pressure,dpressure,Siso[3][3],S_aniso[3][3],Sbar[3][3],P4[3][3][3][3],P4bar[3][3][3][3];
  double invC_odot[3][3][3][3],C4iso[3][3][3][3],C4vol[3][3][3][3],C4bar[3][3][3][3],term2;

  //PDL
  double Ic4bar;
  int fiberNum=0;
  double lambda;
  double k1=0.298e0*1e6;
  double k2=1.525e0;
  double term4,term4_2,a0[3],a[3];
  // double averageLambda=0e0;
  // for(int i=0;i<3;i++) lambda_ave[ic][i]=0e0;
  double F_initial[3][3],Ftmp[3][3];

  //------two point---------
  Gauss gauss(numOfGaussPoint),gauss0(numOfGaussPoint-1);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      x_current[p][i] = x0[element[ic].node[p]][i]+U_tmp[element[ic].node[p]][i];
      x_ref[p][i]     = x0[element[ic].node[p]][i];
    }
  }

  //--- Selective reduced integration ---
  for(int i1=0;i1<numOfGaussPoint-1;i1++){
    for(int i2=0;i2<numOfGaussPoint-1;i2++){
      for(int i3=0;i3<numOfGaussPoint-1;i3++){

  ShapeFunction::C3D8_dNdr(dNdr,gauss0.point[i1],gauss0.point[i2],gauss0.point[i3]);
  // ShapeFunction::C3D8_dNdr(dNdr,0e0,0e0,0e0);
  calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
  detJ = mathTool::calcDeterminant_3x3(dxdr);
  calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
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
    cout << "centroid Jacobian is nagtive. Exit..." << endl;
    for(int p=0;p<numOfNodeInElm;p++){
      printf("%e %e %e\n",x_current[p][0],x_current[p][1],x_current[p][2]);
    }
    exit(1);
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      C[i][j]=0e0;
      for(int k=0;k<3;k++) C[i][j]+=F[k][i]*F[k][j];
    }
  }

  mathTool::calcInverseMatrix_3x3(invC,C);

  //S_vol or hydrostatic pressure
  pressure=bulkModulus*(J-1e0/J);
  dpressure=bulkModulus*(1e0+1e0/(J*J));
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      sigma[i][j] = pressure*I2[i][j];
    }
  }
  //calc_internal force vector
  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) Qu[ic][p][i] += sigma[i][j] * dNdx[p][j] * detJ * gauss0.weight[i1] * gauss0.weight[i2] * gauss0.weight[i3];
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
        for(int l=0;l<3;l++) C4vol[i][j][k][l]=J*(pressure+J*dpressure)*invC[i][j]*invC[k][l]-2e0*J*pressure*invC_odot[i][j][k][l];
      }
    }
  }

  tensorPushForward_4order(elasticityTensor_current,C4vol,F,J);

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++){
          tangentCoefficient[i][j][k][l] = elasticityTensor_current[i][k][j][l]+sigma[k][l]*I2[i][j];
        }
      }
    }
  }

  for(int p=0;p<numOfNodeInElm;p++){
    for(int q=0;q<numOfNodeInElm;q++){
      for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
          for(int k=0;k<3;k++){
            for(int l=0;l<3;l++){
              // K[ic][p][q][i][j] += dNdx[p][k]*(tangentCoefficient[i][j][k][l]*dNdx[q][l]) * detJ * 2e0 * 2e0 * 2e0;
              K[ic][p][q][i][j] += dNdx[p][k]*(tangentCoefficient[i][j][k][l]*dNdx[q][l]) * detJ * gauss0.weight[i1] * gauss0.weight[i2] * gauss0.weight[i3];
            }
          }
        }
      }
    }
  }
  }}}
  //------------

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
            B[i][j]=0e0;
            C[i][j]=0e0;
            for(int k=0;k<3;k++){
              B[i][j]+=F[i][k]*F[j][k];
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
        for(int i=0;i<3;i++) a0[i]=fiberDirection_elm[ic][i];

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
        // averageLambda+=lambda;
        // for(int i=0;i<3;i++) lambda_ave[ic][i]+=a[i];

        // term4=k1*exp(k2*(Ic4bar-1e0));
        if(Ic4bar<1e0){
          term4=0e0;
          term4_2=0e0;
        }else{
          term4=k1*(Ic4bar-1e0)*exp(k2*(Ic4bar-1e0)*(Ic4bar-1e0));
          term4_2=k1*(1e0+2e0*k2*(Ic4bar-1e0)*(Ic4bar-1e0))*exp(k2*(Ic4bar-1e0)*(Ic4bar-1e0));
        }

        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) Sbar[i][j]+=2e0*term4*a0[i]*a0[j];
          // for(int j=0;j<3;j++) S_aniso[i][j]=term4*a0[i]*a0[j];
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

        //transformation
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            sigma[i][j] = 0e0;
            for(int k=0;k<3;k++){
              // for(int l=0;l<3;l++) sigma[i][j] += F[i][k]*(Siso[k][l]+S_aniso[k][l])*F[j][l]/J;
              for(int l=0;l<3;l++) sigma[i][j] += F[i][k]*Siso[k][l]*F[j][l]/J;
            }
          }
        }

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
              for(int l=0;l<3;l++) P4bar[i][j][k][l]=invC_odot[i][j][k][l]-1e0/3e0*(invC[i][j]*invC[k][l]);
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
                // C4iso[i][j][k][l]+=term4_2*a0[i]*a0[j]*a0[k]*a0[l];
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
  // double tmp=sqrt(lambda_ave[ic][0]*lambda_ave[ic][0]+lambda_ave[ic][1]*lambda_ave[ic][1]+lambda_ave[ic][2]*lambda_ave[ic][2]);
  // for(int j=0;j<3;j++) lambda_ave[ic][j]/=tmp;
  // for(int j=0;j<3;j++) lambda_ave[ic][j]*=averageLambda/8e0;
  // if(averageLambda/8e0>1.1e0) cout << "test" << endl;
  // printf("volume=%e\n",volume);
  Allocation::free2d(x_current);
  Allocation::free2d(x_ref);
  Allocation::free2d(dNdr);
  Allocation::free2d(dNdx);
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
   PDL model and parameters: Natali et al., J. Biomech. Eng., 1996
 */
void Fem::calcStressTensor_PDL_element_spatialForm_hexa_SRI(const int &ic,const DOUBLEARRAY2 &U_tmp,
const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option)
{
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  double sigma[3][3],S[3][3],C[3][3],invC[3][3],F[3][3],B[3][3];
  double elasticityTensor_ref[3][3][3][3],elasticityTensor_current[3][3][3][3],tangentCoefficient[3][3][3][3];

  DOUBLEARRAY2 x_current=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 x_ref=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 dNdr=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 dNdx=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);

  double c1=4.3e-3*1e6;
  double alpha1=2.11e0;
  double r=37.25e0;
  double Kv=3.48e-4*1e6;
  double bulkModulus=Kv/(2e0+r*(r+1e0));

  //nearly-incompressible material
  double pressure,dpressure,Siso[3][3],Sbar[3][3],P4[3][3][3][3],P4bar[3][3][3][3];
  double invC_odot[3][3][3][3],C4iso[3][3][3][3],C4vol[3][3][3][3],C4bar[3][3][3][3],term2;

  //PDL
  double Ic1_tilda,Ic4;
  int fiberNum=0;
  double sigma_aniso[3][3],S_aniso[3][3];
  double lambda;
  double c3=4.07e-3*1e6;
  double alpha3=3.39e0;
  double term4,term4_2,a0[3],a[3];
  double averageLambda=0e0;
  for(int i=0;i<3;i++) lambda_ave[ic][i]=0e0;
  double F_initial[3][3],Ftmp[3][3];

  //------two point---------
  Gauss gauss(numOfGaussPoint),gauss0(numOfGaussPoint-1);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      x_current[p][i] = x0[element[ic].node[p]][i]+U_tmp[element[ic].node[p]][i];
      x_ref[p][i]     = x0[element[ic].node[p]][i];
    }
  }

  //--- Selective reduced integration ---
  for(int i1=0;i1<numOfGaussPoint-1;i1++){
    for(int i2=0;i2<numOfGaussPoint-1;i2++){
      for(int i3=0;i3<numOfGaussPoint-1;i3++){

  ShapeFunction::C3D8_dNdr(dNdr,gauss0.point[i1],gauss0.point[i2],gauss0.point[i3]);
  // ShapeFunction::C3D8_dNdr(dNdr,0e0,0e0,0e0);
  calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
  detJ = mathTool::calcDeterminant_3x3(dxdr);
  calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
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
    cout << "centroid Jacobian is nagtive. Exit..." << endl;
    for(int p=0;p<numOfNodeInElm;p++){
      printf("%e %e %e\n",x_current[p][0],x_current[p][1],x_current[p][2]);
    }
    exit(1);
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      C[i][j]=0e0;
      for(int k=0;k<3;k++) C[i][j]+=F[k][i]*F[k][j];
    }
  }

  mathTool::calcInverseMatrix_3x3(invC,C);

  //S_vol or hydrostatic pressure
  pressure=bulkModulus*(2e0*(J-1e0)-r*pow(J,-r-1e0)+r);
  dpressure=bulkModulus*(2e0-r*(-r-1e0)*pow(J,-r-2e0));
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      sigma[i][j] = pressure*I2[i][j];
    }
  }
  //calc_internal force vector
  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) Qu[ic][p][i] += sigma[i][j] * dNdx[p][j] * detJ * gauss0.weight[i1] * gauss0.weight[i2] * gauss0.weight[i3];
      // for(int j=0;j<3;j++) Qu[ic][p][i] += sigma[i][j] * dNdx[p][j] * detJ * 2e0 * 2e0 * 2e0;
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
        for(int l=0;l<3;l++) C4vol[i][j][k][l]=J*(pressure+J*dpressure)*invC[i][j]*invC[k][l]-2e0*J*pressure*invC_odot[i][j][k][l];
      }
    }
  }

  tensorPushForward_4order(elasticityTensor_current,C4vol,F,J);

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++){
          tangentCoefficient[i][j][k][l] = elasticityTensor_current[i][k][j][l]+sigma[k][l]*I2[i][j];
        }
      }
    }
  }

  for(int p=0;p<numOfNodeInElm;p++){
    for(int q=0;q<numOfNodeInElm;q++){
      for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
          for(int k=0;k<3;k++){
            for(int l=0;l<3;l++){
              // K[ic][p][q][i][j] += dNdx[p][k]*(tangentCoefficient[i][j][k][l]*dNdx[q][l]) * detJ * 2e0 * 2e0 * 2e0;
              K[ic][p][q][i][j] += dNdx[p][k]*(tangentCoefficient[i][j][k][l]*dNdx[q][l]) * detJ * gauss0.weight[i1] * gauss0.weight[i2] * gauss0.weight[i3];
            }
          }
        }
      }
    }
  }
  }}}
  //------------

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
            B[i][j]=0e0;
            C[i][j]=0e0;
            for(int k=0;k<3;k++){
              B[i][j]+=F[i][k]*F[j][k];
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
        Ic1_tilda=pow(J,-2e0/3e0)*(C[0][0]+C[1][1]+C[2][2]);

        //S_bar
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) Sbar[i][j]=2e0*c1*exp(alpha1*(Ic1_tilda-3e0))*I2[i][j];
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

        //sigma an-isotropic term
        for(int i=0;i<3;i++) a0[i]=fiberDirection_elm[ic][i];

        Ic4=0e0;
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) Ic4 += a0[i]*C[i][j]*a0[j];
        }
        lambda=sqrt(Ic4);
        for(int i=0;i<3;i++){
          a[i]=0e0;
          for(int j=0;j<3;j++) a[i] += F[i][j] * a0[j];
        }
        for(int i=0;i<3;i++) a[i] = a[i]/lambda;
        averageLambda+=lambda;
        for(int i=0;i<3;i++) lambda_ave[ic][i]+=a[i];

        term4=2e0*(c3/alpha3*exp(alpha3*(Ic4-1e0))-c3*(Ic4-1e0));
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) S_aniso[i][j]=term4*a0[i]*a0[j];
        }

        //transformation
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            sigma[i][j] = 0e0;
            for(int k=0;k<3;k++){
              for(int l=0;l<3;l++) sigma[i][j] += F[i][k]*(Siso[k][l]+S_aniso[k][l])*F[j][l]/J;
            }
          }
        }

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
              for(int l=0;l<3;l++) P4bar[i][j][k][l]=invC_odot[i][j][k][l]-1e0/3e0*(invC[i][j]*invC[k][l]);
            }
          }
        }

        //------------------------------specific routine-----------------------------
        //isotropic term
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
              for(int l=0;l<3;l++) C4bar[i][j][k][l]=4e0*pow(J,-4e0/3e0)*alpha1*c1*exp(alpha1*(Ic1_tilda-3e0))*I2[i][j]*I2[k][l];
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
              }
            }
          }
        }

        //ansiotropic term
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
              for(int l=0;l<3;l++) C4iso[i][j][k][l]+=4e0*c3*(exp(alpha3*(Ic4-1e0))-1e0)*a0[i]*a0[j]*a0[k]*a0[l];
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
  double tmp=sqrt(lambda_ave[ic][0]*lambda_ave[ic][0]+lambda_ave[ic][1]*lambda_ave[ic][1]+lambda_ave[ic][2]*lambda_ave[ic][2]);
  for(int j=0;j<3;j++) lambda_ave[ic][j]/=tmp;
  for(int j=0;j<3;j++) lambda_ave[ic][j]*=averageLambda/8e0;
  // if(averageLambda/8e0>1.1e0) cout << "test" << endl;
  // printf("volume=%e\n",volume);
  Allocation::free2d(x_current);
  Allocation::free2d(x_ref);
  Allocation::free2d(dNdr);
  Allocation::free2d(dNdx);
}
