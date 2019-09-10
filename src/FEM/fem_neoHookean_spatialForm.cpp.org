/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "fem.h"

using namespace std;

// #################################################################
/**
 * @brief calc stress tensor of NeoHookean material
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] numOfNodeInElm   number of node in each element
 * @param [in] numOfGaussPoint  number of Gauss point set in each element
 * @param [in] option           true or faluse: calculate tangential stiffness matrix or not.
 */
void Fem::calcStressTensor_NeoHookean_element_spatialForm_hexa(const int &ic,const DOUBLEARRAY2 &U_tmp,
const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option)
{
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  double sigma[3][3],S[3][3],C[3][3],invC[3][3],F[3][3],B[3][3];
  double C4[3][3][3][3],c4[3][3][3][3];

  DOUBLEARRAY2 x_current=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 x_ref=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 dNdr=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 dNdx=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);

  double c1=1.95e0;
  double bulkModulus=c1*1e3;

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
                S[i][j]=2e0*c1*(I2[i][j]-invC[i][j])+bulkModulus*log(J)*invC[i][j];
          }
        }

        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            sigma[i][j] = 0e0;
            for(int k=0;k<3;k++){
              for(int l=0;l<3;l++) sigma[i][j] += F[i][k]*S[k][l]*F[j][l]/J;
            }
          }
        }


        //calc_internal force vector
        for(int p=0;p<numOfNodeInElm;p++){
          for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
              Qu[ic][p][i] += sigma[i][j] * dNdx[p][j] * detJ * gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
            }
          }
        }

        if(option==false) continue;

        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
              for(int l=0;l<3;l++){
                c4[i][j][k][l]=bulkModulus/J*I2[i][j]*I2[k][l]+2e0/J*(2e0*c1-bulkModulus*log(J))*I4[i][j][k][l];
              }
            }
          }
        }

        //calc_tangential_stiffness_matrix
        for(int p=0;p<numOfNodeInElm;p++){
          for(int q=0;q<numOfNodeInElm;q++){
            for(int i=0;i<3;i++){
              for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                  for(int l=0;l<3;l++){
                    K[ic][p][q][i][j] += dNdx[p][k]*(c4[i][k][j][l]+sigma[k][l]*I2[i][j])*dNdx[q][l] 
                    * detJ * gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
                  }
                }
              }
            }
          }
        }

      }
    }
  }
  Allocation::free2d(x_current);
  Allocation::free2d(x_ref);
  Allocation::free2d(dNdr);
  Allocation::free2d(dNdx);
  // printf("volume=%e\n",volume);
}

// #################################################################
/**
 * @brief calc stress tensor of NeoHookean with using Fbar method (de Souza Neto E et al., Int. J. Solids and Struct. 1996; 33(20-22):3277-3296)
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] numOfNodeInElm   number of node in each element
 * @param [in] numOfGaussPoint  number of Gauss point set in each element
 * @param [in] option           true or faluse: calculate tangential stiffness matrix or not.
 */
void Fem::calcStressTensor_NeoHookean_element_spatialForm_hexa_Fbar(const int &ic,const DOUBLEARRAY2 &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option)
{
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  double sigma[3][3],S[3][3],C[3][3],invC[3][3],F[3][3],B[3][3];
  double C4[3][3][3][3],c4[3][3][3][3],tangentCoefficient[3][3][3][3];

  DOUBLEARRAY2 x_current=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 x_ref=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 dNdr=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 dNdx=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);

  double c1=1.95e0;
  double bulkModulus=c1*1e3;

  //------two point---------
  Gauss gauss(numOfGaussPoint);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      x_current[p][i] = x0[element[ic].node[p]][i]+U_tmp[element[ic].node[p]][i];
      x_ref[p][i]     = x0[element[ic].node[p]][i];
    }
  }

  //--- Fbar ---
  DOUBLEARRAY2 dNdx0=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  double q_Fbar[3][3][3][3];
  ShapeFunction::C3D8_dNdr(dNdr,0e0,0e0,0e0);
  calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
  calc_dNdx(dNdx0,dNdr,dxdr,numOfNodeInElm);
  calc_dXdr(dXdr,dNdr,x_ref,numOfNodeInElm);
  mathTool::calcInverseMatrix_3x3(drdX,dXdr);
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      F[i][j]=0e0;
      for(int k=0;k<3;k++) F[i][j] += dxdr[i][k]*drdX[k][j];
    }
  }
  double J0 = mathTool::calcDeterminant_3x3(F);
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

        //Fbar
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) F[i][j] = F[i][j] * pow(J0/J,1e0/3e0);
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
                S[i][j]=2e0*c1*(I2[i][j]-invC[i][j])+bulkModulus*log(J)*invC[i][j];
          }
        }

        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            sigma[i][j] = 0e0;
            for(int k=0;k<3;k++){
              for(int l=0;l<3;l++) sigma[i][j] += F[i][k]*S[k][l]*F[j][l]/J;
            }
          }
        }


        //calc_internal force vector
        for(int p=0;p<numOfNodeInElm;p++){
          for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
              Qu[ic][p][i] += sigma[i][j] * dNdx[p][j] * detJ * gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
            }
          }
        }

        if(option==false) continue;

        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
              for(int l=0;l<3;l++){
                c4[i][j][k][l]=bulkModulus/J*I2[i][j]*I2[k][l]+2e0/J*(2e0*c1-bulkModulus*log(J))*I4[i][j][k][l];
              }
            }
          }
        }

        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
              for(int l=0;l<3;l++){
                tangentCoefficient[i][j][k][l] = c4[i][k][j][l]+sigma[k][l]*I2[i][j];
              }
            }
          }
        }

        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
              for(int l=0;l<3;l++){
                q_Fbar[i][j][k][l]=0e0;
                for(int m=0;m<3;m++){
                  for(int n=0;n<3;n++){
                    q_Fbar[i][j][k][l] += 1e0/3e0*tangentCoefficient[i][k][m][n]*I2[m][n]*I2[j][l];
                  }
                }
                q_Fbar[i][j][k][l] -= 2e0/3e0*(sigma[i][k]*I2[j][l]);
              }
            }
          }
        }

        //calc_tangential_stiffness_matrix
        for(int p=0;p<numOfNodeInElm;p++){
          for(int q=0;q<numOfNodeInElm;q++){
            for(int i=0;i<3;i++){
              for(int j=0;j<3;j++){
                for(int k=0;k<3;k++){
                  for(int l=0;l<3;l++){
                    K[ic][p][q][i][j] += dNdx[p][k]*(tangentCoefficient[i][j][k][l]*dNdx[q][l]+q_Fbar[i][j][k][l]*(dNdx0[q][l]-dNdx[q][l])) * detJ * gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
                  }
                }
              }
            }
          }
        }

      }
    }
  }
  Allocation::free2d(x_current);
  Allocation::free2d(x_ref);
  Allocation::free2d(dNdr);
  Allocation::free2d(dNdx);
}



