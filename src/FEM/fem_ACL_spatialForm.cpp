
/**
 * @file fem_ACL_spatialForm.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "fem.h"

using namespace std;

// #################################################################
/**
 * @brief calc stress tensor of SantVenant material
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] option           true or faluse: calculate tangential stiffness matrix or not.
 */
void Fem::calcStressTensor_ACL_element_spatialForm(const int ic,const DOUBLEARRAY2 &U_tmp,const bool option)
{
  int numOfNodeInElm=element[ic].node.size();
  DOUBLEARRAY2 x_current=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 x_ref=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 dNdr=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);
  DOUBLEARRAY2 dNdx=Allocation::allocate2dDOUBLE(numOfNodeInElm,3);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      x_current[p][i] = x0[element[ic].node[p]][i]+U_tmp[element[ic].node[p]][i];
      x_ref[p][i]     = x0[element[ic].node[p]][i];
    }
  }

  Gauss g(1),g2(2);
  GaussTetra gTet(1),gTet2(2);
  GaussTriangle gTri(1),gTri2(2);

  switch(element[ic].meshType){
    case VTK_TETRA:
      ShapeFunction::C3D4_dNdr(dNdr,gTet.point[0][0],gTet.point[0][1],gTet.point[0][2],gTet.point[0][3]);
      ACL_inGaussIntegral(dNdr,x_current,x_ref,dNdx,numOfNodeInElm,gTet.weight[0]*1e0/6e0,ic,option);
      break;
    case VTK_HEXAHEDRON:
      for(int i1=0;i1<2;i1++){
        for(int i2=0;i2<2;i2++){
          for(int i3=0;i3<2;i3++){
            ShapeFunction::C3D8_dNdr(dNdr,g.point[i1],g.point[i2],g.point[i3]);
            ACL_inGaussIntegral(dNdr,x_current,x_ref,dNdx,numOfNodeInElm,g.weight[i1]*g.weight[i2]*g.weight[i3],ic,option);
          }
        }
      }
      break;
    case VTK_QUADRATIC_TETRA:
      for(int i1=0;i1<4;i1++){
        ShapeFunction::C3D10_dNdr(dNdr,gTet2.point[i1][0],gTet2.point[i1][1],gTet2.point[i1][2],gTet2.point[i1][3]);
        ACL_inGaussIntegral(dNdr,x_current,x_ref,dNdx,numOfNodeInElm,gTet2.weight[i1]*1e0/6e0,ic,option);
      }
      break;
    case VTK_QUADRATIC_HEXAHEDRON:
      break;
    case VTK_WEDGE:
      for(int i1=0;i1<2;i1++){
        ShapeFunction::C3D6_dNdr(dNdr,gTri.point[0][0],gTri.point[0][1],gTri.point[0][2],g.point[i1]);
        ACL_inGaussIntegral(dNdr,x_current,x_ref,dNdx,numOfNodeInElm,5e-1*gTri.weight[0]*g.weight[i1],ic,option);
      }
      break;
    default:
      cout << "undefined mesh type" << endl;
      exit(1);
  }
  Allocation::free2d(x_current);
  Allocation::free2d(x_ref);
  Allocation::free2d(dNdr);
  Allocation::free2d(dNdx);
}

// #################################################################
/**
 * @brief calc stress tensor of ACL with using Fbar method
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] numOfNodeInElm   number of node in each element
 * @param [in] numOfGaussPoint  number of Gauss point set in each element
 * @param [in] option           true or faluse: calculate tangential stiffness matrix or not.
 * @detail
   ACL model: Weiss et al., Comput. Meth. Appl. Mech. Eng., 1996
   ACL model parameter: Pena et al., J. Biomech., 2006
   Fbar method: de Souza Neto E et al., Int. J. Solids and Struct. 1996; 33(20-22):3277-3296
 */
void Fem::ACL_inGaussIntegral(const DOUBLEARRAY2 &dNdr,const DOUBLEARRAY2 &x_current,const DOUBLEARRAY2 &x_ref,
DOUBLEARRAY2 &dNdx,const int numOfNodeInElm,const double weight,const int ic,const bool option)
{
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  double sigma[3][3],S[3][3],C[3][3],invC[3][3],F[3][3],B[3][3];
  double elasticityTensor_ref[3][3][3][3],elasticityTensor_current[3][3][3][3],tangentCoefficient[3][3][3][3];

  double c1=1.95e0;
  double bulkModulus=c1*1e3;

  //nearly-incompressible material
  double pressure,dpressure,Siso[3][3],Saniso[3][3],Sbar[3][3],P4[3][3][3][3],P4bar[3][3][3][3];
  double invC_odot[3][3][3][3],C4iso[3][3][3][3],C4vol[3][3][3][3],C4bar[3][3][3][3],term2;

  //ACL
  int fiberNum=0;
  double sigma_aniso[3][3];
  double lambda,lambda_bar=1.046e0;
  double c3=0.0139e0;
  double c4=116.22e0;
  double c5=535.039e0;
  double c6=c3*(exp(c4*(lambda_bar-1e0))-1e0)-c5*lambda_bar;
  double term4,term4_2,a0[3],a[3];
  double averageLambda=0e0;
  for(int i=0;i<3;i++) lambda_ave[ic][i]=0e0;

  calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
  calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
  detJ = mathTool::calcDeterminant_3x3(dxdr);
  volume += detJ * weight;

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
  //S_vol or hydrostatic pressure
  pressure=bulkModulus*log(J)/J;
  dpressure=bulkModulus*(1e0-log(J))/(J*J);

  //sigma an-isotropic term
  for(int i=0;i<3;i++) a0[i]=fiberDirection[ic][fiberNum][i];
  fiberNum++;

  lambda=0e0;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) lambda += a0[i]*C[i][j]*a0[j];
  }
  lambda=sqrt(lambda);

  if(lambda<1e0){
    term4=0e0;
  }else if(lambda<lambda_bar){
    term4=c3*(exp(c4*(lambda-1e0))-1e0);
  }else{
    term4=c5*lambda+c6;
  }
  for(int i=0;i<3;i++){
    a[i]=0e0;
    for(int j=0;j<3;j++) a[i] += F[i][j] * a0[j];
  }
  for(int i=0;i<3;i++) a[i] = a[i]/lambda;

  averageLambda+=lambda;
  for(int i=0;i<3;i++) lambda_ave[ic][i]+=a[i];

  //S_bar from anisotropic term
  // for(int i=0;i<3;i++){
  //   for(int j=0;j<3;j++) Sbar[i][j]+=term4*a0[i]*a0[j]/(lambda*lambda);
  // }
  //------------------------------end specific routine-----------------------------

        //S_iso
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      Siso[i][j]=0e0;
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) Siso[i][j] += pow(J,-2e0/3e0)*P4[i][j][k][l]*Sbar[k][l];
      }
    }
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      sigma[i][j] = 0e0;
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) sigma[i][j] += F[i][k]*Siso[k][l]*F[j][l]/J;
      }
      sigma[i][j] += term4*a[i]*a[j]/J; //aniso
      sigma[i][j] += pressure*I2[i][j]; //volume
    }
  }

        //calc_internal force vector
  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) Qu[ic][p][i] += sigma[i][j] * dNdx[p][j] * detJ * weight;
    }
  }

  if(option==false) return;

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
    //anisotropic term
  if(lambda<1e0){
    term4_2=0e0;
  }else if(lambda<lambda_bar){
    term4_2=c3*c4*exp(c4*(lambda-1e0));
  }else{
    term4_2=c5;
  }
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) C4bar[i][j][k][l]=pow(J,-4e0/3e0)*(term4_2-2e0*term4/lambda)/(lambda*lambda*lambda)*a0[i]*a0[j]*a0[k]*a0[l];
      }
    }
  }
        //------------------------------end specific routine-----------------------------

  term2=0e0;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) term2+=pow(J,-2e0/3e0)*(Sbar[i][j]*C[i][j]);
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++){
          C4iso[i][j][k][l]=0e0;
          // for(int m=0;m<3;m++){
          //   for(int n=0;n<3;n++){
          //     for(int p=0;p<3;p++){
          //       for(int q=0;q<3;q++) C4iso[i][j][k][l]+=P4[i][j][m][n]*C4bar[m][n][p][q]*P4[k][l][p][q];
          //     }
          //   }
          // }
          C4iso[i][j][k][l] += 2e0/3e0*term2*P4bar[i][j][k][l]-2e0/3e0*(invC[i][j]*Siso[k][l]+Siso[i][j]*invC[k][l]);
        }
      }
    }
  }
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++){
          elasticityTensor_ref[i][j][k][l]=C4iso[i][j][k][l]+C4vol[i][j][k][l];
          elasticityTensor_ref[i][j][k][l]+=(lambda*term4_2-2e0*term4)*a[i]*a[j]*a[k]*a[l]/J;
        }
      }
    }
  }
  tensorPushForward_4order(elasticityTensor_current,elasticityTensor_ref,F,J);

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++){
          tangentCoefficient[i][j][k][l] = elasticityTensor_current[i][k][j][l]+sigma[k][l]*I2[i][j];
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
              K[ic][p][q][i][j] += dNdx[p][k]*tangentCoefficient[i][j][k][l]*dNdx[q][l] * detJ * weight;
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
}
// #################################################################
/**
 * @brief calc stress tensor of ACL with using Fbar method
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] numOfNodeInElm   number of node in each element
 * @param [in] numOfGaussPoint  number of Gauss point set in each element
 * @param [in] option           true or faluse: calculate tangential stiffness matrix or not.
 * @detail
   ACL model: Weiss et al., Comput. Meth. Appl. Mech. Eng., 1996
   ACL model parameter: Pena et al., J. Biomech., 2006
   Fbar method: de Souza Neto E et al., Int. J. Solids and Struct. 1996; 33(20-22):3277-3296
 */
void Fem::calcStressTensor_ACL_element_spatialForm_hexa_Fbar(const int &ic,const DOUBLEARRAY2 &U_tmp,
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

  double c1=1.95e0;
  double bulkModulus=c1*1e3;

  //nearly-incompressible material
  double pressure,dpressure,Siso[3][3],Sbar[3][3],P4[3][3][3][3],P4bar[3][3][3][3];
  double invC_odot[3][3][3][3],C4iso[3][3][3][3],C4vol[3][3][3][3],C4bar[3][3][3][3],term2;

  //ACL
  int fiberNum=0;
  double sigma_aniso[3][3];
  double lambda,lambda_bar=1.046e0;
  double c3=0.0139e0;
  double c4=116.22e0;
  double c5=535.039e0;
  double c6=c3*(exp(c4*(lambda_bar-1e0))-1e0)-c5*lambda_bar;
  double term4,term4_2,a0[3],a[3];
  double averageLambda=0e0;
  for(int i=0;i<3;i++) lambda_ave[ic][i]=0e0;
  double F_initial[3][3],Ftmp[3][3];
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
  if(J0<0e0){
    cout << "centroid Jacobian is nagtive. Exit..." << endl;
    for(int p=0;p<numOfNodeInElm;p++){
      printf("%e %e %e\n",x_current[p][0],x_current[p][1],x_current[p][2]);
    }
    exit(1);
  }
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

        //-------------initial stretch-------------------
        if(bundleElement[ic]==-1){  //AMB
          for(int i=0;i<3;i++) a0[i]=fiberDirection[ic][fiberNum][i];
          calc_F_initial(F_initial,a0,AMBinitialStretch);
          for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
              Ftmp[i][j]=0e0;
              for(int k=0;k<3;k++) Ftmp[i][j]+=F[i][k]*F_initial[k][j];
            }
          }
          for(int i=0;i<3;i++){
            for(int j=0;j<3;j++) F[i][j]=Ftmp[i][j];
          }
        }
        if(bundleElement[ic]==1){ //PLB
          for(int i=0;i<3;i++) a0[i]=fiberDirection[ic][fiberNum][i];
          calc_F_initial(F_initial,a0,PLBinitialStretch);
          for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
              Ftmp[i][j]=0e0;
              for(int k=0;k<3;k++) Ftmp[i][j]+=F[i][k]*F_initial[k][j];
            }
          }
          for(int i=0;i<3;i++){
            for(int j=0;j<3;j++) F[i][j]=Ftmp[i][j];
          }
        }
        //-----------------------------------------------

        J = mathTool::calcDeterminant_3x3(F);

        if(J<0e0){
          cout << "Jacobian is nagtive. Exit..." << endl;
          exit(1);
        }

        //-----------Fbar
        // for(int i=0;i<3;i++){
        //   for(int j=0;j<3;j++) F[i][j] = F[i][j] * pow(J0/J,1e0/3e0);
        // }
        // J = J0;
        //----------------

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
        //S_vol or hydrostatic pressure
        pressure=bulkModulus*log(J)/J;
        dpressure=bulkModulus*(1e0-log(J))/(J*J);

        //sigma an-isotropic term
        for(int i=0;i<3;i++) a0[i]=fiberDirection[ic][fiberNum][i];
        fiberNum++;

        lambda=0e0;
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) lambda += a0[i]*C[i][j]*a0[j];
        }
        lambda=sqrt(lambda);

        if(lambda<1e0){
          term4=0e0;
        }else if(lambda<lambda_bar){
          term4=c3*(exp(c4*(lambda-1e0))-1e0);
        }else{
          term4=c5*lambda+c6;
        }
        for(int i=0;i<3;i++){
          a[i]=0e0;
          for(int j=0;j<3;j++) a[i] += F[i][j] * a0[j];
        }
        for(int i=0;i<3;i++) a[i] = a[i]/lambda;

        averageLambda+=lambda;
        for(int i=0;i<3;i++) lambda_ave[ic][i]+=a[i];

        //S_bar from anisotropic term
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) Sbar[i][j]+=term4*a0[i]*a0[j]/(lambda*lambda);
        }
        //------------------------------end specific routine-----------------------------

        //S_iso
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            Siso[i][j]=0e0;
            for(int k=0;k<3;k++){
              for(int l=0;l<3;l++) Siso[i][j] += pow(J,-2e0/3e0)*P4[i][j][k][l]*Sbar[k][l];
            }
          }
        }

        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            sigma[i][j] = 0e0;
            for(int k=0;k<3;k++){
              for(int l=0;l<3;l++) sigma[i][j] += F[i][k]*Siso[k][l]*F[j][l]/J;
            }
            sigma[i][j] += pressure*I2[i][j];
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
        //anisotropic term
        if(lambda<1e0){
          term4_2=0e0;
        }else if(lambda<lambda_bar){
          term4_2=c3*c4*exp(c4*(lambda-1e0));
        }else{
          term4_2=c5;
        }
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
              for(int l=0;l<3;l++) C4bar[i][j][k][l]=pow(J,-4e0/3e0)*(term4_2-2e0*term4/lambda)/(lambda*lambda*lambda)*a0[i]*a0[j]*a0[k]*a0[l];
            }
          }
        }
        //------------------------------end specific routine-----------------------------

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
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
              for(int l=0;l<3;l++){
                elasticityTensor_ref[i][j][k][l]=C4iso[i][j][k][l]+C4vol[i][j][k][l];
              }
            }
          }
        }
        tensorPushForward_4order(elasticityTensor_current,elasticityTensor_ref,F,J);

        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
              for(int l=0;l<3;l++){
                tangentCoefficient[i][j][k][l] = elasticityTensor_current[i][k][j][l]+sigma[k][l]*I2[i][j];
              }
            }
          }
        }

        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
              for(int l=0;l<3;l++){
                q_Fbar[i][j][k][l]=0e0;
                //--------------------Fbar
                // for(int m=0;m<3;m++){
                //   for(int n=0;n<3;n++){
                //     q_Fbar[i][j][k][l] += 1e0/3e0*tangentCoefficient[i][k][m][n]*I2[m][n]*I2[j][l];
                //   }
                // }
                // q_Fbar[i][j][k][l] -= 2e0/3e0*(sigma[i][k]*I2[j][l]);
                //----------------------------
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
                    K[ic][p][q][i][j] += dNdx[p][k]*(tangentCoefficient[i][j][k][l]*dNdx[q][l]+q_Fbar[i][j][k][l]*(dNdx0[q][l]-dNdx[q][l])) 
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

// #################################################################
/**
 * @brief push forward routine for 4th order tensor
 * @param [out] c4    elasticity tensor in current coordinates
 * @param [in]  C4    elasticity tensor in reference coordinates
 * @param [in]  F     deformation gradient tensor
 * @param [in]  J     Jacobian (volume change ratio)
 */
void Fem::tensorPushForward_4order(double (&c4)[3][3][3][3],const double (&C4)[3][3][3][3],const double (&F)[3][3],const double J)
{
  for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        for(int k=0;k<3;k++){
          for(int l=0;l<3;l++){
            c4[i][j][k][l]=0e0;
            for(int p=0;p<3;p++){
              for(int q=0;q<3;q++){
                for(int r=0;r<3;r++){
                  for(int s=0;s<3;s++) c4[i][j][k][l]+=F[i][p]*F[j][q]*F[k][r]*F[l][s]*C4[p][q][r][s]/J;
                }
            }
          }
        }
      }
    }
  }
}

void Fem::calc_F_initial(double (&F_initial)[3][3],const double (&a0)[3],const double lambda)
{
  double Fhat[3][3],R[3][3],a1[3]={1e0,1e0,1e0},a2[3]={1e0,1e0,1e0},a3[3];
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) Fhat[i][j]=0e0;
  }
  Fhat[0][0]=1e0/sqrt(lambda);
  Fhat[1][1]=1e0/sqrt(lambda);
  Fhat[2][2]=lambda;

  double tmp;
  mathTool::crossProduct(a0,a1,a2,tmp);
  for(int i=0;i<3;i++) a2[i]/=tmp;

  mathTool::crossProduct(a2,a0,a1,tmp);
  for(int i=0;i<3;i++) a1[i]/=tmp;

  for(int i=0;i<3;i++){
    R[i][0]=a1[i];
    R[i][1]=a2[i];
    R[i][2]=a0[i];
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      F_initial[i][j]=0e0;
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++){
          F_initial[i][j]+=R[i][k]*Fhat[k][l]*R[j][l];
        }
      }
    }
  }
}
// #################################################################
/**
 * @brief calc stress tensor of ACL with using Fbar method
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] numOfNodeInElm   number of node in each element
 * @param [in] numOfGaussPoint  number of Gauss point set in each element
 * @param [in] option           true or faluse: calculate tangential stiffness matrix or not.
 * @detail
   ACL model: Weiss et al., Comput. Meth. Appl. Mech. Eng., 1996
   ACL model parameter: Pena et al., J. Biomech., 2006
   Fbar method: de Souza Neto E et al., Int. J. Solids and Struct. 1996; 33(20-22):3277-3296
 */
void Fem::calcStressTensor_ACL_element_spatialForm_hexa_SRI(const int &ic,const DOUBLEARRAY2 &U_tmp,
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

  double c1=1.95e0;
  double bulkModulus=c1*1e3;

  //nearly-incompressible material
  double pressure,dpressure,Siso[3][3],Sbar[3][3],P4[3][3][3][3],P4bar[3][3][3][3];
  double invC_odot[3][3][3][3],C4iso[3][3][3][3],C4vol[3][3][3][3],C4bar[3][3][3][3],term2;

  //ACL
  int fiberNum=0;
  double sigma_aniso[3][3];
  double lambda,lambda_bar=1.046e0;
  double c3=0.0139e0;
  double c4=116.22e0;
  double c5=535.039e0;
  double c6=c3*(exp(c4*(lambda_bar-1e0))-1e0)-c5*lambda_bar;
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
  pressure=bulkModulus*log(J)/J;
  dpressure=bulkModulus*(1e0-log(J))/(J*J);
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

        //-------------initial stretch-------------------
        if(bundleElement[ic]==-1){  //AMB
          for(int i=0;i<3;i++) a0[i]=fiberDirection[ic][fiberNum][i];
          calc_F_initial(F_initial,a0,AMBinitialStretch);
          for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
              Ftmp[i][j]=0e0;
              for(int k=0;k<3;k++) Ftmp[i][j]+=F[i][k]*F_initial[k][j];
            }
          }
          for(int i=0;i<3;i++){
            for(int j=0;j<3;j++) F[i][j]=Ftmp[i][j];
          }
        }
        if(bundleElement[ic]==1){ //PLB
          for(int i=0;i<3;i++) a0[i]=fiberDirection[ic][fiberNum][i];
          calc_F_initial(F_initial,a0,PLBinitialStretch);
          for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
              Ftmp[i][j]=0e0;
              for(int k=0;k<3;k++) Ftmp[i][j]+=F[i][k]*F_initial[k][j];
            }
          }
          for(int i=0;i<3;i++){
            for(int j=0;j<3;j++) F[i][j]=Ftmp[i][j];
          }
        }
        //-----------------------------------------------

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
        for(int i=0;i<3;i++) a0[i]=fiberDirection[ic][fiberNum][i];
        fiberNum++;

        lambda=0e0;
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) lambda += a0[i]*C[i][j]*a0[j];
        }
        lambda=sqrt(lambda);

        if(lambda<1e0){
          term4=0e0;
        }else if(lambda<lambda_bar){
          term4=c3*(exp(c4*(lambda-1e0))-1e0);
        }else{
          term4=c5*lambda+c6;
        }
        for(int i=0;i<3;i++){
          a[i]=0e0;
          for(int j=0;j<3;j++) a[i] += F[i][j] * a0[j];
        }
        for(int i=0;i<3;i++) a[i] = a[i]/lambda;

        averageLambda+=lambda;
        for(int i=0;i<3;i++) lambda_ave[ic][i]+=a[i];

        //S_bar from anisotropic term
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) Sbar[i][j]+=term4*a0[i]*a0[j]/(lambda*lambda);
        }
        //------------------------------end specific routine-----------------------------

        //S_iso
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            Siso[i][j]=0e0;
            for(int k=0;k<3;k++){
              for(int l=0;l<3;l++) Siso[i][j] += pow(J,-2e0/3e0)*P4[i][j][k][l]*Sbar[k][l];
            }
          }
        }

        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            sigma[i][j] = 0e0;
            for(int k=0;k<3;k++){
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
        //anisotropic term
        if(lambda<1e0){
          term4_2=0e0;
        }else if(lambda<lambda_bar){
          term4_2=c3*c4*exp(c4*(lambda-1e0));
        }else{
          term4_2=c5;
        }
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
              for(int l=0;l<3;l++) C4bar[i][j][k][l]=pow(J,-4e0/3e0)*(term4_2-2e0*term4/lambda)/(lambda*lambda*lambda)*a0[i]*a0[j]*a0[k]*a0[l];
            }
          }
        }
        //------------------------------end specific routine-----------------------------

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
