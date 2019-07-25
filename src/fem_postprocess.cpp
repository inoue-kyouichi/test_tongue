
/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */


#include "fem.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Eigen>

using namespace Eigen;
using namespace std;


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
void Fem::postProcess_PDL_element_spatialForm_hexa_SRI(const int &ic,const DOUBLEARRAY2 &U_tmp,
const int &numOfNodeInElm,const int &numOfGaussPoint)
{
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  double sigma[3][3],S[3][3],C[3][3],invC[3][3],F[3][3],B[3][3],A[3][3],invA[3][3];
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

  double AEigen[3],AEigenVector[3][3],sigmaEigen[3],sigmaEigenVector[3][3];
  for(int i=0;i<3;i++){
    AEigen_Ave[ic][i]=0e0;
    sigmaEigen_Ave[ic][i]=0e0;
    for(int j=0;j<3;j++){
      AEigenVector_Ave[ic][i][j]=0e0;
      sigmaEigenVector_Ave[ic][i][j]=0e0;
    }
  }

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
    }
  }
  }}}

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

        mathTool::calcInverseMatrix_3x3(invA,B);
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) A[i][j]=5e-1*(I2[i][j]-invA[i][j]);
        }
        calcEigen(A,AEigen,AEigenVector);
        calcEigen(sigma,sigmaEigen,sigmaEigenVector);

        for(int i=0;i<3;i++){
          AEigen_Ave[ic][i]+=AEigen[i];
          sigmaEigen_Ave[ic][i]+=sigmaEigen[i];
          // printf("%e %e %e\n",AEigenVector[i][0],AEigenVector[i][1],AEigenVector[i][2]);
          for(int j=0;j<3;j++){
            AEigenVector_Ave[ic][i][j]+=AEigenVector[i][j];
            sigmaEigenVector_Ave[ic][i][j]+=sigmaEigenVector[i][j];
          }
        }

      }
    }
  }
  double tmp=sqrt(lambda_ave[ic][0]*lambda_ave[ic][0]+lambda_ave[ic][1]*lambda_ave[ic][1]+lambda_ave[ic][2]*lambda_ave[ic][2]);
  for(int j=0;j<3;j++) lambda_ave[ic][j]/=tmp;
  for(int j=0;j<3;j++) lambda_ave[ic][j]*=averageLambda/8e0;
  for(int j=0;j<3;j++){
    AEigen_Ave[ic][j]/=8e0;
    sigmaEigen_Ave[ic][j]/=8e0;
  }
  normalize(AEigen_Ave,AEigenVector_Ave,ic);
  normalize(sigmaEigen_Ave,sigmaEigenVector_Ave,ic);
  // if(averageLambda/8e0>1.1e0) cout << "test" << endl;
  // printf("volume=%e\n",volume);
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
void Fem::postProcess_ACL_element_spatialForm_hexa_SRI(const int &ic,const DOUBLEARRAY2 &U_tmp,
const int &numOfNodeInElm,const int &numOfGaussPoint)
{
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  double sigma[3][3],S[3][3],C[3][3],invC[3][3],F[3][3],B[3][3],A[3][3],invA[3][3];
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
  double Ftmp[3][3],F_initial[3][3];

  double averageLambda=0e0;
  for(int i=0;i<3;i++) lambda_ave[ic][i]=0e0;
  double AEigen[3],AEigenVector[3][3],sigmaEigen[3],sigmaEigenVector[3][3];
  for(int i=0;i<3;i++){
    AEigen_Ave[ic][i]=0e0;
    sigmaEigen_Ave[ic][i]=0e0;
    for(int j=0;j<3;j++){
      AEigenVector_Ave[ic][i][j]=0e0;
      sigmaEigenVector_Ave[ic][i][j]=0e0;
    }
  }

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
    }
  }
  }}}

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

        mathTool::calcInverseMatrix_3x3(invA,B);
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) A[i][j]=5e-1*(I2[i][j]-invA[i][j]);
        }
        calcEigen(A,AEigen,AEigenVector);
        calcEigen(sigma,sigmaEigen,sigmaEigenVector);

        for(int i=0;i<3;i++){
          AEigen_Ave[ic][i]+=AEigen[i];
          sigmaEigen_Ave[ic][i]+=sigmaEigen[i];
          // printf("%e %e %e\n",AEigenVector[i][0],AEigenVector[i][1],AEigenVector[i][2]);
          for(int j=0;j<3;j++){
            AEigenVector_Ave[ic][i][j]+=AEigenVector[i][j];
            sigmaEigenVector_Ave[ic][i][j]+=sigmaEigenVector[i][j];
          }
        }

      }
    }
  }
  double tmp=sqrt(lambda_ave[ic][0]*lambda_ave[ic][0]+lambda_ave[ic][1]*lambda_ave[ic][1]+lambda_ave[ic][2]*lambda_ave[ic][2]);
  for(int j=0;j<3;j++) lambda_ave[ic][j]/=tmp;
  for(int j=0;j<3;j++) lambda_ave[ic][j]*=averageLambda/8e0;
  for(int j=0;j<3;j++){
    AEigen_Ave[ic][j]/=8e0;
    sigmaEigen_Ave[ic][j]/=8e0;
  }
  normalize(AEigen_Ave,AEigenVector_Ave,ic);
  normalize(sigmaEigen_Ave,sigmaEigenVector_Ave,ic);
  // if(averageLambda/8e0>1.1e0) cout << "test" << endl;
  // printf("volume=%e\n",volume);
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
void Fem::postProcess_ACL_element_spatialForm_hexa_Fbar(const int &ic,const DOUBLEARRAY2 &U_tmp,
const int &numOfNodeInElm,const int &numOfGaussPoint)
{
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  double sigma[3][3],S[3][3],C[3][3],invC[3][3],F[3][3],B[3][3],A[3][3],invA[3][3];
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
  double Ftmp[3][3],F_initial[3][3];

  double averageLambda=0e0;
  for(int i=0;i<3;i++) lambda_ave[ic][i]=0e0;
  double AEigen[3],AEigenVector[3][3],sigmaEigen[3],sigmaEigenVector[3][3];
  for(int i=0;i<3;i++){
    AEigen_Ave[ic][i]=0e0;
    sigmaEigen_Ave[ic][i]=0e0;
    for(int j=0;j<3;j++){
      AEigenVector_Ave[ic][i][j]=0e0;
      sigmaEigenVector_Ave[ic][i][j]=0e0;
    }
  }

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

        mathTool::calcInverseMatrix_3x3(invA,B);
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) A[i][j]=5e-1*(I2[i][j]-invA[i][j]);
        }
        calcEigen(A,AEigen,AEigenVector);
        calcEigen(sigma,sigmaEigen,sigmaEigenVector);

        for(int i=0;i<3;i++){
          AEigen_Ave[ic][i]+=AEigen[i];
          sigmaEigen_Ave[ic][i]+=sigmaEigen[i];
          // printf("%e %e %e\n",AEigenVector[i][0],AEigenVector[i][1],AEigenVector[i][2]);
          for(int j=0;j<3;j++){
            AEigenVector_Ave[ic][i][j]+=AEigenVector[i][j];
            sigmaEigenVector_Ave[ic][i][j]+=sigmaEigenVector[i][j];
          }
        }

      }
    }
  }
  double tmp=sqrt(lambda_ave[ic][0]*lambda_ave[ic][0]+lambda_ave[ic][1]*lambda_ave[ic][1]+lambda_ave[ic][2]*lambda_ave[ic][2]);
  for(int j=0;j<3;j++) lambda_ave[ic][j]/=tmp;
  for(int j=0;j<3;j++) lambda_ave[ic][j]*=averageLambda/8e0;
  for(int j=0;j<3;j++){
    AEigen_Ave[ic][j]/=8e0;
    sigmaEigen_Ave[ic][j]/=8e0;
  }
  normalize(AEigen_Ave,AEigenVector_Ave,ic);
  normalize(sigmaEigen_Ave,sigmaEigenVector_Ave,ic);
  // if(averageLambda/8e0>1.1e0) cout << "test" << endl;
  // printf("volume=%e\n",volume);
  Allocation::free2d(x_current);
  Allocation::free2d(x_ref);
  Allocation::free2d(dNdr);
  Allocation::free2d(dNdx);
}
// #################################################################
/**
 * @brief normalize vectors
 */
void Fem::normalize(const DOUBLEARRAY2 &AEigen,DOUBLEARRAY3 &AEigenVector_Ave,const int ic)
{
  double tmp=0e0;
  for(int i=0;i<3;i++){
    tmp=     AEigenVector_Ave[ic][i][0]*AEigenVector_Ave[ic][i][0]
            +AEigenVector_Ave[ic][i][1]*AEigenVector_Ave[ic][i][1]
            +AEigenVector_Ave[ic][i][2]*AEigenVector_Ave[ic][i][2];
    if(tmp<1e-15) tmp=1e0;
    tmp=sqrt(tmp);
    for(int j=0;j<3;j++) AEigenVector_Ave[ic][i][j]=AEigenVector_Ave[ic][i][j]/tmp*AEigen[ic][i];
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
