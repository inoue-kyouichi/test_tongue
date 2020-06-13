/**
 * @file fem_PDL_spatialForm2018.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "PDL_analysis.h"

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
void humanPDL::PeriodontalLigament::calcStressTensor_PDL_element_2018(const int &ic,ARRAY2D<double> &U_tmp)
{
  const int numOfNodeInElm = element[ic].node.size();
  double volume=0e0;
  double stress[3][3];

  GaussTetra gTet2(2);
  Gauss g(1),g2(2);
  ARRAY2D<double> x_current(numOfNodeInElm,3);
  ARRAY2D<double> x_ref(numOfNodeInElm,3);
  ARRAY2D<double> dNdr(numOfNodeInElm,3);
  ARRAY2D<double> dNdx(numOfNodeInElm,3);

  double term4,term4_2,a[3];

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      x_current(p,i) = x0(element[ic].node[p],i)+U_tmp(element[ic].node[p],i);
      x_ref(p,i)     = x0(element[ic].node[p],i);
    }
  }

  double F[3][3];

  //--- Selective reduced integration ---
  switch(element[ic].meshType){
    case VTK_HEXAHEDRON:
    for(int i1=0;i1<2;i1++){
      for(int i2=0;i2<2;i2++){
        for(int i3=0;i3<2;i3++){

          double weight = g.weight[i1] * g.weight[i2] * g.weight[i3];
          ShapeFunction3D::C3D8_dNdr(dNdr,g.point[i1],g.point[i2],g.point[i3]);

          double Ic4bar = calcI4bar(F,dNdr,dNdx,x_current,x_ref,numOfNodeInElm,ic);

          if(Ic4bar<1e0){
            calcStressTensor_hyperFoam_element_spatialForm_inGaussIntegral(ic,U,numOfNodeInElm,x_current,x_ref,dNdr,dNdx,weight,stress,true);
          }else{
            calcStressTensor_PDL_element_spatialForm_2018_inGaussIntegral(ic,U,numOfNodeInElm,x_current,x_ref,dNdr,dNdx,weight,stress,true);
          }

        }
      }
    }
    break;
    case VTK_QUADRATIC_TETRA:
    for(int i1=0;i1<4;i1++){
      double weight = gTet2.weight[i1] * 1e0/6e0;
      ShapeFunction3D::C3D10_dNdr(dNdr,gTet2.point[i1][0],gTet2.point[i1][1],gTet2.point[i1][2],gTet2.point[i1][3]);

      double Ic4bar = calcI4bar(F,dNdr,dNdx,x_current,x_ref,numOfNodeInElm,ic);

      if(Ic4bar<1e0){
        calcStressTensor_hyperFoam_element_spatialForm_inGaussIntegral(ic,U,numOfNodeInElm,x_current,x_ref,dNdr,dNdx,weight,stress,true);
      }else{
        calcStressTensor_PDL_element_spatialForm_2018_inGaussIntegral(ic,U,numOfNodeInElm,x_current,x_ref,dNdr,dNdx,weight,stress,true);
      }
    }
    break;
    default:
      cout << "undefine mesh type. Exit..." << endl;
      exit(1);
  }


}

// #################################################################
/**
 * @brief calc I4 bar
 * @param [in] F               deformation gradient tensor
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] numOfNodeInElm   number of node in each element
 * @param [in] numOfGaussPoint  number of Gauss point set in each element
 * @detail
   PDL model and parameters: Ortun-Terrazas et al., J. Mech. Behavior Biomed. Mat., 2018
 */
double humanPDL::PeriodontalLigament::calcI4bar(double (&F)[3][3],ARRAY2D<double> &dNdr,ARRAY2D<double> &dNdx,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,const int numOfNodeInElm,const int ic)
{
  double dxdr[3][3],dXdr[3][3],drdX[3][3];

        FEM_MathTool::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
        double detJ = mathTool::calcDeterminant_3x3(dxdr);
        FEM_MathTool::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
        FEM_MathTool::calc_dXdr(dXdr,dNdr,x_ref,numOfNodeInElm);
        mathTool::calcInverseMatrix_3x3(drdX,dXdr);

        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            F[i][j]=0e0;
            for(int k=0;k<3;k++) F[i][j] += dxdr[i][k]*drdX[k][j];
          }
        }
        double J = mathTool::calcDeterminant_3x3(F);

        if(J<0e0){
          printf("Jacobian is nagtive (element number=%d). Exit...\n",ic);
          for(int p=0;p<numOfNodeInElm;p++){
            printf("%e %e %e\n",x_current(p,0),x_current(p,1),x_current(p,2));
          }
          exit(1);
        }

        double C[3][3];
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            C[i][j]=0e0;
            for(int k=0;k<3;k++) C[i][j]+=F[k][i]*F[k][j];
          }
        }

        //sigma an-isotropic term
        double a0[3];
        for(int i=0;i<3;i++) a0[i]=fiberDirection_elm(ic,i);

        double Ic4bar=0e0;
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) Ic4bar += a0[i]*C[i][j]*a0[j]*pow(J,-2e0/3e0);
        }

  return Ic4bar;
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
void humanPDL::PeriodontalLigament::postProcess_PDL_element_2018(const int &ic,ARRAY2D<double> &U_tmp,const int &numOfNodeInElm,const int &numOfGaussPoint)
{
  double stress[3][3];

  Gauss g(1);
  GaussTetra gTet2(2);
  ARRAY2D<double> x_current(numOfNodeInElm,3);
  ARRAY2D<double> x_ref(numOfNodeInElm,3);
  ARRAY2D<double> dNdr(numOfNodeInElm,3);
  ARRAY2D<double> dNdx(numOfNodeInElm,3);

  //PDL
  double Ic4bar;
  double lambda,averageLambda=0e0,tmp;
  // int fiberNum=0;
  double term4,term4_2,a0[3],a[3];

  for(int j=0;j<3;j++) lambda_ave(ic,j)=0e0;

  //post process
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

  double F[3][3];
  //--- Selective reduced integration ---
  switch(element[ic].meshType){
    case VTK_HEXAHEDRON:
    for(int i1=0;i1<2;i1++){
      for(int i2=0;i2<2;i2++){
        for(int i3=0;i3<2;i3++){

          double weight = g.weight[i1] * g.weight[i2] * g.weight[i3];
          ShapeFunction3D::C3D8_dNdr(dNdr,g.point[i1],g.point[i2],g.point[i3]);
          double Ic4bar = calcI4bar(F,dNdr,dNdx,x_current,x_ref,numOfNodeInElm,ic);

        if(Ic4bar<1e0){
          calcStressTensor_hyperFoam_element_spatialForm_inGaussIntegral(ic,U,numOfNodeInElm,x_current,x_ref,dNdr,dNdx,weight,stress,false);
        }else{
          calcStressTensor_PDL_element_spatialForm_2018_inGaussIntegral(ic,U,numOfNodeInElm,x_current,x_ref,dNdr,dNdx,weight,stress,false);
        }
          summation_postProcess(averageLambda,stress,F,Ic4bar,ic);
        }
      }
    }
    break;
    case VTK_QUADRATIC_TETRA:
    for(int i1=0;i1<4;i1++){
      double weight = gTet2.weight[i1] * 1e0/6e0;
      ShapeFunction3D::C3D10_dNdr(dNdr,gTet2.point[i1][0],gTet2.point[i1][1],gTet2.point[i1][2],gTet2.point[i1][3]);

      double Ic4bar = calcI4bar(F,dNdr,dNdx,x_current,x_ref,numOfNodeInElm,ic);

      if(Ic4bar<1e0){
        calcStressTensor_hyperFoam_element_spatialForm_inGaussIntegral(ic,U,numOfNodeInElm,x_current,x_ref,dNdr,dNdx,weight,stress,false);
      }else{
        calcStressTensor_PDL_element_spatialForm_2018_inGaussIntegral(ic,U,numOfNodeInElm,x_current,x_ref,dNdr,dNdx,weight,stress,false);
      }
      summation_postProcess(averageLambda,stress,F,Ic4bar,ic);
    }
    break;
    default:
      cout << "undefined mesh type. Exit..."  << endl;
    exit(1);
  }

  double lambda_direction[3];

  tmp=sqrt(lambda_ave(ic,0)*lambda_ave(ic,0)+lambda_ave(ic,1)*lambda_ave(ic,1)+lambda_ave(ic,2)*lambda_ave(ic,2));
  for(int j=0;j<3;j++) lambda_ave(ic,j)/=tmp;
  for(int j=0;j<3;j++) lambda_direction[j] = lambda_ave(ic,j);
  for(int j=0;j<3;j++) lambda_ave(ic,j)*=averageLambda/8e0;

  double force[3];
  for(int i=0;i<3;i++){
    force[i] = 0e0;
    for(int j=0;j<3;j++){
      force[i] += stress[i][j] * lambda_direction[j];
    }
  }

  fibreStress(ic) = sqrt(force[0]*force[0]+force[1]*force[1]+force[2]*force[2]);


  for(int j=0;j<3;j++) sigmaEigen_Ave(ic,j)/=8e0;
  normalize(sigmaEigenVector_Ave,ic);

  double stress1st_direction[3];
  for(int j=0;j<3;j++) stress1st_direction[j] = sigmaEigenVector_Ave(ic,0,j);
  
  double angle=0e0;
  for(int j=0;j<3;j++ ) angle += stress1st_direction[j] * lambda_direction[j];
  angle = acos(fabs(angle))* 180e0/PI;
  angleVariation(ic) = angle;

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
void humanPDL::PeriodontalLigament::summation_postProcess(double &averageLambda,const double (&stress)[3][3],const double (&F)[3][3],const double Ic4bar,const int ic)
{
  double sigmaEigen[3],sigmaEigenVector[3][3];
  calcEigen(stress,sigmaEigen,sigmaEigenVector);
  for(int i=0;i<3;i++){
    sigmaEigen_Ave(ic,i)+=sigmaEigen[i];
    for(int j=0;j<3;j++){
      sigmaEigenVector_Ave(ic,i,j)+=sigmaEigenVector[i][j];
    }
  }

  //sigma an-isotropic term
  double a0[3],a[3];
  for(int i=0;i<3;i++) a0[i]=fiberDirection_elm(ic,i);
  double lambda=sqrt(Ic4bar);
  averageLambda+=lambda;
  for(int i=0;i<3;i++){
    a[i]=0e0;
    for(int j=0;j<3;j++) a[i] += F[i][j] * a0[j];
  }
  double tmp=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
  for(int j=0;j<3;j++) lambda_ave(ic,j)+=sqrt(Ic4bar)*a[j]/tmp;
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
void humanPDL::PeriodontalLigament::calcStressTensor_PDL_element_spatialForm_2018_inGaussIntegral(const int &ic,ARRAY2D<double> &U_tmp,
const int &numOfNodeInElm,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,ARRAY2D<double> &dNdr,ARRAY2D<double> &dNdx,const double weight,double (&stress)[3][3],const bool mainLoop)
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
  FEM_MathTool::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
  detJ = mathTool::calcDeterminant_3x3(dxdr);
  volume += detJ * weight;

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
      for(int j=0;j<3;j++) Qu[ic](p,i) += sigma[i][j] * dNdx(p,j) * detJ * weight;
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
              Ku[ic](p,q,i,j) += dNdx(p,k)*(tangentCoefficient[i][j][k][l]*dNdx(q,l)) * detJ * weight;
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
void humanPDL::PeriodontalLigament::calcStressTensor_hyperFoam_element_spatialForm_inGaussIntegral(const int &ic,ARRAY2D<double> &U_tmp,
const int &numOfNodeInElm,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,ARRAY2D<double> &dNdr,ARRAY2D<double> &dNdx,const double weight,double (&stress)[3][3],const bool mainLoop)
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
  FEM_MathTool::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
  detJ = mathTool::calcDeterminant_3x3(dxdr);
  // volume += detJ * weight;

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
      for(int j=0;j<3;j++) Qu[ic](p,i) += sigma[i][j] * dNdx(p,j) * detJ * weight;
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
              Ku[ic](p,q,i,j) += dNdx(p,k)*(tangentCoefficient[i][j][k][l]*dNdx(q,l)) * detJ * weight;
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
void humanPDL::PeriodontalLigament::calcLambda(double (&stretch)[3],double (&stretchDirection)[3][3],const double (&C)[3][3])
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
