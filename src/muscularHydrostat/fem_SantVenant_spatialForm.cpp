/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "eye_mechanics.h"

using namespace std;

// #################################################################
/**
 * @brief calc stress tensor of SantVenant material
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] option           true or faluse: calculate tangential stiffness matrix or not.
 */
void EYE_Mechanics::EyeMech::calcStressTensor_SantVenant_element_spatialForm(const int ic,ARRAY2D<double> &U_tmp,const bool option)
{
  int numOfNodeInElm=element[ic].node.size();
  ARRAY2D<double> x_current(numOfNodeInElm,3);
  ARRAY2D<double> x_ref(numOfNodeInElm,3);
  ARRAY2D<double> dNdr(numOfNodeInElm,3);
  ARRAY2D<double> dNdx(numOfNodeInElm,3);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      x_current(p,i) = x0(element[ic].node[p],i)+U_tmp(element[ic].node[p],i);
      x_ref(p,i)     = x0(element[ic].node[p],i);
    }
  }

  Gauss g(1),g2(2);
  GaussTetra gTet(1),gTet2(2);
  GaussTriangle gTri(1),gTri2(2);

  switch(element[ic].meshType){
    case VTK_TETRA:
      ShapeFunction3D::C3D4_dNdr(dNdr,gTet.point[0][0],gTet.point[0][1],gTet.point[0][2],gTet.point[0][3]);
      SantVenant_inGaussIntegral(dNdr,x_current,x_ref,dNdx,numOfNodeInElm,gTet.weight[0]*1e0/6e0,ic,option);
      break;
    case VTK_HEXAHEDRON:
      for(int i1=0;i1<2;i1++){
        for(int i2=0;i2<2;i2++){
          for(int i3=0;i3<2;i3++){
            ShapeFunction3D::C3D8_dNdr(dNdr,g.point[i1],g.point[i2],g.point[i3]);
            SantVenant_inGaussIntegral(dNdr,x_current,x_ref,dNdx,numOfNodeInElm,g.weight[i1]*g.weight[i2]*g.weight[i3],ic,option);
          }
        }
      }
      break;
    case VTK_QUADRATIC_TETRA:
      for(int i1=0;i1<4;i1++){
        ShapeFunction3D::C3D10_dNdr(dNdr,gTet2.point[i1][0],gTet2.point[i1][1],gTet2.point[i1][2],gTet2.point[i1][3]);
        SantVenant_inGaussIntegral(dNdr,x_current,x_ref,dNdx,numOfNodeInElm,gTet2.weight[i1]*1e0/6e0,ic,option);
      }
      break;
    case VTK_QUADRATIC_HEXAHEDRON:
      break;
    case VTK_WEDGE:
      for(int i1=0;i1<2;i1++){
        ShapeFunction3D::C3D6_dNdr(dNdr,gTri.point[0][0],gTri.point[0][1],gTri.point[0][2],g.point[i1]);
        SantVenant_inGaussIntegral(dNdr,x_current,x_ref,dNdx,numOfNodeInElm,5e-1*gTri.weight[0]*g.weight[i1],ic,option);
      }
      break;
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
double EYE_Mechanics::EyeMech::SantVenant_inGaussIntegral(ARRAY2D<double> &dNdr,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,
ARRAY2D<double> &dNdx,const int numOfNodeInElm,const double weight,const int ic,const bool option)
{
  double C4[3][3][3][3],c4[3][3][3][3];
  double YoungModulus = 1e3; //MPa
  double PoissonRatio = 0.3e0; //[-]
  double lambda=YoungModulus * PoissonRatio/((1e0+PoissonRatio)*(1e0-2e0*PoissonRatio));
  double mu=5e-1*YoungModulus/(1e0+PoissonRatio);

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) C4[i][j][k][l] = lambda*I2[i][j]*I2[k][l] + mu*2e0*I4[i][j][k][l];
      }
    }
  }

  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  FEM_MathTool::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
  FEM_MathTool::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
  double detJ = mathTool::calcDeterminant_3x3(dxdr);
  double volume = detJ * weight;

  FEM_MathTool::calc_dXdr(dXdr,dNdr,x_ref,numOfNodeInElm);
  mathTool::calcInverseMatrix_3x3(drdX,dXdr);

  double F[3][3];
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      F[i][j]=0e0;
      for(int k=0;k<3;k++) F[i][j] += dxdr[i][k]*drdX[k][j];
    }
  }
  double J = mathTool::calcDeterminant_3x3(F);

  double sigma[3][3],S[3][3],E[3][3];
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      E[i][j]=0e0;
      for(int k=0;k<3;k++) E[i][j]+=5e-1*F[k][i]*F[k][j];
    }
    E[i][i]-=5e-1;
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      S[i][j]=0e0;
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) S[i][j]+=C4[i][j][k][l]*E[k][l];
      }
    }
  }
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      sigma[i][j]=0e0;
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) sigma[i][j]+=F[i][k]*S[k][l]*F[j][l]/J;
      }
    }
  }

        //calc_internal force vector
  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        Qu[ic](p,i) += sigma[i][j] * dNdx(p,j) * detJ * weight;
      }
    }
  }

  if(option==false) return volume;

  FEM_MathTool::tensorPushForward_4order(c4,C4,F,J);

  //calc_tangential_stiffness_matrix
  for(int p=0;p<numOfNodeInElm;p++){
    for(int q=0;q<numOfNodeInElm;q++){
      for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
          for(int k=0;k<3;k++){
            for(int l=0;l<3;l++){
              Ku[ic](p,q,i,j) += dNdx(p,k)*(c4[i][k][j][l]+sigma[k][l]*I2[i][j])*dNdx(q,l) * detJ * weight;
            }
          }
        }
      }
    }
  }
  return volume;
}
