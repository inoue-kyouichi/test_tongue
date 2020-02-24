/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "PDL.h"

using namespace std;

// #################################################################
/**
 * @brief calc stress tensor of SantVenant material
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] option           true or faluse: calculate tangential stiffness matrix or not.
 */
void PeriodontalLigament::calcStressTensor_LinearElastic_element_spatialForm(const int ic,const bool option)
{
  double young = 66.7e-3; //(MPa)
  double poisson = 0.49e0;
  double PDL_lambda = young * poisson / ((1e0+poisson) * (1e0-2e0*poisson));
  double PDL_mu = 5e-1 * young / (1e0+poisson);

  young = 2.1e-3; //(MPa)
  poisson = 0.45e0;
  double PULP_lambda = young * poisson / ((1e0+poisson) * (1e0-2e0*poisson));
  double PULP_mu = 5e-1 * young / (1e0+poisson);

  int numOfNodeInElm=element[ic].node.size();
  DOUBLEARRAY2D x_ref(numOfNodeInElm,3);
  DOUBLEARRAY2D dNdr(numOfNodeInElm,3);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      x_ref(p,i) = x0(element[ic].node[p],i);
    }
  }

  Gauss g(1),g2(2);
  GaussTetra gTet(1),gTet2(2);
  GaussTriangle gTri(1),gTri2(2);

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
      for(int i1=0;i1<4;i1++){
        ShapeFunction3D::C3D10_dNdr(dNdr,gTet2.point[i1][0],gTet2.point[i1][1],gTet2.point[i1][2],gTet2.point[i1][3]);
        if(element[ic].materialType==PDL){
          LinearElastic_inGaussIntegral(dNdr,x_ref,numOfNodeInElm,gTet2.weight[i1]*1e0/6e0,ic,PDL_lambda,PDL_mu,option);
        }else if(element[ic].materialType==PULP){
          LinearElastic_inGaussIntegral(dNdr,x_ref,numOfNodeInElm,gTet2.weight[i1]*1e0/6e0,ic,PULP_lambda,PULP_mu,option);
        }
      }
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
double PeriodontalLigament::LinearElastic_inGaussIntegral(DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &x_ref,
const int numOfNodeInElm,const double weight,const int ic,const double lambda,const double mu,const bool option)
{
  double detJ,volume,J;
  double dXdr[3][3],C4[3][3][3][3];
  DOUBLEARRAY2D dNdx(numOfNodeInElm,3);

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

  if(option==false) return volume;

  //calc_tangential_stiffness_matrix
  for(int p=0;p<numOfNodeInElm;p++){
    for(int q=0;q<numOfNodeInElm;q++){
      for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
          for(int k=0;k<3;k++){
            for(int l=0;l<3;l++){
              K(ic,p,q,i,j) += dNdx(p,k)*C4[i][k][j][l]*dNdx(q,l) * detJ * weight;
            }
          }
        }
      }
    }
  }

  return volume;
}