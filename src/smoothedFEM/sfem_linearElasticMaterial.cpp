/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "SFEM.h"

using namespace std;

// #################################################################
/**
 * @brief calc stress tensor of SantVenant material
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] option           true or faluse: calculate tangential stiffness matrix or not.
 */
void SmoothedFEM::SFEM::calcStressTensor_linearElasticMaterial_element_SFEM()
{
  double C4[3][3][3][3];
  double YoungModulus = 1e6; //MPa
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
  double D[6][6];
  constitutiveMatrix(D,C4);

  std::vector<ARRAY2D<double>> B(numOfElm);
  for(int ic=0;ic<numOfElm;ic++) B[ic].allocate(6,12);

  if(element[0].meshType!=VTK_TETRA){
    cout << "Error: please use Tetra element in SFEM routine. Exit..." << endl;
    exit(1);
  }

  int numOfNodeInElm=element[0].node.size();
  ARRAY2D<double> x_ref(numOfNodeInElm,3);
  ARRAY2D<double> dNdr(numOfNodeInElm,3);
  GaussTetra gTet(1);
  ShapeFunction3D::C3D4_dNdr(dNdr,gTet.point[0][0],gTet.point[0][1],gTet.point[0][2],gTet.point[0][3]);

  ARRAY1D<double> volume(numOfElm);

  //calc B matrix
  for(int ic=0;ic<numOfElm;ic++){

    for(int p=0;p<numOfNodeInElm;p++){
      for(int i=0;i<3;i++){
        x_ref(p,i) = x0(element[ic].node[p],i);
      }
    }

    volume(ic) = linearElasticMaterial_SFEM(B[ic],dNdr,x_ref,numOfNodeInElm,gTet.weight[0]/6e0);
  }

  //calc K matrix
  for(int ic=0;ic<numOfElm;ic++){
    
    double K[12][12];
    for(int i=0;i<12;i++){
      for(int j=0;j<12;j++){
        K[i][j] = 0e0;
        for(int k=0;k<6;k++){
          for(int l=0;l<6;l++) K[i][j] += B[ic](k,i) * D[k][l] * B[ic](l,j);
        }
      }
    }

    for(int p=0;p<4;p++){
      for(int q=0;q<4;q++){
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++){
            Ku[ic](p,q,i,j) = K[i+3*p][j+3*q] * volume(ic);
          }
        }
      }
    }

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
double SmoothedFEM::SFEM::linearElasticMaterial_SFEM(ARRAY2D<double> &B,ARRAY2D<double> &dNdr,ARRAY2D<double> &x_ref,const int numOfNodeInElm,const double weight)
{
  ARRAY2D<double> dNdX(numOfNodeInElm,3);

  double dXdr[3][3],drdX[3][3];
  FEM_MathTool::calc_dXdr(dXdr,dNdr,x_ref,numOfNodeInElm);
  FEM_MathTool::calc_dNdX(dNdX,dNdr,dXdr,numOfNodeInElm);
  double detJ = mathTool::calcDeterminant_3x3(dXdr);
  double volume = detJ * weight;

  //Voigt form

  for(int i=0;i<6;i++){
    for(int j=0;j<12;j++) B(i,j) = 0e0;
  }

  for(int p=0;p<4;p++){
    B(0,0+3*p) = dNdX(p,0); 
    B(1,1+3*p) = dNdX(p,1); 
    B(2,2+3*p) = dNdX(p,2);
    B(3,0+3*p) = dNdX(p,1); B(3,1+3*p) = dNdX(p,0);
    B(4,1+3*p) = dNdX(p,2); B(4,2+3*p) = dNdX(p,1);
    B(5,0+3*p) = dNdX(p,2); B(5,2+3*p) = dNdX(p,0);
  }
  return volume;
}

// #################################################################
/**
 * @brief tentative
 */
void SmoothedFEM::SFEM::constitutiveMatrix(double (&D)[6][6],double (&C4)[3][3][3][3])
{
  D[0][0] = C4[0][0][0][0];
  D[0][1] = C4[0][0][1][1];
  D[0][2] = C4[0][0][2][2];
  D[0][3] = C4[0][0][0][1];
  D[0][4] = C4[0][0][0][2];
  D[0][5] = C4[0][0][1][2];

  D[1][1] = C4[1][1][1][1];
  D[1][2] = C4[1][1][2][2];
  D[1][3] = C4[1][1][0][1];
  D[1][4] = C4[1][1][0][2];
  D[1][5] = C4[1][1][1][2];

  D[2][2] = C4[2][2][2][2];
  D[2][3] = C4[2][2][0][1];
  D[2][4] = C4[2][2][0][2];
  D[2][5] = C4[2][2][1][2];

  D[3][3] = C4[0][1][0][1];
  D[3][4] = C4[0][1][0][2];
  D[3][5] = C4[0][1][1][2];

  D[4][4] = C4[0][2][0][2];
  D[4][5] = C4[0][2][1][2];

  D[5][5] = C4[1][2][1][2];

  for(int i=0;i<6;i++){
    for(int j=i+1;j<6;j++) D[j][i] = D[i][j];
  }
}