/**
 * @file fem_base.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "fem_base_mathTool.h"
#include "math_tools.h"

using namespace std;

// #################################################################
/**
 * @brief calc dudX
 * @param [in] stress
 */
void FEM_MathTool::calc_dudX(double (&dudX)[3][3],ARRAY2D<double> &dNdX,ARRAY2D<double> &u,const int &numOfNodeInElm)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      dudX[i][j] = 0e0;
      for(int p=0;p<numOfNodeInElm;p++) dudX[i][j] += dNdX(p,j) * u(p,i);
    }
  }
}

// #################################################################
/**
 * @brief calc dxdr
 * @param [in] stress
 */
void FEM_MathTool::calc_dxdr(double (&dxdr)[3][3],ARRAY2D<double> &dNdr,ARRAY2D<double> &x1,const int &numOfNodeInElm)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      dxdr[i][j] = 0e0;
      for(int p=0;p<numOfNodeInElm;p++){
        dxdr[i][j] += dNdr(p,j) * x1(p,i);
      }
    }
  }
}
// #################################################################
/**
 * @brief calc dXdr
 * @param [in] stress
 */
void FEM_MathTool::calc_dXdr(double (&dXdr)[3][3],ARRAY2D<double> &dNdr,ARRAY2D<double> &x0,const int &numOfNodeInElm)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      dXdr[i][j] = 0e0;
      for(int p=0;p<numOfNodeInElm;p++){
        dXdr[i][j] += dNdr(p,j) * x0(p,i);
      }
    }
  }
}

// #################################################################
/**
 * @brief calc dNdx
 * @param [in] stress
 */
void FEM_MathTool::calc_dNdx(ARRAY2D<double> &dNdx,ARRAY2D<double> &dNdr,const double (&dxdr)[3][3],const int &numOfNodeInElm)
{
  double drdx[3][3];

  mathTool::calcInverseMatrix_3x3(drdx,dxdr);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      dNdx(p,i) = 0e0;
      for(int j=0;j<3;j++) dNdx(p,i) += dNdr(p,j) * drdx[j][i];
    }
  }
}

// #################################################################
/**
 * @brief calc dNdX
 * @param [in] stress
 */
void FEM_MathTool::calc_dNdX(ARRAY2D<double> &dNdX,ARRAY2D<double> &dNdr,const double (&dXdr)[3][3],const int &numOfNodeInElm)
{
  double drdX[3][3];

  mathTool::calcInverseMatrix_3x3(drdX,dXdr);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      dNdX(p,i) = 0e0;
      for(int j=0;j<3;j++){
        dNdX(p,i) += dNdr(p,j) * drdX[j][i];
      }
    }
  }
}

// #################################################################
/**
 * @brief push forward routine for 4th order tensor
 * @param [out] c4    elasticity tensor in current coordinates
 * @param [in]  C4    elasticity tensor in reference coordinates
 * @param [in]  F     deformation gradient tensor
 * @param [in]  J     Jacobian (volume change ratio)
 */
void FEM_MathTool::tensorPushForward_4order(double (&c4)[3][3][3][3],const double (&C4)[3][3][3][3],const double (&F)[3][3],const double J)
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
