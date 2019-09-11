/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "fem.h"

using namespace std;

// #################################################################
/**
 * @brief calc surface boundary force
 */
void Fem::calc_surfaceBoundaryForce()
{

  #pragma omp parallel for
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) externalForce(i,j) = 0e0;
  }

  #pragma omp parallel for
  for(int ic=0;ic<numOfNeumann;ic++){

    for(int p=0;p<belement[ic].node.size();p++){
      for(int i=0;i<3;i++) BFe(ic,p,i) = 0e0;
    }

    switch(belement[ic].meshType){
      case VTK_TRIANGLE:
        // calc_Traction_element_triangle(ic,3,3);
        cout << "未実装です" << endl;
        exit(1);
        break;
      case VTK_QUAD:
        calc_Traction_element_quad(ic,4,2);
        break;
      default:
        cout << "undefined boundary mesh type" << endl;
        exit(1);
      }
  }

  for(int ic=0;ic<numOfNeumann;ic++){
    for(int p=0;p<belement[ic].node.size();p++){
      for(int i=0;i<3;i++) externalForce(belement[ic].node[p],i) += BFe(ic,p,i);
    }
  }
}

// #################################################################
/**
 * @brief calc traction force in each surface element
 * @param [in] ic element number
 * @param [in] numOfNodeInBdElm number of node in boundary element
 * @param [in] numOfGaussPoint number of gauss point set in each element
 */
void Fem::calc_Traction_element_quad(const int &ic,const int numOfNodeInBdElm,const int numOfGaussPoint)
{
  double dXdr1[3],dXdr2[3],dXdr3[3],Jacobian,normalVector[3];

  DOUBLEARRAY2D X(numOfNodeInBdElm,3);
  DOUBLEARRAY1D N(numOfNodeInBdElm);
  DOUBLEARRAY2D dNdr(numOfNodeInBdElm,2);
  Gauss gauss(numOfGaussPoint);

  for(int p=0;p<numOfNodeInBdElm;p++){
    for(int i=0;i<3;i++) X(p,i) = x0(belement[ic].node[p],i);
  }

  for(int i1=0;i1<numOfGaussPoint;i1++){
    for(int i2=0;i2<numOfGaussPoint;i2++){

      ShapeFunction2D::C2D4_N(N,gauss.point[i1],gauss.point[i2]);
      ShapeFunction2D::C2D4_dNdr(dNdr,gauss.point[i1],gauss.point[i2]);

      for(int i=0;i<3;i++){
        dXdr1[i] = 0e+0;
        dXdr2[i] = 0e+0;
        for(int p=0;p<numOfNodeInBdElm;p++){
          dXdr1[i] += dNdr(p,0) * X(p,i);
          dXdr2[i] += dNdr(p,1) * X(p,i);
        }
      }

      mathTool::crossProduct(dXdr1,dXdr2,dXdr3,Jacobian);

      for(int i=0;i<3;i++) normalVector[i] = dXdr3[i] / Jacobian;

      for(int i=0;i<3;i++){
        for(int p=0;p<numOfNodeInBdElm;p++){
          BFe(ic,p,i) += N(p) * bn(ic,i) * Jacobian * gauss.weight[i1] * gauss.weight[i2];
        }
      }
    }
  }
}
// #################################################################
/**
 * @brief calc boundary conditions
 * @param [in] stress
 */
void Fem::calc_TractionByPressure_element(const int &ic,const Element &boundaryElement)
{
  int numOfNodeInBdElm=boundaryElement.node.size();
  double dXdr1[3],dXdr2[3],dXdr3[3],Jacobian,normalVector[3];

  DOUBLEARRAY2D dNdr(numOfNodeInBdElm,2);
  DOUBLEARRAY2D X(numOfNodeInBdElm,3);

  for(int p=0;p<numOfNodeInBdElm;p++){
    for(int i=0;i<3;i++) X(p,i) = x(boundaryElement.node[p],i);
  }

  Gauss g(2);
  GaussTriangle gTri(1),gTri2(2);

  switch(boundaryElement.meshType){
    case VTK_TRIANGLE:
      ShapeFunction2D::C2D3_dNdr(dNdr,gTri.point[0][0],gTri.point[0][1],gTri.point[0][2]);
      calcBFe_inGaussIntegral(BFe,dNdr,X,numOfNodeInBdElm,gTri.weight[0],ic);
      break;
    case VTK_QUADRATIC_TRIANGLE:
      for(int i1=0;i1<3;i1++){
        ShapeFunction2D::C2D6_dNdr(dNdr,gTri2.point[i1][0],gTri2.point[i1][1],gTri.point[i1][2]);
        calcBFe_inGaussIntegral(BFe,dNdr,X,numOfNodeInBdElm,gTri2.weight[i1],ic);
      }
      break;
    default:
      cout << "undefined surface element" << endl;
      exit(1);
  }
}


// #################################################################
/**
 * @brief calc BFe (internal force in each element during Gauss integral)
 * @param [in] stress
 */
void Fem::calcBFe_inGaussIntegral(DOUBLEARRAY3D &BFe,DOUBLEARRAY2D &dNdr,
  DOUBLEARRAY2D &X,const int numOfNodeInBdElm,const double weight,const int ic)
{
  double dXdr1[3],dXdr2[3],dXdr3[3],Jacobian,normalVector[3];

    for(int i=0;i<3;i++){
      dXdr1[i] = 0e0;
      dXdr2[i] = 0e0;
      for(int p=0;p<numOfNodeInBdElm;p++){
        dXdr1[i] += dNdr(p,0) * X(p,i);
        dXdr2[i] += dNdr(p,1) * X(p,i);
      }
    }

    mathTool::crossProduct(dXdr1,dXdr2,dXdr3,Jacobian);

    for(int j=0;j<3;j++) normalVector[j] = dXdr3[j]/Jacobian;

    for(int p=0;p<numOfNodeInBdElm;p++){
      for(int j=0;j<3;j++) BFe(ic,p,j) += -boundaryPressure * normalVector[j] * (5e-1* Jacobian) * weight;
    }

}


