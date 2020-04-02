/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "fem.h"

using namespace std;


// #################################################################
/**
 * @brief calc external surface boundary force
 */
void Fem::inputSurfaceBoundary(TextParser &tp)
{
  string str,base_label,label,inputDir;
  string N_file,NmeshType_file;

  base_label = "/Domain";
  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }

  label = base_label + "/boundaryElementFile";
  if ( !tp.getInspectedValue(label, N_file)){
    cout << label << " is not found" << endl;
  }
  label = base_label + "/boundaryMeshTypeFile";
  if ( !tp.getInspectedValue(label, NmeshType_file)){
    cout << label << " is not found" << endl;
  }

  N_file = inputDir + "/" + N_file;
  NmeshType_file = inputDir + "/" + NmeshType_file;

  fileIO::read_geometry_meshType(belement,numOfNeumann,NmeshType_file);
  fileIO::read_geometry_element(belement,numOfNeumann,N_file);

  externalSurfaceForce.allocate(numOfNode,3);

  BFe.resize(belement.size());
  for(int i=0;i<belement.size();i++){
    BFe[i].allocate(belement[i].node.size(),3);
  }
}

// #################################################################
/**
 * @brief calc external surface boundary force
 */
void Fem::calc_externalSurfaceForce_prescribedTraction(std::vector<ElementType> &belement,ARRAY2D<double> &Traction)
{
  #pragma omp parallel for
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) externalSurfaceForce(i,j) = 0e0;
  }

  #pragma omp parallel for
  for(int ic=0;ic<belement.size();ic++){

    for(int p=0;p<belement[ic].node.size();p++){
      for(int i=0;i<3;i++) BFe[ic](p,i) = 0e0;
    }

    double tractionForce[3];
    for(int j=0;j<3;j++) tractionForce[j] = Traction(ic,j);

    switch(belement[ic].meshType){
      case VTK_TRIANGLE:
        // calc_Traction_element_triangle(ic,3,1,belement[ic],tractionForce);
        cout << "未実装です" << endl;
        exit(1);
        break;
      case VTK_QUAD:
        calc_Traction_element_quad(ic,4,2,belement[ic],tractionForce);
        break;
      default:
        cout << "undefined boundary mesh type" << endl;
        exit(1);
      }
  }

  for(int ic=0;ic<belement.size();ic++){
    for(int p=0;p<belement[ic].node.size();p++){
      for(int i=0;i<3;i++) externalSurfaceForce(belement[ic].node[p],i) += BFe[ic](p,i);
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
void Fem::calc_Traction_element_quad(const int ic,const int numOfNodeInBdElm,const int numOfGaussPoint,ElementType &belement,const double (&TractionForce)[3])
{
  double dXdr1[3],dXdr2[3],dXdr3[3],Jacobian,normalVector[3];

  ARRAY2D<double> X(numOfNodeInBdElm,3);
  ARRAY1D<double> N(numOfNodeInBdElm);
  ARRAY2D<double> dNdr(numOfNodeInBdElm,2);
  Gauss gauss(1);

  for(int p=0;p<numOfNodeInBdElm;p++){
    for(int i=0;i<3;i++) X(p,i) = x0(belement.node[p],i);
  }

  for(int i1=0;i1<numOfGaussPoint;i1++){
    for(int i2=0;i2<numOfGaussPoint;i2++){

      switch(belement.meshType){
        case VTK_QUAD:
          ShapeFunction2D::C2D4_N(N,gauss.point[i1],gauss.point[i2]);
          ShapeFunction2D::C2D4_dNdr(dNdr,gauss.point[i1],gauss.point[i2]);
          break;
        default:
          cout << "undefined boundary mesh type . Exit..." << endl;
          exit(1);
      }

      for(int i=0;i<3;i++){
        dXdr1[i] = 0e+0;
        dXdr2[i] = 0e+0;
        for(int p=0;p<numOfNodeInBdElm;p++){
          dXdr1[i] += dNdr(p,0) * X(p,i);
          dXdr2[i] += dNdr(p,1) * X(p,i);
        }
      }

      mathTool::crossProduct(dXdr1,dXdr2,dXdr3,Jacobian);
      // for(int i=0;i<3;i++) normalVector[i] = dXdr3[i] / Jacobian;

      for(int i=0;i<3;i++){
        for(int p=0;p<numOfNodeInBdElm;p++){
          BFe[ic](p,i) += N(p) * TractionForce[i] * Jacobian * gauss.weight[i1] * gauss.weight[i2];
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
void Fem::calc_Traction_element_triangle(const int ic,const int numOfNodeInBdElm,const int numOfGaussPoint,ElementType &belement,const double (&TractionForce)[3])
{
  double dXdr1[3],dXdr2[3],dXdr3[3],Jacobian,normalVector[3];

  ARRAY2D<double> dNdr(numOfNodeInBdElm,2);
  ARRAY2D<double> X(numOfNodeInBdElm,3);
  ARRAY1D<double> N(numOfNodeInBdElm);

  for(int p=0;p<numOfNodeInBdElm;p++){
    for(int i=0;i<3;i++) X(p,i) = x0(belement.node[p],i);
  }

  GaussTriangle gTri(numOfGaussPoint);

  for(int i1=0;i1<numOfGaussPoint;i1++){
    switch(belement.meshType){
      case VTK_TRIANGLE:
        ShapeFunction2D::C2D3_N(N,gTri.point[0][0],gTri.point[0][1],gTri.point[0][2]);
        ShapeFunction2D::C2D3_dNdr(dNdr,gTri.point[0][0],gTri.point[0][1],gTri.point[0][2]);
        break;
      // case VTK_QUADRATIC_TRIANGLE:
      //     ShapeFunction2D::C2D6_N(N,gTri.point[i1][0],gTri.point[i1][1],gTri.point[i1][2]);
      //     ShapeFunction2D::C2D6_dNdr(dNdr,gTri.point[i1][0],gTri.point[i1][1],gTri.point[i1][2]);
      //   break;
      default:
        cout << "undefined surface element" << endl;
        exit(1);
    }

    for(int i=0;i<3;i++){
      dXdr1[i] = 0e0;
      dXdr2[i] = 0e0;
      for(int p=0;p<numOfNodeInBdElm;p++){
        dXdr1[i] += dNdr(p,0) * X(p,i);
        dXdr2[i] += dNdr(p,1) * X(p,i);
      }
    }

    mathTool::crossProduct(dXdr1,dXdr2,dXdr3,Jacobian);

    // for(int j=0;j<3;j++) normalVector[j] = dXdr3[j]/Jacobian;

    for(int p=0;p<numOfNodeInBdElm;p++){
      for(int j=0;j<3;j++) BFe[ic](p,j) += N(p) * TractionForce[j] * (5e-1* Jacobian) * gTri.weight[i1];
    }

  }
}