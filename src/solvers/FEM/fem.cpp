/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "fem.h"
#include <ostream>
#include <fstream>

using namespace std;


// #################################################################
/**
 * @brief fem solid analysis routine
 * @param [inout] PARDISO   PARDISO class
 */
// void Fem::exportRestartData(const int loop)
// {
//   FILE *fp;
//   string output = "Restart_"+to_string(dataNumber)+"/U_" + to_string(loop) + ".dat";
//   if ((fp = fopen(output.c_str(), "w")) == NULL) {
//     cout << "file open error" << endl;
//     exit(1);
//   }
//   for(int i=0;i<numOfNode;i++){
//     fprintf(fp,"%e %e %e\n",U(i,0),U(i,1),U(i,2));
//   }
//   fclose(fp);
// }

// #################################################################
/**
* @brief set Inner force
 */
void Fem::setInnerForce()
{
  for(int ic=0;ic<numOfElm;ic++){
    for(int p=0;p<element[ic].node.size();p++){
      for(int i=0;i<3;i++){
        innerForce(element[ic].node[p],i) += Qu[ic](p,i);
      }
    }
  }
}

// #################################################################
/**
 * @brief assembly right-hand side
 */
void Fem::set_rhs_statics()
{
 // #pragma omp parallel for
  for(int ic=0;ic<numOfNode;ic++){
    for(int j=0;j<3;j++){
      RHS(ic,j) = (double)ibd(ic,j)*( - innerForce(ic,j));
    }
  }
}

// #################################################################
/**
 * @brief initialize stress tensors
 */
void Fem::stress_tensor_initialize()
{

  #pragma omp parallel for
  for(int ic=0;ic<numOfElm;ic++){
    for(int p=0;p<element[ic].node.size();p++){
      for(int q=0;q<element[ic].node.size();q++){
        for(int i=0;i<3;i++){
          for(int j=0;j<3;j++) K[ic](p,q,i,j) = 0e0;
        }
      }
      for(int i=0;i<3;i++) Qu[ic](p,i) = 0e0;
    }
  }

  #pragma omp parallel for
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) innerForce(i,j) = 0e0;
  }
}

// #################################################################
/**
 * @brief corrector scheme. Normally, line seach method should be used rather than this routine.
 * @param [in] u          displacement vector
 * @param [in] relaxation  relaxation parameters
 */
void Fem::corrector_statistics(const double *u,const double relaxation)
{
  #pragma omp parallel for
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) U(i,j) += u[i+j*numOfNode]*relaxation;
  }
}

// #################################################################
/**
 * @brief calc mass matrix (not completed)
 * @param [in] stress
 */
void Fem::calc_MassMatrix()
{
  int numOfNodeInElm;
  double detJ,dXdr[3][3];

  DOUBLEARRAY1D N;
  DOUBLEARRAY2D dNdr,X;

  for(int ic=0;ic<numOfElm;ic++){

    int numOfGaussPoint = element[ic].numOfGaussPoint;

    Gauss gauss(numOfGaussPoint);

    numOfNodeInElm=element[ic].node.size();
    N.allocate(numOfNodeInElm);
    dNdr.allocate(numOfNodeInElm,3);
    X.allocate(numOfNodeInElm,3);

    for(int p=0;p<numOfNodeInElm;p++){
      for(int i=0;i<3;i++) X(p,i) = x0(element[ic].node[p],i);
    }

    for(int p=0;p<numOfNodeInElm;p++){
      for(int q=0;q<numOfNodeInElm;q++) Mass[ic](p,q) = 0e0;
     }

    for(int i1=0;i1<numOfGaussPoint;i1++){
      for(int i2=0;i2<numOfGaussPoint;i2++){
        for(int i3=0;i3<numOfGaussPoint;i3++){

          switch(element[ic].meshType){
            case VTK_HEXAHEDRON:
              ShapeFunction3D::C3D8_N(N,gauss.point[i1],gauss.point[i2],gauss.point[i3]);
              ShapeFunction3D::C3D8_dNdr(dNdr,gauss.point[i1],gauss.point[i2],gauss.point[i3]);
              break;
            case VTK_TRIQUADRATIC_HEXAHEDRON:
              ShapeFunction3D::C3D27_N(N,gauss.point[i1],gauss.point[i2],gauss.point[i3]);
              ShapeFunction3D::C3D27_dNdr(dNdr,gauss.point[i1],gauss.point[i2],gauss.point[i3]);
              break;
            default:
              cout << "error in calcMassMatrix" << endl;
          }

          FEM_MathTool::calc_dXdr(dXdr,dNdr,X,numOfNodeInElm);
          detJ = mathTool::calcDeterminant_3x3(dXdr);

          //calc_internal force vector
          for(int p=0;p<numOfNodeInElm;p++){
            for(int q=0;q<numOfNodeInElm;q++){
              Mass[ic](p,q) += rho * N(p) * N(q) * detJ * gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
            }
          }

        }
      }
    }
  }

}

// #################################################################
/**
 * @brief calc volume hexahedron
 * @param [in]  ic element number
 * @param [out]  elementVolume volume in each element
 * @param [in]  numOfNodeInElm number of node in each element
 * @param [in]  numOfGaussPoint number of gauss point set in each element
 * @param [in]  option true:current configuration, false:reference configuration
 */
void Fem::calcVolume_hexa(const int &ic,DOUBLEARRAY1D &elementVolume,const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option)
{
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];

  DOUBLEARRAY2D x_current(numOfNodeInElm,3);
  DOUBLEARRAY2D x_ref(numOfNodeInElm,3);
  DOUBLEARRAY2D dNdr(numOfNodeInElm,3);
  DOUBLEARRAY2D dNdx(numOfNodeInElm,3);

  //------two point---------
  Gauss gauss(numOfGaussPoint);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      x_current(p,i) = x0(element[ic].node[p],i)+U(element[ic].node[p],i);
      x_ref(p,i)     = x0(element[ic].node[p],i);
    }
  }

  for(int i1=0;i1<numOfGaussPoint;i1++){
    for(int i2=0;i2<numOfGaussPoint;i2++){
      for(int i3=0;i3<numOfGaussPoint;i3++){

        switch(element[ic].meshType){
          case VTK_HEXAHEDRON:
            ShapeFunction3D::C3D8_dNdr(dNdr,gauss.point[i1],gauss.point[i2],gauss.point[i3]);
            break;
          case VTK_TRIQUADRATIC_HEXAHEDRON:
            ShapeFunction3D::C3D27_dNdr(dNdr,gauss.point[i1],gauss.point[i2],gauss.point[i3]);
            break;
          default:
            cout << "error in calcVolume_hexa" << endl;
            cout << element[ic].meshType << endl;
            break;
        }
        if(option==1){
          FEM_MathTool::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
        }else if(option==0){
          FEM_MathTool::calc_dXdr(dxdr,dNdr,x_ref,numOfNodeInElm);
        }
        FEM_MathTool::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
        detJ = mathTool::calcDeterminant_3x3(dxdr);

        volume += detJ * gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
      }
    }
  }
  elementVolume(ic)=volume;
  // printf("volume=%e\n",volume);
}
