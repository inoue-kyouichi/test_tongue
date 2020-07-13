#include "SignedDistanceFunction.h"
using namespace std;

// #################################################################
/**
 * @brief calc prediction matrix
 * @param [in] ic element number
 */
double SignedDistanceFunction::calcEikonalEquation()
{
  VECTOR1D<ARRAY2D<double>> LHS,RHS;
  LHS.resize(FEM.numOfElm);
  RHS.resize(FEM.numOfElm);
  for(int ic=0;ic<FEM.numOfElm;ic++){
    LHS[ic].allocate(4,4);
    RHS[ic].allocate(4,4);
  }

  #pragma omp parallel for
  for(int i=0;i<FEM.numOfElm;i++){
    int numOfNodeInElm = FEM.element[i].node.size();
    for(int j=0;j<numOfNodeInElm;j++){
      for(int k=0;k<numOfNodeInElm;k++){
        LHS[i](j,k) = 0e0;
        RHS[i](j,k) = 0e0;
      }
    }
  }

  #pragma omp parallel for
  for(int ic=0;ic<FEM.numOfElm;ic++){
    calcPredictionMatrix(LHS,RHS,ic);
  }

  ARRAY1D<double> A(FEM.numOfNode);

  #pragma omp parallel for
  for(int ic=0;ic<FEM.numOfNode;ic++) A(ic) = 0e0;

  for(int ic=0;ic<FEM.numOfElm;ic++){
    for(int j=0;j<4;j++){
      for(int k=0;k<4;k++) A(FEM.element[ic].node[j]) += LHS[ic](j,k);
    }
  }

  ARRAY1D<double> b(FEM.numOfNode);

  #pragma omp parallel for
  for(int ic=0;ic<FEM.numOfNode;ic++) b(ic) = 0e0;

  for(int i=0;i<FEM.numOfElm;i++){
    for(int q=0;q<FEM.element[i].node.size();q++){
      for(int j=0;j<FEM.element[i].node.size();j++) b(FEM.element[i].node[q]) += RHS[i](q,j)*SDF(FEM.element[i].node[j]);
    }
  }

  for(int ic=0;ic<FEM.numOfElm;ic++){
    calcSourceTerm(b,ic);
  }

  double residual=0e0;

  #pragma omp parallel for reduction(max:residual)
  for(int ic=0;ic<FEM.numOfNode;ic++){
    if(mask(ic)==0) continue;
    double tmp = fabs(b(ic)/A(ic)-SDF(ic));
    if(tmp>residual) residual = tmp;
    SDF(ic) = b(ic)/A(ic);
  }

  return residual;
}

// #################################################################
/**
 * @brief calc prediction matrix
 * @param [in] ic element number
 */
void SignedDistanceFunction::calcSourceTerm(ARRAY1D<double> &b,const int ic)
{
  int numOfNodeInElm = FEM.element[ic].node.size();

  ARRAY1D<int> e(numOfNodeInElm);
  for(int i=0;i<numOfNodeInElm;i++) e(i) = FEM.element[ic].node[i];

  ARRAY2D<double> x_current(numOfNodeInElm,3);
  ARRAY1D<double> N(numOfNodeInElm);
  ARRAY2D<double> dNdr(numOfNodeInElm,3);
  ARRAY2D<double> dNdx(numOfNodeInElm,3);

  for(int i=0;i<numOfNodeInElm;i++){
    for(int j=0;j<3;j++){
      x_current(i,j) = FEM.x(FEM.element[ic].node[i],j);
    }
  }

  GaussTetra gTet(1);
  ShapeFunction3D::C3D4_N(N,gTet.point[0][0],gTet.point[0][1],gTet.point[0][2],gTet.point[0][3]);
  ShapeFunction3D::C3D4_dNdr(dNdr,gTet.point[0][0],gTet.point[0][1],gTet.point[0][2],gTet.point[0][3]);

  double dxdr[3][3],drdx[3][3];

  FEM_MathTool::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
  FEM_MathTool::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
  double detJ = mathTool::calcDeterminant_3x3(dxdr);

  for(int i=0;i<numOfNodeInElm;i++){
      b(FEM.element[ic].node[i]) += N(i) * detJ * gTet.weight[0]/6e0; //Galerking term
  }


  double advel[3]={0e0,0e0,0e0};
  for(int i=0;i<numOfNodeInElm;i++){
    for(int j=0;j<3;j++) advel[j] += dNdx(i,j)*SDF(e(i));
  }

  double dadv = sqrt(advel[0]*advel[0]+advel[1]*advel[1]+advel[2]*advel[2]);
  for(int j=0;j<3;j++) advel[j] /= (dadv+1e-10);

  double tau = SUPG_stabilizationParameter(dxdr,advel,FEM.element[ic].meshType);

  for(int i=0;i<numOfNodeInElm;i++){
    for(int j=0;j<3;j++){
      // b(FEM.element[ic].node[i]) += tau*advel[j]*dNdx(i,j) * gTet.weight[0]/6e0;   //SUPG term
    }
    if(b(FEM.element[ic].node[i])<0e0) printf("%e %e %e %e\n",b(FEM.element[ic].node[i]),advel[0],advel[1],advel[2]);
  }

}

// #################################################################
/**
 * @brief calc prediction matrix
 * @param [in] ic element number
 */
void SignedDistanceFunction::calc_dt()
{
  double minLength=1e10;
  ARRAY2D<double> dNdr(4,3),x_ref(4,3);
  GaussTetra gTet(1);
  double dxdr[3][3];

  for(int ic=0;ic<FEM.numOfElm;ic++){

    for(int i=0;i<4;i++){
      for(int j=0;j<3;j++) x_ref(i,j) = FEM.x(FEM.element[ic].node[i],j);
    }

    ShapeFunction3D::C3D4_dNdr(dNdr,gTet.point[0][0],gTet.point[0][1],gTet.point[0][2],gTet.point[0][3]);
    FEM_MathTool::calc_dxdr(dxdr,dNdr,x_ref,4);
    double detJ = mathTool::calcDeterminant_3x3(dxdr);
    double volume = detJ * 1e0/6e0;

    double diameter = 2e0 * pow(3e0*volume/(4e0*PI),1e0/3e0);
    if(diameter<minLength) minLength = diameter;
  }

  dt = minLength*1e-1;
}


// #################################################################
/**
 * @brief calc prediction matrix
 * @param [in] ic element number
 */
void SignedDistanceFunction::calcPredictionMatrix(VECTOR1D<ARRAY2D<double>> &LHS,VECTOR1D<ARRAY2D<double>> &RHS,const int ic)
{
  int numOfNodeInElm = FEM.element[ic].node.size();
  ARRAY1D<int> e(numOfNodeInElm);
  ARRAY2D<double> x_current(numOfNodeInElm,3);
  ARRAY1D<double> N(numOfNodeInElm);
  ARRAY2D<double> dNdr(numOfNodeInElm,3);
  ARRAY2D<double> dNdx(numOfNodeInElm,3);

  for(int i=0;i<numOfNodeInElm;i++){
    for(int j=0;j<3;j++){
      x_current(i,j) = FEM.x(FEM.element[ic].node[i],j);
    }
  }

  GaussTetra gTet(1);
  switch(FEM.element[ic].meshType){
    case VTK_TETRA:
      ShapeFunction3D::C3D4_N(N,gTet.point[0][0],gTet.point[0][1],gTet.point[0][2],gTet.point[0][3]);
      ShapeFunction3D::C3D4_dNdr(dNdr,gTet.point[0][0],gTet.point[0][1],gTet.point[0][2],gTet.point[0][3]);
      Prediction_Galerkin_inGaussIntegral(LHS,RHS,N,dNdr,x_current,dNdx,numOfNodeInElm,gTet.weight[0]/6e0,ic);
      Prediction_SUPG_inGaussIntegral(LHS,RHS,N,dNdr,x_current,dNdx,numOfNodeInElm,gTet.weight[0]*1e0/6e0,ic);
      break;
    default:
      cout << "undefined mesh type" << endl;
      exit(1);
  }
}

// #################################################################
/**
 * @brief calc prediction matrix (normal Galerkin term)
 * @param [in] ic element number
 */
void SignedDistanceFunction::Prediction_Galerkin_inGaussIntegral(VECTOR1D<ARRAY2D<double>> &LHS,VECTOR1D<ARRAY2D<double>> &RHS,ARRAY1D<double> &N,ARRAY2D<double> &dNdr,ARRAY2D<double> &x_current,ARRAY2D<double> &dNdx,const int numOfNodeInElm,const double weight,const int ic)
{
  ARRAY1D<int> e(numOfNodeInElm);
  for(int i=0;i<numOfNodeInElm;i++) e(i) = FEM.element[ic].node[i];

  double dxdr[3][3],drdx[3][3];

  FEM_MathTool::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
  FEM_MathTool::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
  double detJ = mathTool::calcDeterminant_3x3(dxdr);
  double volume = detJ * weight;

  ARRAY2D<double> Mass(numOfNodeInElm,numOfNodeInElm);
  for(int i=0;i<numOfNodeInElm;i++){
    for(int j=0;j<numOfNodeInElm;j++){
      Mass(i,j) = N(i)*N(j);
    }
  }

  double advel[3]={0e0,0e0,0e0};
  for(int i=0;i<numOfNodeInElm;i++){
    for(int j=0;j<3;j++) advel[j] += dNdx(i,j)*SDF(e(i));
  }

  double dadv = sqrt(advel[0]*advel[0]+advel[1]*advel[1]+advel[2]*advel[2]);
  for(int j=0;j<3;j++) advel[j] /= (dadv+1e-10);

  ARRAY2D<double> Adv(numOfNodeInElm,numOfNodeInElm);
  for(int i=0;i<numOfNodeInElm;i++){
    for(int q=0;q<numOfNodeInElm;q++){
      Adv(i,q) = 0e0;
      for(int j=0;j<3;j++){
        Adv(i,q) += N(i)*advel[j]*dNdx(q,j);
      }
    }
  }

  for(int i=0;i<numOfNodeInElm;i++){
    for(int j=0;j<numOfNodeInElm;j++){
      LHS[ic](i,j) += (Mass(i,j)/dt) * detJ * weight;
      RHS[ic](i,j) += (Mass(i,j)/dt-Adv(i,j)) * detJ * weight;
    }
  }

}

// #################################################################
/**
 * @brief calc prediction matrix (SUPG term)
 * @param [in] ic element number
 */
void SignedDistanceFunction::Prediction_SUPG_inGaussIntegral(VECTOR1D<ARRAY2D<double>> &LHS,VECTOR1D<ARRAY2D<double>> &RHS,ARRAY1D<double> &N,ARRAY2D<double> &dNdr,ARRAY2D<double> &x_current,ARRAY2D<double> &dNdx,const int numOfNodeInElm,const double weight,const int ic)
{
  ARRAY1D<int> e(numOfNodeInElm);
  for(int i=0;i<numOfNodeInElm;i++) e(i) = FEM.element[ic].node[i];

  double dxdr[3][3],drdx[3][3];

  FEM_MathTool::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
  mathTool::calcInverseMatrix_3x3(drdx,dxdr);
  FEM_MathTool::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
  double detJ = mathTool::calcDeterminant_3x3(dxdr);
  double volume = detJ * weight;

  double advel[3]={0e0,0e0,0e0};
  for(int i=0;i<numOfNodeInElm;i++){
    for(int j=0;j<3;j++) advel[j] += dNdx(i,j)*SDF(e(i));
  }

  double dadv = sqrt(advel[0]*advel[0]+advel[1]*advel[1]+advel[2]*advel[2]);
  for(int j=0;j<3;j++) advel[j] /= (dadv+1e-10);

  double tau = SUPG_stabilizationParameter(dxdr,advel,FEM.element[ic].meshType);

  //SUPG mass
  ARRAY2D<double> Mass(numOfNodeInElm,numOfNodeInElm);
  for(int i=0;i<numOfNodeInElm;i++){
    for(int q=0;q<numOfNodeInElm;q++){
      Mass(i,q) = 0e0;
      for(int j=0;j<3;j++){
        Mass(i,q) += tau*advel[j]*dNdx(i,j)*N(q);
      }
    }
  }

  //SUPG advection term
  ARRAY1D<double> vdNdx(numOfNodeInElm);
  for(int i=0;i<numOfNodeInElm;i++){
    vdNdx(i) = 0e0;
    for(int j=0;j<3;j++){
      vdNdx(i) += advel[j] * dNdx(i,j);
    }
  }

  ARRAY2D<double> Adv(numOfNodeInElm,numOfNodeInElm);
  for(int i=0;i<numOfNodeInElm;i++){
    for(int q=0;q<numOfNodeInElm;q++){
        Adv(i,q) = tau*vdNdx(i)*vdNdx(q);
    }
  }

  for(int i=0;i<numOfNodeInElm;i++){
    for(int j=0;j<numOfNodeInElm;j++){
      LHS[ic](i,j) += (Mass(i,j)/dt) * detJ * weight;
      RHS[ic](i,j) += (Mass(i,j)/dt-Adv(i,j)) * detJ * weight;
    }
  }

}

// #################################################################
/**
 * @brief calc stabilization parameter (SUPG term)
 * @param [in] ic element number
 */
double SignedDistanceFunction::SUPG_stabilizationParameter(const double (&dxdr)[3][3],const double (&advel)[3],VTKCellType &CellType)
{
  double tau = 0e0;
  double term1 = (2e0/dt)*(2e0/dt);

  double drdx[3][3];
  mathTool::calcInverseMatrix_3x3(drdx,dxdr);
  double G[3][3];
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      G[i][j] = 0e0;
      for(int k=0;k<3;k++) G[i][j] += drdx[k][i] * drdx[k][j];
    }
  }

  double term2 = 0e0;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      term2 += advel[i] * G[i][j] * advel[j];
    }
  }

  tau = pow(term1+term2,-5e-1);

  return tau;
}
