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
void Fem::femSolidAnalysis(PARDISO_solver &PARDISO,RigidBody &RBdy)
{
  int output_iter=1;
  int increment=maxIteration;
  string output;

  FILE *fp;
  string outputFile_tooth = outputDir + "/toothDisplacement.dat";
  if ((fp = fopen(outputFile_tooth.c_str(), "w")) == NULL) {
    exit(1);
  }
  fprintf(fp,"loop F Ux Uy Uz Qx Qy Qz Dispx Dispy Dispz\n");
  fclose(fp);

  for(int j=0;j<3;j++) initialMomentArm[j] = FUpoint[j]-RBdy.xg[j];
  printf("%e %e %e\n",initialMomentArm[0],initialMomentArm[1],initialMomentArm[2]);

  for(int ic=0;ic<numOfElm;ic++) calcVolume_hexa(ic,volume0,8,2,0);

  for(int loop=1;loop<=maxIteration;loop++){

    for(int j=0;j<3;j++) FU[j] = (double)loop/(double)maxIteration * FU_input[j];

    if(NRscheme(PARDISO,RBdy)==1){
      printf("loop=%d\n",loop);
      exit(1);
    }

    for(int i=0;i<numOfNode;i++){
      for(int j=0;j<3;j++) x(i,j) = x0(i,j) + U(i,j);
    }
    for(int i=0;i<numOfElm;i++){
      calcVolume_hexa(i,volume,8,2,1);
      volumeChangeRatio(i)=volume(i)/volume0(i);
    }

    // exportRestartData(loop);
    for(int ic=0;ic<numOfElm;ic++) postProcess_PDL_element_spatialForm_hexa_SRI(ic,U,8,2);

    // exportRestartData(loop);

    // if(loop%output_iter==0){
    //   output = outputDir+"/PDL_"+to_string(loop)+".vtu";
    //   fileIO::export_vtu(x,element,numOfNode,numOfElm,U,volumeChangeRatio,lambda_ave,sigmaEigen_Ave,AEigen_Ave,sigmaEigenVector_Ave,AEigenVector_Ave,innerForce,output);
    // }

    RBdy.updateShape();
    output = outputDir + "/tooth_"+to_string(dataNumber)+"_"+to_string(loop)+".ply";
    RBdy.exportPLY(output);

    double angle[3];
    calc_thetaFromRotationMatrix(angle,RBdy.R);

    double momentArm[3];
    for(int i=0;i<3;i++){
      momentArm[i]=0e0;
      for(int j=0;j<3;j++) momentArm[i] += RBdy.R[i][j] * initialMomentArm[j];
      momentArm[i]-=initialMomentArm[i];
    }

    if ((fp = fopen(outputFile_tooth.c_str(), "a")) == NULL) {
      exit(1);
    }
    double FroceMagnitude= (double)loop/(double)maxIteration * sqrt(FU_input[0]*FU_input[0]+FU_input[1]*FU_input[1]+FU_input[2]*FU_input[2]);
    fprintf(fp,"%d %e %e %e %e %e %e %e %e %e %e\n",loop,FroceMagnitude,RBdy.U[0],RBdy.U[1],RBdy.U[2],angle[0],angle[1],angle[2],momentArm[0],momentArm[1],momentArm[2]);
    fclose(fp);

    output = outputDir + "/PDL_"+to_string(dataNumber)+"_"+to_string(loop)+".vtu";
    fileIO::export_vtu(x,element,numOfNode,numOfElm,U,volumeChangeRatio,lambda_ave,sigmaEigen_Ave,AEigen_Ave,sigmaEigenVector_Ave,AEigenVector_Ave,innerForce,output);
    //fileIO::export_vtu(x,element,numOfNode,numOfElm,U,volumeChangeRatio,lambda_ave,output);
  }
}


// #################################################################
/**
 * @brief calc theta (Nour-Omid and Rankin, Compt. Methods Appl. Mech. Eng., 1991.)
 * @param [out] ql      rotation angle in local coordinates
 * @param [in] Rbar     rotation matrix [reference to current coorindate (node level)]
 * @param [in] ic element number
 */
void Fem::calc_thetaFromRotationMatrix(double (&ql)[3],const double (&R)[3][3])
{
  double trR,Ra[3][3],tau;

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      Ra[i][j] = R[i][j]-R[j][i];
    }
  }

  ql[0] = Ra[2][1];
  ql[1] = Ra[0][2];
  ql[2] = Ra[1][0];
  tau=5e-1*sqrt(ql[0]*ql[0]+ql[1]*ql[1]+ql[2]*ql[2]);
  if(tau<1e-15){
    for(int i=0;i<3;i++) ql[i]=0e0;
  }else{
    for(int i=0;i<3;i++) ql[i]=5e-1*asin(tau)/tau*ql[i];
  }
}


// #################################################################
/**
 * @brief fem solid analysis routine
 * @param [inout] PARDISO   PARDISO class
 */
void Fem::exportRestartData(const int loop)
{
  FILE *fp;
  string output = "Restart_"+to_string(dataNumber)+"/U_" + to_string(loop) + ".dat";
  if ((fp = fopen(output.c_str(), "w")) == NULL) {
    cout << "file open error" << endl;
    exit(1);
  }
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",U(i,0),U(i,1),U(i,2));
  }
  fclose(fp);
}

// #################################################################
/**
 * @brief fem solid analysis routine
 * @param [inout] PARDISO   PARDISO class
 */
int Fem::NRscheme(PARDISO_solver &PARDISO,RigidBody &RBdy)
{
    double residual,residual0,norm,norm0;
    string output;

    for(int ic=1;ic<=NRiteration;ic++){

      calcStressTensor();  //calc K and Q
      calcTemporalFw(RBdy);

      //elastic body-rigid body interaction
      for(int i=0;i<numOfCP;i++){
        for(int j=0;j<3;j++) innerForce(CP(i),j)+=LAMBDA(i,j);
      }

      set_rhs_statics();

      PARDISO.set_CSR_value(K,element,numOfNode,numOfElm,inb);

      //-----rigid body interaction term-------
      rigidBodyInteraction(RBdy);
      PARDISO.set_CSR_value_rigidBodyInteraction(numOfNode,iCP,Rb,Kqq,numOfCP);
      PARDISO.set_CSR_dirichlet_boundary_condition(numOfNode,ibd);

      for(int i=0;i<numOfNode;i++){
        for(int j=0;j<3;j++){
          PARDISO.b[i+j*numOfNode]=RHS(i,j);
        }
      }
      for(int i=0;i<numOfCP;i++){
        for(int j=0;j<3;j++){
          PARDISO.b[3*numOfNode+i+j*numOfCP]=-Qlambda(i,j);
        }
      }
      for(int j=0;j<3;j++){
        PARDISO.b[3*numOfNode+3*numOfCP+j]=FU[j]-QU[j];
      }
      for(int j=0;j<3;j++){
        PARDISO.b[3*numOfNode+3*numOfCP+3+j]=Fw[j]-Qw[j];
      }
      //---------------------------------------

      PARDISO.main(3*numOfNode+3*numOfCP+6,OMPnumThreads);
      norm = PARDISO.vector_norm(3*numOfNode+3*numOfCP+6,PARDISO.x);
      if(isnan(norm)){
        cout << "norm is nan " << endl;
        exit(1);
      }
      if(ic==1) norm0 = norm;
      corrector_statics(PARDISO.x,relaxation,RBdy);

      for(int i=0;i<numOfNode;i++){
        for(int j=0;j<3;j++) x(i,j) = x0(i,j) + U(i,j);
      }

      // output = outputDir + "/test_NR_" + to_string(ic) + ".vtu";
      // fileIO::export_vtu(x,element,numOfNode,numOfElm,U,volumeChangeRatio,lambda_ave,output);

      // RBdy.updateShape();
      // output = outputDir + "/rigidBody_NR_" + to_string(ic) + ".ply";
      // RBdy.exportPLY(output);

      cout << "NewtonRaphson_iteration = " << ic << endl;
      // cout << " Normalized residual = " <<scientific<< residual/residual0 << " normalized norm = " << scientific<< norm/norm0  <<  endl;
      printf("Tooth displacement=(%e %e %e)\n",RBdy.U[0],RBdy.U[1],RBdy.U[2]);

      if(norm/norm0<NRtolerance) break;
      // if(test!=1 && ic>50) break;
    }
    return 0;
}

void Fem::calcTemporalFw(RigidBody &RBdy)
{
  double momentArm[3];
  for(int i=0;i<3;i++){
    momentArm[i]=0e0;
    for(int j=0;j<3;j++) momentArm[i] += RBdy.R[i][j] * initialMomentArm[j];
  }
  printf("momentArm=%e %e %e\n",momentArm[0],momentArm[1],momentArm[2]);

  double tmp;
  mathTool::crossProduct(momentArm,FU,Fw,tmp);
}


// #################################################################
/**
 * @brief calc stress tensor
 */
void Fem::calcStressTensor()
{
  double elementVolume,volume;

  stress_tensor_initialize();

  volume = 0e0;
  int test;

  // #pragma omp parallel for
  // for(int ic=0;ic<numOfElm;ic++){
  //   test=calcStressTensor_PDL_element_fibreStretch(ic,U,8,2);
  //   if(test==1){
  //     calcStressTensor_PDL_element_spatialForm_hexa_2018(ic,U,8,2,true);
  //   }else{
  //     calcStressTensor_hyperFoam_element_spatialForm_hexa(ic,U,8,2,true);
  //   }
  // }
  // for(int ic=0;ic<numOfElm;ic++) calcStressTensor_ACL_element_spatialForm_hexa_Fbar(ic,U,8,2,true);
  //for(int ic=0;ic<numOfElm;ic++) calcStressTensor_PDL_element_spatialForm_hexa_SRI(ic,U,8,2,true);
  //for(int ic=0;ic<numOfElm;ic++) calcStressTensor_hyperFoam_element_spatialForm_hexa(ic,U,8,2,true);
  // for(int ic=0;ic<numOfElm;ic++) calcStressTensor_PDL_element_spatialForm_hexa_SRI_2018(ic,U,8,2,true);
  // for(int ic=0;ic<numOfElm;ic++) calcStressTensor_PDL_element_spatialForm_hexa_2018(ic,U,8,2,true);
  for(int ic=0;ic<numOfElm;ic++) calcStressTensor_PDL_element_2018(ic,U,8,2);
  totalVolume = volume;

  for(int ic=0;ic<numOfElm;ic++){
    for(int p=0;p<element[ic].node.size();p++){
      for(int i=0;i<3;i++){
        innerForce(element[ic].node[p],i) += Qu(ic,p,i);
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
      // RHS[ic][j] = (double)ibd[ic][j]*(BF[ic][j] - GF[ic][j]);
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
          for(int j=0;j<3;j++) K(ic,p,q,i,j) = 0e0;
        }
      }
      for(int i=0;i<3;i++) Qu(ic,p,i) = 0e0;
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

  //------two point---------
  Gauss gauss(numOfGaussPoint);
  //------------------------

  for(int ic=0;ic<numOfElm;ic++){
    numOfNodeInElm=element[ic].node.size();
    N.allocate(numOfNodeInElm);
    dNdr.allocate(numOfNodeInElm,3);
    X.allocate(numOfNodeInElm,3);

    for(int p=0;p<numOfNodeInElm;p++){
      for(int i=0;i<3;i++) X(p,i) = x0(element[ic].node[p],i);
    }

    for(int p=0;p<numOfNodeInElm;p++){
      for(int q=0;q<numOfNodeInElm;q++) Mass(ic,p,q) = 0e0;
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

          calc_dXdr(dXdr,dNdr,X,numOfNodeInElm);
          detJ = mathTool::calcDeterminant_3x3(dXdr);

          //calc_internal force vector
          for(int p=0;p<numOfNodeInElm;p++){
            for(int q=0;q<numOfNodeInElm;q++){
              Mass(ic,p,q) += rho * N(p) * N(q) * detJ * gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
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
          calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
        }else if(option==0){
          calc_dXdr(dxdr,dNdr,x_ref,numOfNodeInElm);
        }
        calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
        detJ = mathTool::calcDeterminant_3x3(dxdr);

        volume += detJ * gauss.weight[i1] * gauss.weight[i2] * gauss.weight[i3];
      }
    }
  }
  elementVolume(ic)=volume;
  // printf("volume=%e\n",volume);
}
