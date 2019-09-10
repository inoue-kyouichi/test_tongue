/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "rigidElasticInteraction.h"
#include <ostream>
#include <fstream>

using namespace std;

// #################################################################
/**
 * @brief fem solid analysis routine
 * @param [inout] PARDISO   PARDISO class
 */
void RigidElasticInteraction::femSolidAnalysis(PARDISO_solver &PARDISO,RigidBody &RBdy)
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
void RigidElasticInteraction::calc_thetaFromRotationMatrix(double (&ql)[3],const double (&R)[3][3])
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
int RigidElasticInteraction::NRscheme(PARDISO_solver &PARDISO,RigidBody &RBdy)
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

void RigidElasticInteraction::calcTemporalFw(RigidBody &RBdy)
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
 * @brief domain information from tp file
 */
void RigidElasticInteraction::inputRigidBodyInterface()
{
  string str,base_label,label,inputDir;

  base_label = "/Domain";

  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }
  string file;
  label = base_label + "/interface";
  if ( !tp.getInspectedValue(label, file)){
    cout << label << " is not found" << endl;
    exit(0);
  }

  file=inputDir + "/" + file;

  numOfCP = fileIO::CountNumbersOfTextLines(file);
  CP.allocate(numOfCP);
  iCP.allocate(numOfNode);

  FILE *fp;
  if ((fp = fopen(file.c_str(), "r")) == NULL) {
    cout << "file open error" << endl;
    exit(1);
  }
  for(int i=0;i<numOfCP;i++) fscanf(fp,"%d\n",&CP(i));
  fclose(fp);

  for(int i=0;i<numOfNode;i++) iCP(i)=-1;
  for(int i=0;i<numOfCP;i++) iCP(CP(i))=i;

  base_label = "/RigidBody";
  label = base_label + "/Force";
  if ( !tp.getInspectedVector(label,FU_input,3)){
    cout << label << " is not found" << endl;
    exit(0);
  }

  label = base_label + "/ForcePoint";
  if ( !tp.getInspectedVector(label,FUpoint,3)){
    cout << label << " is not found" << endl;
    exit(0);
  }
  // label = base_label + "/Moment";
  // if ( !tp.getInspectedVector(label,Fw,3)){
  //   cout << label << " is not found" << endl;
  //   exit(0);
  // }
}
