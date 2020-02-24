/**
 * @file rigidElasticInteraction.cpp
 * @brief RigidElasticInteraction class
 * @author T. Otani
 */

#include "rigidElasticInteraction.h"

using namespace std;

// #################################################################
/**
 * @brief rigid body-elastic body interaction problem
 */
void RigidElasticInteraction::mainLoop()
{
  int output_iter=1;
  string output;

  //------------------------------------------------------------
  FILE *fp;
  string outputFile_tooth = outputDir + "/toothDisplacement.dat";
  if ((fp = fopen(outputFile_tooth.c_str(), "w")) == NULL) {
    exit(1);
  }
  fprintf(fp,"loop F Ux Uy Uz Qx Qy Qz ARMx ARMy ARMz\n");
  fclose(fp);
  //------------------------------------------------------------

  for(int j=0;j<3;j++) initialMomentArm[j] = FUpoint[j]-RBdy.xg[j];
  // for(int ic=0;ic<numOfElm;ic++) calcVolume_hexa(ic,volume0,8,2,false);

  //linear elastic material only
  // stress_tensor_initialize();
  // for(int ic=0;ic<numOfElm;ic++) calcStressTensor_LinearElastic_element_spatialForm(ic,true);

  for(int loop=1;loop<=maxIteration;loop++){

    for(int j=0;j<3;j++) FU[j] = (double)loop/(double)maxIteration * FU_input[j];

    if(NRscheme()==1){
      printf("NR scheme is wrong. Loop=%d\n",loop);
      exit(1);
    }

    for(int i=0;i<numOfNode;i++){
      for(int j=0;j<3;j++) x(i,j) = x0(i,j) + U(i,j);
    }
    // for(int i=0;i<numOfElm;i++){
    //   calcVolume_hexa(i,volume,8,2,true);
    //   volumeChangeRatio(i)=volume(i)/volume0(i);
    // }


    // for(int ic=0;ic<numOfElm;ic++) postProcess_LinearElastic_element_spatialForm(ic,true);
    for(int ic=0;ic<numOfElm;ic++) postProcess_PDL_element_2018(ic,U,8,2);
    // exportRestartData(loop);

    RBdy.updateShape();
    output = outputDir + "/tooth_"+to_string(dataNumber)+"_"+to_string(loop)+".ply";
    RBdy.exportPLY(output);

    double angle[3];
    mathTool::calc_thetaFromRotationMatrix(angle,RBdy.R);

    double momentArm[3];
    for(int i=0;i<3;i++){
      momentArm[i]=0e0;
      for(int j=0;j<3;j++) momentArm[i] += RBdy.R[i][j] * initialMomentArm[j];
      momentArm[i]-=initialMomentArm[i];
    }

    if ((fp = fopen(outputFile_tooth.c_str(), "a")) == NULL) {
      exit(1);
    }
    // double ForceMagnitude= (double)loop/(double)maxIteration * sqrt(FU_input[0]*FU_input[0]+FU_input[1]*FU_input[1]+FU_input[2]*FU_input[2]);
    // fprintf(fp,"%d %e %e %e %e %e %e %e %e %e %e\n",loop,ForceMagnitude,RBdy.U[0],RBdy.U[1],RBdy.U[2],angle[0],angle[1],angle[2],momentArm[0],momentArm[1],momentArm[2]);
    // fclose(fp);

    output = outputDir + "/PDL_"+to_string(dataNumber)+"_"+to_string(loop)+".vtu";
    // fileIO::export_vtu_Mises(x,element,numOfNode,numOfElm,U,Mises,output);
    fileIO::export_vtu(x,element,numOfNode,numOfElm,U,volumeChangeRatio,lambda_ave,sigmaEigen_Ave,sigmaEigenVector_Ave,output);
    exit(1);
  }
}

// #################################################################
/**
 * @brief fem solid analysis routine
 */
int RigidElasticInteraction::NRscheme()
{
  double residual,residual0,norm,norm0;
  string output;

  for(int ic=1;ic<=NRiteration;ic++){

    calcStressTensor();  //calc K and Q

    #pragma omp parallel for
    for(int i=0;i<numOfNode;i++){
      for(int j=0;j<3;j++) innerForce(i,j) = 0e0;
    }

    //linear elastic material only
    for (int ia = 0; ia < numOfElm; ia++){
      for (int p = 0; p < element[ia].node.size(); p++){
       for (int q = 0; q < element[ia].node.size(); q++){
        for (int i = 0; i < 3; i++){
          for (int j = 0; j < 3; j++) innerForce(element[ia].node[p],i) += K(ia,p,q,i,j) * U(element[ia].node[q],j);
        }
      }
    }
  }

    calcTemporalFw();

    //elastic body-rigid body interaction
    for(int i=0;i<numOfCP;i++){
      for(int j=0;j<3;j++) innerForce(CP(i),j)+=LAMBDA(i,j);
    }

    set_rhs_statics();
    PARDISO.set_CSR_value3D(K,element,numOfNode,numOfElm,inb);

    //-----rigid body interaction term-------
    calcRigidBodyInteractionTerm(RBdy);
    PARDISO.set_CSR_value_rigidBodyInteraction(numOfNode,iCP,Rb,Kqq,numOfCP);
    PARDISO.set_CSR_dirichlet_boundary_condition(numOfNode,ibd);

    //rhs term
    for(int i=0;i<numOfNode;i++){
      for(int j=0;j<3;j++) PARDISO.b[i+j*numOfNode]=RHS(i,j);
    }
    for(int i=0;i<numOfCP;i++){
      for(int j=0;j<3;j++) PARDISO.b[3*numOfNode+i+j*numOfCP]=-Qlambda(i,j);
    }
    for(int j=0;j<3;j++){
      PARDISO.b[3*numOfNode+3*numOfCP+j]  = FU[j]-QU[j];
      PARDISO.b[3*numOfNode+3*numOfCP+3+j]= Fw[j]-Qw[j];
    }

    PARDISO.main(3*numOfNode+3*numOfCP+6,OMPnumThreads);
    norm = PARDISO.vector_norm(3*numOfNode+3*numOfCP+6,PARDISO.x);

    if(isnan(norm)){
      cout << "norm is nan " << endl;
      return 1;
    }
    if(ic==1) norm0 = norm;
    corrector_statics(PARDISO.x,relaxation);

    for(int i=0;i<numOfNode;i++){
      for(int j=0;j<3;j++) x(i,j) = x0(i,j) + U(i,j);
    }

    //for debug
    // output = outputDir + "/test_NR_" + to_string(ic) + ".vtu";
    // fileIO::export_vtu(x,element,numOfNode,numOfElm,U,volumeChangeRatio,lambda_ave,output);

    //for debug
    // RBdy.updateShape();
    // output = outputDir + "/rigidBody_NR_" + to_string(ic) + ".ply";
    // RBdy.exportPLY(output);
    printf("NR iter.=%d norm/norm0=%e\n",ic,norm/norm0);
    //printf("NR iter=%d Tooth displacement=(%e %e %e)\n",ic,RBdy.U[0],RBdy.U[1],RBdy.U[2]);

    if(norm/norm0<NRtolerance) break;
    // if(test!=1 && ic>50) break;
  }
  return 0;
}

// #################################################################
/**
 * @brief temporal fw.
 */
void RigidElasticInteraction::calcTemporalFw()
{
  double momentArm[3];
  for(int i=0;i<3;i++){
    momentArm[i]=0e0;
    for(int j=0;j<3;j++) momentArm[i] += RBdy.R[i][j] * initialMomentArm[j];
  }
  //printf("momentArm=%e %e %e\n",momentArm[0],momentArm[1],momentArm[2]);

  double tmp;
  mathTool::crossProduct(momentArm,FU,Fw,tmp);
}

// #################################################################
/**
 * @brief corrector scheme.
 * @param [in] u           displacement vector
 * @param [in] relaxation  relaxation parameters
 */
void RigidElasticInteraction::corrector_statics(const double *u,const double relaxation)
{
  double w[3];

  #pragma omp parallel for
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) U(i,j) += u[i+j*numOfNode]*relaxation;
  }

  #pragma omp parallel for
  for(int i=0;i<numOfCP;i++){
    for(int j=0;j<3;j++) LAMBDA(i,j) += u[3*numOfNode+i+j*numOfCP]*relaxation;
  }

  for(int j=0;j<3;j++) RBdy.U[j] += u[3*numOfNode+3*numOfCP+j]*relaxation;
  for(int j=0;j<3;j++) w[j]       = u[3*numOfNode+3*numOfCP+3+j]*relaxation;

  RBdy.updateRotationMatrix_spatialForm(w);
}

// #################################################################
/**
 * @brief calc b0
 * @param [in] RBdy          rigid body class
 */
void RigidElasticInteraction::initialize_rigidBodyInteraction()
{
  initialize();
  RBdy.initialize(tp);
  inputRigidBodyInterface();

  LAMBDA.allocate(numOfCP,3);
  Rb.allocate(numOfCP,3,3);
  b0.allocate(numOfCP,3);
  b.allocate(numOfCP,3);
  Qlambda.allocate(numOfCP,3);

  for(int i=0;i<numOfCP;i++){
    for(int j=0;j<3;j++) LAMBDA(i,j)=0e0;
  }

  for(int ic=0;ic<numOfCP;ic++){
    for(int j=0;j<3;j++) b0(ic,j)=x(CP(ic),j)-RBdy.xg[j];
  }

  omp_set_num_threads(OMPnumThreads);

  mkdir(outputDir.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);

  for(int i=0;i<numOfNode;i++){
      if(ibd(i,0)==0 && iCP(i)!=-1) cout << i << endl;
  }

  //CSR setting
  PARDISO.initialize(3*numOfNode,3*numOfCP);
  PARDISO.CSR_initialize(inb,numOfNode,iCP,CP,numOfCP,3);

  string output = outputDir + "/" + fileName + "_boundary" + ".vtu";
  fileIO::export_vtu_boundary(x,element,numOfNode,numOfElm,ibd,bd,fiberDirection_elm,output);
  output = outputDir + "/start.ply";
  RBdy.exportPLY(output);
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