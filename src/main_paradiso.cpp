//##################################################################################
//
// FEM solid analysis
//
// Copyright (c) 2016-8 Mechanical and Bioengineering Systems Lab.,
//                      Department of Mechanical Science and Bioengineering,
//                      Graduate School of Engineering Science,
//                      Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   main_pardiso.cpp
 * @author T. Otani
 * @detail Fem solid analysis with displacement control.
 */

#include "fem.h"
#include "glog/logging.h"

using namespace std;

int main(int argc,char *argv[])
{
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();

  Fem FEM;
  RigidBody RBdy;
  PARDISO_solver PARDISO;
  string output;

  if(argc!=2){
    cout << "Invalid input" << endl;
    return -1;
  }

  //read tp file
  std::string input_file = argv[1];
  int ierror;
  if ((ierror = FEM.tp.read(input_file)) != TP_NO_ERROR) {
    printf("\tError at reading '%s' file\n", input_file.c_str());
    return 1;
  }

  FEM.initialize();
  RBdy.initialize(FEM.tp);
  FEM.preprocess_rigidBodyInteraction(RBdy);

  omp_set_num_threads(FEM.OMPnumThreads);

  mkdir(FEM.outputDir.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);

  //CSR setting
  PARDISO.initialize(3*FEM.numOfNode,3*FEM.numOfCP);
  PARDISO.CSR_initialize(FEM.inb,FEM.numOfNode,FEM.iCP,FEM.CP,FEM.numOfCP,3);

  output = FEM.outputDir + "/" + FEM.fileName + "_boundary" + ".vtu";
  fileIO::export_vtu_boundary(FEM.x,FEM.element,FEM.numOfNode,FEM.numOfElm,FEM.ibd,FEM.bd,FEM.fiberDirection_elm,output);
  output = FEM.outputDir + "/start.ply";
  RBdy.exportPLY(output);

  cout << "---------main loop start----------" << endl;

  FEM.femSolidAnalysis(PARDISO,RBdy);

  return 0;
}
