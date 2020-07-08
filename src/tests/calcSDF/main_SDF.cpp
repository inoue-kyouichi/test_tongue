//##################################################################################
//
// test: calc surface curvature
//
// Copyright (c) 2020  Mechanical and Bioengineering Systems Lab.,
//                     Department of Mechanical Science and Bioengineering,
//                     Graduate School of Engineering Science,
//                     Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   main_curvature.cpp
 * @author T. Otani
 * @detail Fem solid analysis with displacement control.
 */

#include "SignedDistanceFunction.h"
#include "glog/logging.h"


void inputMaterialInfo(TextParser &tp,Fem &FEM);

using namespace std;

int main(int argc,char *argv[])
{
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();

  int numOfNode,numOfElm;
  SignedDistanceFunction SDF;
  TextParser tp;

  //read triangle surface
  if(argc!=2){
    cout << "please set input file name. Exit..." << endl;
    exit(1);
  }
  std::string input_file = argv[1];
  int ierror;
  if ((ierror = tp.read(input_file)) != TP_NO_ERROR) {
    printf("\tError at reading '%s' file\n", input_file.c_str());
    return 1;
  }

  omp_set_num_threads(2);

  SDF.FEMorg.inputDomainInfo(tp);
  SDF.set1stTetra();

  SDF.SDF.allocate(SDF.FEM.numOfNode);
  inputMaterialInfo(tp,SDF.FEM);

  SDF.calcElementsPerNodes();
  SDF.calcKnownNodes();
  //prprocess completed.

  SDF.calcSDF();


  std::string output_file = "test.vtu";
  SDF.exportVTU(output_file);
  return 0;
}

// #################################################################
/**
 * @brief domain information from tp file
 */
void inputMaterialInfo(TextParser &tp,Fem &FEM)
{
  string str,base_label,label,inputDir;

  base_label = "/Domain";

  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }
  string file4;

  label = base_label + "/materialTypeFile";
  if ( !tp.getInspectedValue(label,file4)){
    cout << label << " is not found" << endl;
    exit(0);
  }
  file4=inputDir+"/"+file4;

  fileIO::read_geometry_materialType(FEM.element,FEM.numOfElm,file4);
}