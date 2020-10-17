/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "muscularHydrostat.h"
#include <ostream>
#include <fstream>

using namespace std;

void muscularHydrostat::Muscle::preprocess()
{
  inputDomainInfo(tp);
  allocate();

  fiberDirection_elm.allocate(numOfElm,3);
  inputFiberInfo(tp);

  inputDirichletInfo(tp);

  inputSolverInfo(tp);
  inputOutputInfo(tp);

  omp_set_num_threads(OMPnumThreads);

  mkdir(outputDir.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);

  //CSR setting
  PARDISO.initialize(3*numOfNode);
  PARDISO.CSR_initialize(inb,numOfNode,3);
}

// #################################################################
/**
 * @brief solver information from TP file
 */
void muscularHydrostat::Muscle::inputSolverInfo(TextParser &tp)
{
  string str,base_label,label;
  int tmp;

  base_label = "/Solver";

  label = base_label + "/maxIteration";
  if ( !tp.getInspectedValue(label, maxIteration)){
    cout << "maxiteration is not set" << endl;
    exit(0);
  }

  label = base_label + "/NR_iteration";
  if ( !tp.getInspectedValue(label, NRiteration)){
    cout << "NLiteration is not set" << endl;
    exit(0);
  }

  label = base_label + "/NR_tolerance";
  if ( !tp.getInspectedValue(label, tmp)){
    cout << "NRtolerance is not set" << endl;
    exit(0);
  }
  NRtolerance = pow(1e-1,tmp);

  label = base_label + "/Restart";
  if ( !tp.getInspectedValue(label, Restart)){
    cout << label <<" is not set" << endl;
    exit(0);
  }

  label = base_label + "/OMPnumThreads";
  if ( !tp.getInspectedValue(label, OMPnumThreads)){
    cout << "OMPnumThreads is not set" << endl;
    exit(0);
  }

  label = base_label + "/relaxation";
  if ( !tp.getInspectedValue(label,relaxation)){
    cout << "NLiteration is not set" << endl;
    exit(0);
  }
}

// #################################################################
/**
 * @brief tentative
 */
void muscularHydrostat::Muscle::inputFiberInfo(TextParser &tp)
{
  string str,base_label,label,inputDir;
  base_label = "/Domain";

  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }
  string D_file,Dvalue_file;
  label = base_label + "/fiberFile";
  if ( !tp.getInspectedValue(label, D_file)){
    cout << label << " is not found" << endl;
    exit(1);
  }
  D_file=inputDir+"/"+D_file;

  ifstream file(D_file);
  if(!file){
    cout << "Error:Input "<< D_file << " not found" << endl;
    exit(1);
  }

  string tmp,dtmp;
  for(int i=0;i<numOfElm;i++){
    getline(file,str);
    istringstream stream(str);
    for(int j=0;j<3;j++){
      getline(stream,tmp,' ');
      fiberDirection_elm(i,j) = stod(tmp);
    }
  }
}

// #################################################################
/**
 * @brief tentative
 */
void muscularHydrostat::Muscle::inputMaterialInfo(TextParser &tp)
{
  string str,base_label,label,inputDir;
  base_label = "/Domain";

  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }
  string D_file,Dvalue_file;
  label = base_label + "/materialTypeFile";
  if ( !tp.getInspectedValue(label, D_file)){
    cout << label << " is not found" << endl;
    exit(1);
  }
  D_file=inputDir+"/"+D_file;

  ifstream file(D_file);
  if(!file){
    cout << "Error:Input "<< D_file << " not found" << endl;
    exit(1);
  }

  string tmp,dtmp;
  for(int i=0;i<numOfElm;i++){
    getline(file,str);
    element[i].materialType = stoi(str);
  }

}

// #################################################################
/**
 * @brief Dirichlet information from TP file
 */
void muscularHydrostat::Muscle::inputDirichletInfo(TextParser &tp)
{
  string str,base_label,label,inputDir;
  base_label = "/Domain";

  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }
  string D_file,Dvalue_file;
  label = base_label + "/dirichletFile";
  if ( !tp.getInspectedValue(label, D_file)){
    cout << label << " is not found" << endl;
  }
  D_file=inputDir+"/"+D_file;

  int numOfData = fileIO::CountNumbersOfTextLines(D_file);

  int node;
  string tmp,dtmp;

  ifstream file(D_file);
  if(!file){
    cout << "Error:Input "<< D_file << " not found" << endl;
    exit(1);
  }

  ibd.allocate(numOfNode,3);
  bd.allocate(numOfNode,3);

  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++){
      bd(i,j) = 0e0;
      ibd(i,j) = 1;
    }
  }

  for(int i=0;i<numOfData;i++){
    getline(file,str);
    istringstream stream(str);
    getline(stream,tmp,' ');
    int number = stoi(tmp);
    for(int j=0;j<3;j++){
      getline(stream,tmp,' ');
      ibd(number,j) = stoi(tmp);
    }
  }

}


// #################################################################
/**
 * @brief output information from TP file
 */
void muscularHydrostat::Muscle::inputOutputInfo(TextParser &tp)
{
  string str,base_label,label;

  base_label = "/Output";
  label = base_label + "/outputDir";
  if ( !tp.getInspectedValue(label, outputDir)){
    cout << "outputDir is not set" << endl;
    exit(0);
  }

  label = base_label + "/fileName";
  if ( !tp.getInspectedValue(label, fileName)){
    cout << "fileName is not set" << endl;
    exit(0);
  }
}