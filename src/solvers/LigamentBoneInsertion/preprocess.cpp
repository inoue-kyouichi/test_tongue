/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "LigamentBoneInsertion.h"
#include <ostream>
#include <fstream>

using namespace std;


void InsertionSite::preprocess()
{
  initialize(tp);
  inputSolverInfo(tp);
  inputOutputInfo(tp);

  omp_set_num_threads(OMPnumThreads);

  mkdir(outputDir.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);

  //CSR setting
  PARDISO.initialize(3*numOfNode);
  PARDISO.CSR_initialize(inb,numOfNode,3);

  string output = outputDir + "/" + fileName + "_boundary" + ".vtu";
  fileIO::export_vtu_boundary(x,element,numOfNode,numOfElm,ibd,bd,output);
}


// #################################################################
/**
 * @brief solver information from TP file
 */
void InsertionSite::inputSolverInfo(TextParser &tp)
{
  string str,base_label,label;
  int tmp;

  base_label = "/Solver";

  label = base_label + "/dataNumber";
  if ( !tp.getInspectedValue(label, dataNumber)){
    cout << "maxiteration is not set" << endl;
    exit(0);
  }

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
 * @brief output information from TP file
 */
void InsertionSite::inputOutputInfo(TextParser &tp)
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