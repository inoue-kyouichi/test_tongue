/**
 * @file minimumComplianceProblem.cpp
 * @brief minimumComplianceProblem class
 * @author T. Otani
 */

#include "minimumComplianceProblem.h"

using namespace std;

// #################################################################
/**
 * @brief initialize and preprocess
 */
void MinimumComplianceProblem::femSolidAnalysisClass::mainLoop()
{
  int output_iter = 1;
  int increment = 1;
  string output;

  #pragma omp parallel for
  for (int i = 0; i < elasticBody.numOfNode; i++)
  {
      for (int j = 0; j < 3; j++)
          elasticBody.U(i, j) = 0e0;
  }

  if (NRscheme() == 1)
  {
      printf("calculation error. Exit...");
      exit(1);
  }

  #pragma omp parallel for
  for (int i = 0; i < elasticBody.numOfNode; i++)
  {
      for (int j = 0; j < 3; j++)
          elasticBody.x(i, j) = elasticBody.x0(i, j) + elasticBody.U(i, j);
  }
}

// #################################################################
/**
 * @brief initialize and preprocess
 */
bool MinimumComplianceProblem::femSolidAnalysisClass::NRscheme()
{
  double residual, residual0, norm, norm0;

  for (int ic = 1; ic <= NRiteration; ic++)
  {
    elasticBody.calcStressTensor(); //calc K and Q

    //left hand side
    PARDISO.set_CSR_value3D(elasticBody.Ku, elasticBody.element, elasticBody.numOfNode, elasticBody.numOfElm, elasticBody.inb);
    PARDISO.set_CSR_dirichlet_boundary_condition3D(elasticBody.numOfNode,elasticBody.ibd);

    //right hand side
    calcExternalSurfaceForce();
    elasticBody.set_rhs_statics();
    for(int i = 0; i < elasticBody.numOfNode; i++){
      for(int j = 0; j < 3; j++){
        PARDISO.b[i + j * elasticBody.numOfNode] = elasticBody.RHS(i, j);
      }
    }

    PARDISO.main(3 * elasticBody.numOfNode, OMPnumThreads);

    norm = PARDISO.vector_norm(elasticBody.numOfNode*3, PARDISO.x);
    if (ic == 1) norm0 = norm;

    corrector_statics(elasticBody.U,PARDISO.x,elasticBody.numOfNode,relaxation);
    //residual=line_search(PARDISO.x);

    // if(isnan(residual)){
    //   cout << "residual is nan. Exit..." << endl;
    //   return 1;
    // }

    #pragma omp parallel for
    for (int i = 0; i < elasticBody.numOfNode; i++)
    {
      for (int j = 0; j < 3; j++) elasticBody.x(i, j) = elasticBody.x0(i, j) + elasticBody.U(i, j);
    }
    // string output = outputDir + "/test_NR_" + to_string(ic) + ".vtu";
    // fileIO::export_vtu(elasticBody.x,elasticBody.element,elasticBody.numOfNode,elasticBody.numOfElm,elasticBody.U,output);

    cout << "NewtonRaphson_iteration = " << ic << endl;
    printf("norm=%e and normalized norm=%e\n", norm, norm / norm0);
    if (norm / norm0 < NRtolerance) break;
  }
  return 0;
}

void MinimumComplianceProblem::femSolidAnalysisClass::calcExternalSurfaceForce()
{
  for(int ic=0;ic<elasticBody.numOfNode;ic++){
    for(int j=0;j<3;j++) elasticBody.externalForce(ic,j) = 0e0;
    for(int j=0;j<3;j++) elasticBody.externalSurfaceForce(ic,j) = 0e0;
  }  

  ARRAY2D<double> Traction(elasticBody.belement.size(),3);
  for(int ic=0;ic<elasticBody.belement.size();ic++){
    Traction(ic,0) = 0e0;
    Traction(ic,1) = 0e0;
    Traction(ic,2) = 1e5;
  }

  elasticBody.calc_externalSurfaceForce_prescribedTraction(elasticBody.belement,Traction);

  for(int i=0;i<elasticBody.numOfNode;i++){
    for(int j=0;j<3;j++) elasticBody.externalForce(i,j) += elasticBody.externalSurfaceForce(i,j);
  }
}

// #################################################################
/**
 * @brief corrector scheme.
 * @param [in] u           displacement vector
 * @param [in] relaxation  relaxation parameters
 */
void MinimumComplianceProblem::femSolidAnalysisClass::corrector_statics(ARRAY2D<double> &U, const double *u, const int numOfNode, const double relaxation)
{
  #pragma omp parallel for
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) U(i,j) += u[i+j*numOfNode]*relaxation;
  }
}

// #################################################################
/**
 * @brief initialize and preprocess
 */
void MinimumComplianceProblem::femSolidAnalysisClass::preprocess()
{
  elasticBody.initialize(tp);
  elasticBody.inputMaterialParameters(tp);
  elasticBody.inputSurfaceBoundary(tp);

  inputSolverInfo(tp);
  inputOutputInfo(tp);

  omp_set_num_threads(OMPnumThreads);
  mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

  //CSR setting
  PARDISO.initialize(3*elasticBody.numOfNode);
  PARDISO.CSR_initialize(elasticBody.inb, elasticBody.numOfNode, 3);

  string output = outputDir + "/" + fileName + "_boundary" + ".vtu";
  fileIO::export_vtu_boundary(elasticBody.x, elasticBody.element, elasticBody.numOfNode, elasticBody.numOfElm, elasticBody.ibd, elasticBody.bd, output);
}

// #################################################################
/**
 * @brief solver information from TP file
 */
void MinimumComplianceProblem::femSolidAnalysisClass::inputSolverInfo(TextParser &tp)
{
  string str,base_label,label;
  int tmp;

  base_label = "/Solver";

  // label = base_label + "/maxIteration";
  // if ( !tp.getInspectedValue(label, maxIteration)){
  //   cout << "maxiteration is not set" << endl;
  //   exit(0);
  // }

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
void MinimumComplianceProblem::femSolidAnalysisClass::inputOutputInfo(TextParser &tp)
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
