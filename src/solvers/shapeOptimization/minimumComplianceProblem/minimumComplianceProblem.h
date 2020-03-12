#ifndef _MINIMUM_COMPLIANCE_H_
#define _MINIMUM_COMPLIANCE_H_

//##################################################################################
//
// FEM Base
//
// Copyright (c) 2019 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   fem.h
 * @brief  FEMBase Definition Header
 * @author T. Otani
 */
#include "StVenant_KirchhoffMaterial.h"
#include "pardiso_solver.h"

namespace MinimumComplianceProblem{
///main problem class
class femSolidAnalysisClass{

 public:
   StVenantKirchhoffMaterial elasticBody;
  void preprocess();
  void mainLoop();

  TextParser tp;
  std::string outputDir,fileName;

 private:

  PARDISO_solver PARDISO;

  int dataNumber;
  int maxIteration,NRiteration;
  double NRtolerance;
  int Restart;
  int OMPnumThreads;
  double relaxation;

  bool NRscheme();
  void calcExternalSurfaceForce();
  void corrector_statics(DOUBLEARRAY2D &U, const double *u, const int numOfNode, const double relaxation);
  void inputSolverInfo(TextParser &tp);
  void inputOutputInfo(TextParser &tp);
};

//shape optimization class which has above mainProblem class
class minimumComplianceProblem{
 public:
  femSolidAnalysisClass mainProblem;
  void preprocess();
  void mainLoop();

  // TextParser tp;
  // std::string outputDir,fileName;

 private:

  // PARDISO_solver PARDISO;

  // int dataNumber;
  // int maxIteration,NRiteration;
  // double NRtolerance;
  // int Restart;
  // double relaxation;

  DOUBLEARRAY2D adjointV;

//   void femSolidAnalysis();
//   void NRscheme();
// //   void inputRigidBodyInterface();
// //   int NRscheme();
// //   void set_rhs_statics();
//   void inputSolverInfo(TextParser &tp);
//   void inputOutputInfo(TextParser &tp);
};
}

#endif