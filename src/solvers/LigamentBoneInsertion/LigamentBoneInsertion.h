#ifndef _LIGAMENT_BONE_INSERTION_H_
#define _LIGAMENT_BONE_INSERTION_H_


//##################################################################################
//
// Ligament Bone insertion
//
// Copyright (c) 2020 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   LigamentBoneInsertion.h
 * @brief  PDL Header
 * @author T. Otani
 */

#include "fem.h"
#include "pardiso_solver.h"

class InsertionSite : public Fem {
 public:
  TextParser tp;
  std::string outputDir,fileName;

 private:
  int dataNumber;
  int Restart;
  int OMPnumThreads;
  int maxIteration,NRiteration;
  double NRtolerance;
  double relaxation;
  PARDISO_solver PARDISO;

 public:
  double bulkModulusRatio;
  DOUBLEARRAY3D fiberDirection;
  DOUBLEARRAY2D lambda_ave;

  DOUBLEARRAY2D fiberDirection_elm;

  void preprocess();
  void femSolidAnalysis();
  void calcStressTensor();

  // void postProcess_PDL_element_spatialForm_hexa_SRI(const int &ic,DOUBLEARRAY2D &U_tmp,
  // const int &numOfNodeInElm,const int &numOfGaussPoint);

 private:
  bool NRscheme();
  void set_rhs_statics();
  void calcStressTensor_SantVenant_element_spatialForm(const int ic,DOUBLEARRAY2D &U_tmp,const bool option);
  double SantVenant_inGaussIntegral(DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &x_current,DOUBLEARRAY2D &x_ref,
        DOUBLEARRAY2D &dNdx,const int numOfNodeInElm,const double weight,const int ic,const bool option);

  void inputSolverInfo(TextParser &tp);
  void inputOutputInfo(TextParser &tp);

};

#endif