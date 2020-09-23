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

typedef enum{
  ligament = 0,
  interface = 1,
  bone = 2,
  numOfMaterials
}MaterialType;

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
  ARRAY3D<double> fiberDirection;
  ARRAY2D<double> lambda_ave;

  ARRAY2D<double> fiberDirection_elm;

  void preprocess();
  void femSolidAnalysis();
  void calcStressTensor();

  // void postProcess_PDL_element_spatialForm_hexa_SRI(const int &ic,DOUBLEARRAY2D &U_tmp,
  // const int &numOfNodeInElm,const int &numOfGaussPoint);

 private:
  bool NRscheme();
  void set_rhs_statics();
  void calcStressTensor_SantVenant_element_spatialForm(const int ic,ARRAY2D<double> &U_tmp,const bool option);
  double SantVenant_inGaussIntegral(ARRAY2D<double> &dNdr,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,
        ARRAY2D<double> &dNdx,const int numOfNodeInElm,const double weight,const int ic,const bool option);

  void inputSolverInfo(TextParser &tp);
  void inputOutputInfo(TextParser &tp);

};

#endif