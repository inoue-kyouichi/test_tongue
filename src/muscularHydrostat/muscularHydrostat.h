#ifndef _MUSCULAR_HYDROSTAT_H_
#define _MUSCULAR_HYDROSTAT_H_

//##################################################################################
//
// Eye mechanics simulation
//
// Copyright (c) 2020 Biomechanics Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   eye_Mechanics.h
 * @brief  PDL Header
 * @author T. Otani
 */

#include "fem.h"
#include "pardiso_solver.h"

namespace muscularHydrostat {

class MaterialModel{
  public:
  double contractionCoefficient;
  double initialStretch;
};

typedef enum{
  M0 = 0,
  M1 = 1,
  numOfMaterials
}MaterialType;

class Muscle : public Fem {
 public:
  TextParser tp;
  std::string outputDir,fileName;

  // int numOfBoundaryElm;
  // VECTOR1D<ElementType> boundaryElement;
  VECTOR1D<MaterialModel> Material;

 private:
  int dataNumber;
  int Restart;
  int OMPnumThreads;
  int maxIteration,NRiteration;
  double NRtolerance;
  double relaxation;
  PARDISO_solver PARDISO;
  ARRAY2D<double> fiberDirection_elm;

 public: 
  void preprocess();
  void femSolidAnalysis();
  void calcStressTensor();

 private:
  bool NRscheme();
  void set_rhs_statics();
  void calcStressTensor_SantVenant_element_spatialForm(const int ic,ARRAY2D<double> &U_tmp,const bool option);
  double SantVenant_inGaussIntegral(ARRAY2D<double> &dNdr,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,
        ARRAY2D<double> &dNdx,const int numOfNodeInElm,const double weight,const int ic,const bool option);

  void calcStressTensor_muscle_element_spatialForm(const int ic,ARRAY2D<double> &U_tmp,const bool option);
  void muscle_inGaussIntegral(const int &ic,const int &numOfNodeInElm,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,ARRAY2D<double> &dNdr,const double weight,double (&stress)[3][3],const bool mainLoop);
  void calc_F_initial(double (&F_initial)[3][3],const double (&a0)[3],const double lambda);

  void inputMaterialParameters(TextParser &tp);
  void inputSolverInfo(TextParser &tp);
  void inputOutputInfo(TextParser &tp);
  void inputFiberInfo(TextParser &tp);
  void inputMaterialInfo(TextParser &tp);
  void inputDirichletInfo(TextParser &tp);

  void export_vtu(const std::string &file);

};

}

#endif