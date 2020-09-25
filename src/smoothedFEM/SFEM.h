#ifndef _EYE_MECH_H_
#define _EYE_MECH_H_

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


class FaceType{
 public:
 FaceType(){
   meshType = VTK_TRIANGLE;
   element[0] = -1; element[1] = -1;
 };
 ~FaceType(){};
  VTKCellType meshType;
  VECTOR1D<int> node;
  int element[2];
};

namespace SmoothedFEM {

class SFEM : public Fem {
 public:
  TextParser tp;
  std::string outputDir,fileName;

  int numOfBoundaryElm;
  VECTOR1D<ElementType> boundaryElement;

 private:
  int dataNumber;
  int Restart;
  int OMPnumThreads;
  int maxIteration,NRiteration;
  double NRtolerance;
  double relaxation;
  PARDISO_solver PARDISO;

 public: 
  double boundaryPressure;
  ARRAY2D<double> boundaryForce;
  std::vector<ARRAY4D<double>> pressureStiffnessMatrix;

  std::vector<FaceType> face;
  ARRAY1D<double> faceVolume;

  void preprocess();
  void femSolidAnalysis();
  void femSolidAnalysis_linear();

  void calcStressTensor();

  void setFace();
  void calcFaceVolume();

 private:
  bool NRscheme();
  void set_rhs_statics();

  void calcStressTensor_linearElasticMaterial_element_SFEM(const int ic,const bool option);
  double linearElasticMaterial_SFEM(ARRAY2D<double> &dNdr,ARRAY2D<double> &x_ref,
         const int numOfNodeInElm,const double weight,const int ic,const bool option);

  void calcStressTensor_linearElasticMaterial_element(const int ic,const bool option);
  double linearElasticMaterial_inGaussIntegral(ARRAY2D<double> &dNdr,ARRAY2D<double> &x_ref,
        const int numOfNodeInElm,const double weight,const int ic,const bool option);

  void calcStressTensor_SantVenant_element_spatialForm(const int ic,ARRAY2D<double> &U_tmp,const bool option);
  double SantVenant_inGaussIntegral(ARRAY2D<double> &dNdr,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,
        ARRAY2D<double> &dNdx,const int numOfNodeInElm,const double weight,const int ic,const bool option);

  void calcBoundaryForce();
  void calcBoundaryPressure_spatialForm(const int ic,VECTOR1D<ARRAY2D<double>> &Qb,const double boundaryPressure);
  void boundaryPressure_inGaussIntegral(VECTOR1D<ARRAY2D<double>> &Qb,ARRAY1D<double> &N,ARRAY2D<double> &dNdr,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,const int numOfNodeInElm,const double boundaryPressure,const double weight,const int ic);

  void calcStressTensor_NeoHookean_element_spatialForm(const int ic,const double mu,const double Poisson,ARRAY2D<double> &U_tmp,const bool option);
  void NeoHookean_inGaussIntegral(ARRAY2D<double> &dNdr,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,const int numOfNodeInElm,const double mu,const double lambda,const double weight,const int ic,const bool option);

  void inputSolverInfo(TextParser &tp);
  void inputOutputInfo(TextParser &tp);

  void inputMaterialInfo(TextParser &tp);
  void inputSurfaceInfo(TextParser &tp);
  void inputDirichletInfo(TextParser &tp);

  void constitutiveMatrix(double (&D)[6][6],double (&C4)[3][3][3][3]);

  void export_vtu(const std::string &file);
  void export_vtu_face(const std::string &file);
};

}

#endif