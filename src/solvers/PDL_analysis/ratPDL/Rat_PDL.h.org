#ifndef _RAT_PDL_H_
#define _RAT_PDL_H_


//##################################################################################
//
// PDL
//
// Copyright (c) 2020 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   PDL.h
 * @brief  PDL Header
 * @author T. Otani
 */

#include "fem.h"

typedef enum {
  PULP       = 0,
  PDL        = 1,
  NUMBER_OF_MATERIALTYPE
} MATERIALType;

class Rat_PeriodontalLigament : public Fem{
 public:
  void postProcess_LinearElastic_element_spatialForm(const int ic,const bool option);
  void calcStressTensor_LinearElastic_element_spatialForm(const int ic,const bool option);
  void inputMaterialInfo(TextParser &tp);

 private:
  double LinearElastic_inGaussIntegral(ARRAY2D<double> &dNdr,ARRAY2D<double> &x_ref,const int numOfNodeInElm,const double weight,const int ic,const double lambda,const double mu,const bool option);
  double postProcess_LinearElastic_inGaussIntegral(double (&sigmaEigen)[3],double (&sigmaEigenVector)[3][3],ARRAY2D<double> &u,ARRAY2D<double> &dNdr,ARRAY2D<double> &x_ref,ARRAY2D<double> &dNdx,const int numOfNodeInElm,const double weight,const int ic,const double lambda,const double mu,const bool option);

};

#endif