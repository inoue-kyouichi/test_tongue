#ifndef _RIGIDBODY_H_
#define _RIGIDBODY_H_

//##################################################################################
//
// rigidBody.h
//
// Copyright (c) 2018 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   rididBody.h
 * @brief  rigid body class Header
 * @author T.Otani
 */
#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>

#include "allocation.h"
#include "math_tools.h"

class TriangleSet{
 public:
  int numOfNode,numOfElm;
  INTARRAY2 elm;
  DOUBLEARRAY2 x,x0;
  void readPLY(const std::string &file);
  void translation(const double (&center)[3]);
};

class RigidBody{
 public:
  TriangleSet obj;
  double rho;
  double Mass,J[3][3],coordinates[3][3],coordinates_ref[3][3];
  double xg[3],U[3],w[3];
  double R[3][3];
  void initialize(TextParser &tp);
  void calcMassProperties(const TriangleSet &obj);
  void updateRotationMatrix_spatialForm(const double (&w)[3]);
  void exportPLY(const std::string &file);
  void updateShape();
 private:
};

#endif