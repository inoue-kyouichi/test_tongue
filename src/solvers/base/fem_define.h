#ifndef _FEM_DEFINE_H_
#define _FEM_DEFINE_H_

//##################################################################################
//
// FEM Base
//
// Copyright (c) 2016 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   fem_define.h
 * @brief  FEMBase Definition Header
 * @author T. Otani
 */
#include <string>
#include <vector>
#include "TextParser.h"
#include "vtkCellType.h"

template<typename T>
using VECTOR1D = std::vector<T>;
template<typename T>
using VECTOR2D = std::vector<std::vector<T>>;
template<typename T>
using VECTOR3D = std::vector<std::vector<std::vector<T>>>;
template<typename T>
using VECTOR4D = std::vector<std::vector<std::vector<std::vector<T>>>>;
template<typename T>
using VECTOR5D = std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>>;
template<typename T>
using VECTOR6D = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>>>;
template<typename T>
using VECTOR7D = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>>>>;

// typedef std::vector<int> INTVECTOR1;
// typedef std::vector<std::vector<int>> INTVECTOR2;
// typedef std::vector<std::vector<std::vector<int>>> INTVECTOR3;
// typedef std::vector<std::vector<std::vector<std::vector<int>>>> INTVECTOR4;
// typedef std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> INTVECTOR5;

// typedef std::vector<double> DOUBLEVECTOR1;
// typedef std::vector<std::vector<double>> DOUBLEVECTOR2;
// typedef std::vector<std::vector<std::vector<double>>> DOUBLEVECTOR3;
// typedef std::vector<std::vector<std::vector<std::vector<double>>>> DOUBLEVECTOR4;
// typedef std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> DOUBLEVECTOR5;
// typedef std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> DOUBLEVECTOR6;
// typedef std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>> DOUBLEVECTOR7;

#define FEM_VERS "1.0"
constexpr double GRAVITY = 9.80665; // (m/s2)
constexpr double PI = 3.1415926535897932384626;

// general
#define ON          1
#define OFF         0

//data encode
#define INT         0
#define DOUBLE      1
#define ASCII       0
#define BINARY      1

// IO file format
#define ORIGINAL    0
#define HDF5        1
#define VTK         2
#define VTU         3

//diciclet boundary condiiton
#define WALL     0

typedef enum {
  PULP       = 0,
  PDL        = 1,
  NUMBER_OF_MATERIALTYPE
} MATERIALType;

class Element{
 public:
  VTKCellType meshType;
  int materialType;
  int numOfGaussPoint;
  //MATERIALType  materialType;
  VECTOR1D<int> node;
};

using elementType = std::vector<Element>;


#endif // _FB_DEFINE_H_
