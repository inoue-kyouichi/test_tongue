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
#include "vtkCellType.h"
#include "allocation.h"

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

#define FEM_VERS "1.0"
constexpr double GRAVITY = 9.80665; // (m/s2)
constexpr double PI = 3.1415926535897932384626;



class ElementType{
 public:
  VTKCellType meshType;
  int materialType;
  int numOfGaussPoint;
  VECTOR1D<int> node;
  VECTOR1D<int> face;
};


#endif // _FB_DEFINE_H_
