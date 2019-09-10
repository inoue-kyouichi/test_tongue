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
#include "allocation.h"

#define FEM_VERS "1.0"
#define GRAVITY     9.80665 // (m/s2)
static const double PI = 3.1415926535897932384626;

typedef std::vector<int> INTVECTOR1;
typedef std::vector<std::vector<int>> INTVECTOR2;
typedef std::vector<std::vector<std::vector<int>>> INTVECTOR3;
typedef std::vector<std::vector<std::vector<std::vector<int>>>> INTVECTOR4;
typedef std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> INTVECTOR5;

typedef std::vector<double> DOUBLEVECTOR1;
typedef std::vector<std::vector<double>> DOUBLEVECTOR2;
typedef std::vector<std::vector<std::vector<double>>> DOUBLEVECTOR3;
typedef std::vector<std::vector<std::vector<std::vector<double>>>> DOUBLEVECTOR4;
typedef std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> DOUBLEVECTOR5;
typedef std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> DOUBLEVECTOR6;
typedef std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>> DOUBLEVECTOR7;

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


class Element{
 public:
  VTKCellType meshType;
  INTVECTOR1 node;
};

typedef std::vector<Element> elementType;


#endif // _FB_DEFINE_H_
