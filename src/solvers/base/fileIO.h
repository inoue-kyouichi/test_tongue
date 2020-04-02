#ifndef _FILE_IO_H_
#define _FILE_IO_H_

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
 * @file   fileIO.h
 * @brief  FEMBase Definition Header
 * @author T. Otani
 */

#include "fem_define.h"
#include "allocation.h"
#include <string>

class fileIO{
 public:
  static int CountNumbersOfTextLines(const std::string &filePath);
  static void read_geometry_node(ARRAY2D<double> &x, int &numOfNode,const std::string &file);
  static void read_geometry_meshType(std::vector<ElementType> &element,int &numOfElm,const std::string &file);
  static void read_geometry_materialType(std::vector<ElementType> &element,int &numOfElm,const std::string &file);
  static void read_geometry_element(std::vector<ElementType> &element,const int &numOfElm,const std::string &file);
  static void export_vtu(ARRAY2D<double> &x,const std::vector<ElementType> &element,const int &numOfNode,const int &numOfElm,ARRAY2D<double> &U,const std::string &file);
  static void export_vtu(ARRAY2D<double> &x,const std::vector<ElementType> &element,const int &numOfNode,const int &numOfElm,
              ARRAY2D<double> &U,ARRAY1D<double> &volumeChangeRatio,const std::string &file);
  static void export_vtu(ARRAY2D<double> &x,const std::vector<ElementType> &element,
            const int &numOfNode,const int &numOfElm,
            ARRAY2D<double> &U,ARRAY1D<double> &volumeChangeRatio,ARRAY2D<double> &lambda_ave,
            ARRAY2D<double> &sigmaEigen_ave,ARRAY2D<double> &AEigen_ave,
            ARRAY3D<double> &sigmaEigenVector_ave,ARRAY3D<double> &AEigenVector_ave,ARRAY2D<double> &innerForce, const std::string &file);
  static void export_vtu(ARRAY2D<double> &x,const std::vector<ElementType> &element,
            const int &numOfNode,const int &numOfElm,
            ARRAY2D<double> &U,ARRAY1D<double> &volumeChangeRatio,ARRAY2D<double> &lambda_ave,
            ARRAY2D<double> &sigmaEigen_ave,ARRAY3D<double> &sigmaEigenVector_ave,
            const std::string &file);
  static void export_vtu(ARRAY2D<double> &x,const std::vector<ElementType> &element,
            const int &numOfNode,const int &numOfElm,
            ARRAY2D<double> &U,ARRAY1D<double> &volumeChangeRatio,ARRAY2D<double> &lambda_ave,
            ARRAY1D<double> &bundle,ARRAY1D<int> &bundleElement,
            ARRAY2D<double> &sigmaEigen_ave,ARRAY2D<double> &AEigen_ave,
            ARRAY3D<double> &sigmaEigenVector_ave,ARRAY3D<double> &AEigenVector_ave,ARRAY2D<double> &innerForce, const std::string &file);
  static void export_vtu(ARRAY2D<double> &x,const std::vector<ElementType> &element,const int &numOfNode,const int &numOfElm,
              ARRAY2D<double> &U,ARRAY1D<double> &volumeChangeRatio,ARRAY2D<double> &lambda_ave,const std::string &file);
  static void export_vtu_boundary(ARRAY2D<double> &x,const std::vector<ElementType> &element,
            const int &numOfNode,const int &numOfElm,
            ARRAY2D<int> &ibd,ARRAY2D<double> &bd,ARRAY2D<double> &fiberDirection_elm,const std::string &file);
  static void export_vtu_boundary(ARRAY2D<double> &x,const std::vector<ElementType> &element,
            const int &numOfNode,const int &numOfElm,
            ARRAY2D<int> &ibd,ARRAY2D<double> &bd,const std::string &file);

  static void export_vtu_Mises(ARRAY2D<double> &x,const std::vector<ElementType> &element,const int &numOfNode,const int &numOfElm,ARRAY2D<double> &U,ARRAY1D<double> &Mises,const std::string &file);
};

#endif //_FILE_IO_H_
