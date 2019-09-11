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
  static void read_geometry_node(DOUBLEARRAY2D &x, int &numOfNode,const std::string &file);
  static void read_geometry_meshType(elementType &element,int &numOfElm,const std::string &file);
  static void read_geometry_element(elementType &element,const int &numOfElm,const std::string &file);
  static void export_vtu(DOUBLEARRAY2D &x,const elementType &element,const int &numOfNode,const int &numOfElm,
              DOUBLEARRAY2D &U,DOUBLEARRAY1D &volumeChangeRatio,const std::string &file);
  static void export_vtu(DOUBLEARRAY2D &x,const elementType &element,
            const int &numOfNode,const int &numOfElm,
            DOUBLEARRAY2D &U,DOUBLEARRAY1D &volumeChangeRatio,DOUBLEARRAY2D &lambda_ave,
            DOUBLEARRAY2D &sigmaEigen_ave,DOUBLEARRAY2D &AEigen_ave,
            DOUBLEARRAY3D &sigmaEigenVector_ave,DOUBLEARRAY3D &AEigenVector_ave,DOUBLEARRAY2D &innerForce, const std::string &file);
  static void export_vtu(DOUBLEARRAY2D &x,const elementType &element,
            const int &numOfNode,const int &numOfElm,
            DOUBLEARRAY2D &U,DOUBLEARRAY1D &volumeChangeRatio,DOUBLEARRAY2D &lambda_ave,
            DOUBLEARRAY2D &sigmaEigen_ave,DOUBLEARRAY3D &sigmaEigenVector_ave,
            const std::string &file);
  static void export_vtu(DOUBLEARRAY2D &x,const elementType &element,
            const int &numOfNode,const int &numOfElm,
            DOUBLEARRAY2D &U,DOUBLEARRAY1D &volumeChangeRatio,DOUBLEARRAY2D &lambda_ave,
            DOUBLEARRAY1D &bundle,INTARRAY1D &bundleElement,
            DOUBLEARRAY2D &sigmaEigen_ave,DOUBLEARRAY2D &AEigen_ave,
            DOUBLEARRAY3D &sigmaEigenVector_ave,DOUBLEARRAY3D &AEigenVector_ave,DOUBLEARRAY2D &innerForce, const std::string &file);
  static void export_vtu(DOUBLEARRAY2D &x,const elementType &element,const int &numOfNode,const int &numOfElm,
              DOUBLEARRAY2D &U,DOUBLEARRAY1D &volumeChangeRatio,DOUBLEARRAY2D &lambda_ave,const std::string &file);
  static void export_vtu_boundary(DOUBLEARRAY2D &x,const elementType &element,
            const int &numOfNode,const int &numOfElm,
            INTARRAY2D &ibd,DOUBLEARRAY2D &bd,DOUBLEARRAY2D &fiberDirection_elm,const std::string &file);

};

#endif //_FILE_IO_H_
