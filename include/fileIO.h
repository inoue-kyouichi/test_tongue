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
#include <string>

class fileIO{
 public:
  static int CountNumbersOfTextLines(const std::string &filePath);
  static void read_geometry_node(DOUBLEARRAY2 &x, int &numOfNode,const std::string &file);
  static void read_geometry_meshType(elementType &element,int &numOfElm,const std::string &file);
  static void read_geometry_element(elementType &element,const int &numOfElm,const std::string &file);
  static void export_vtu(const DOUBLEARRAY2 &x,const elementType &element,const int &numOfNode,const int &numOfElm,
              const DOUBLEARRAY2 &U,const DOUBLEARRAY1 &volumeChangeRatio,const std::string &file);
  static void export_vtu(const DOUBLEARRAY2 &x,const elementType &element,
            const int &numOfNode,const int &numOfElm,
            const DOUBLEARRAY2 &U,const DOUBLEARRAY1 &volumeChangeRatio,const DOUBLEARRAY2 &lambda_ave,
            const DOUBLEARRAY2 &sigmaEigen_ave,const DOUBLEARRAY2 &AEigen_ave,
            const DOUBLEARRAY3 &sigmaEigenVector_ave,const DOUBLEARRAY3 &AEigenVector_ave,const DOUBLEARRAY2 &innerForce, const std::string &file);
  static void export_vtu(const DOUBLEARRAY2 &x,const elementType &element,
            const int &numOfNode,const int &numOfElm,
            const DOUBLEARRAY2 &U,const DOUBLEARRAY1 &volumeChangeRatio,const DOUBLEARRAY2 &lambda_ave,
            const DOUBLEARRAY1 &bundle,const INTARRAY1 &bundleElement,
            const DOUBLEARRAY2 &sigmaEigen_ave,const DOUBLEARRAY2 &AEigen_ave,
            const DOUBLEARRAY3 &sigmaEigenVector_ave,const DOUBLEARRAY3 &AEigenVector_ave,const DOUBLEARRAY2 &innerForce, const std::string &file);
  static void export_vtu(const DOUBLEARRAY2 &x,const elementType &element,const int &numOfNode,const int &numOfElm,
              const DOUBLEARRAY2 &U,const DOUBLEARRAY1 &volumeChangeRatio,const DOUBLEARRAY2 &lambda_ave,const std::string &file);
  static void export_vtu_boundary(const DOUBLEARRAY2 &x,const elementType &element,
            const int &numOfNode,const int &numOfElm,
            const INTARRAY2 &ibd,const DOUBLEARRAY2 &bd,const DOUBLEARRAY2 &fiberDirection_elm,const std::string &file);

};

#endif //_FILE_IO_H_
