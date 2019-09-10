#ifndef _DOMAIN_H_
#define _DOMAIN_H_

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
 * @file   domain.h
 * @brief  FEMBase Definition Header
 * @author T. Otani
 */

#include <string>
#include "fem_define.h"

class Domain{
public:
  int numOfNode, numOfElm, numOfDirichlet,numOfNeumann;
  int meshType;
  DOUBLEARRAY2 x,x0;
  elementType element,belement;
  INTARRAY2 ibd;
  DOUBLEARRAY2 bd;
  DOUBLEARRAY2 bn;
  INTVECTOR2 inb;
  INTVECTOR2 ieb;

  void set_geometry(const std::string &file1,const std::string &file2,const std::string &file3);
  void set_dirichlet(const std::string &D_file);
  void set_dirichlet(const std::string &D_file,const std::string &Dvalue_file);
  void set_neumann(const std::string &N_file,const std::string &meshTypeFile,const std::string &NvalueFile);
  void export_vtk(const std::string &file);

private:
  void calc_adjacent_nodes();
  void calc_adjacent_elements();
private:

};

#endif //_DOMAIN_H_