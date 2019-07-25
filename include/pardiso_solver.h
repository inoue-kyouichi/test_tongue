#ifndef _PARDISO_SOLVER_H_
#define _PARDISO_SOLVER_H_

//##################################################################################
//
// LIS solver
//
// Copyright (c) 2018 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   PARDISO_solver.h
 * @brief  FEMBase Definition Header
 * @author T. Otani
 */

#include <string>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "fem_define.h"

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

class PARDISO_solver{
 public:
    PARDISO_solver(){
    };
    ~PARDISO_solver(){
      free(ptr);
      free(index);
      free(value);
      free(b);
      free(x);
    };

  //CSR form
  int nnz;
  MKL_INT *ptr,*index;
  double *value,*b,*x;

  void main(MKL_INT n,const int numOfOMP);
  void initialize(const int &DOF);
  void CSR_initialize(const INTVECTOR2 &inb,const int &numOfNode,const int &dim);

  double vector_norm(const int &nump,const double *x);

  void set_CSR_value(const DOUBLEARRAY5 &K,const elementType &element,const int &numOfNode,
                               const int &numOfElm,const INTVECTOR2 &inb);
  void set_CSR_dirichlet_boundary_condition(const int &numOfNode,const INTARRAY2 &ibd);

private:
  void CSR_ptr_initialize(const INTVECTOR2 &inb,const int &numOfNode,const int &dim);
  void CSR_index_initialize(const INTVECTOR2 &inb,const int &numOfNode,const int &dim);

 public:
  //rigid body interaction
  void initialize(const int &DOF,const int &DOF2);
  void CSR_initialize(const INTVECTOR2 &inb,const int &numOfNode,const INTARRAY1 &iCP,const INTARRAY1 &CP,const int &numOfCP,const int &dim);
  void set_CSR_value_rigidBodyInteraction(const int &numOfNode,const INTARRAY1 &iCP,const DOUBLEARRAY3 &Rb,const double (&Kqq)[3][3],const int numOfCP);
 private:
  //rigid body interaction
  void CSR_ptr_initialize(const INTVECTOR2 &inb,const int &numOfNode,const INTARRAY1 &iCP,const int &numOfCP,const int &dim);
  void CSR_index_initialize(const INTVECTOR2 &inb,const int &numOfNode,const INTARRAY1 &iCP,const INTARRAY1 &CP,const int &numOfCP,const int &dim);

};

#endif //_PARDISO_SOLVER_H_
