//##################################################################################
//
// FEM solid analysis
//
// Copyright (c) 2018 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Enginieering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   pardiso_solver_rigidBodyInteraction.h
 * @author T. Otani
 */

#include "pardiso_solver.h"

using namespace std;

// #################################################################
/**
 * @brief initialize PARDISO class (rigid Body interaction)
 * @param [in] DOF  total DOF
 */
void PARDISO_solver::initialize(const int &DOF,const int &DOF2)
{
  int numOfDOF=DOF+DOF2+3+3;
  b=(double *)malloc(numOfDOF*sizeof(double));
  x=(double *)malloc(numOfDOF*sizeof(double));
}
// #################################################################
/**
 * @brief set CSR matrix format
 * @param [in] inb          numbers of adjacent nodes of each node
 * @param [in] numOfNode    number of nodes
 * @param [in] dim          total DOF
 */
void PARDISO_solver::CSR_initialize(const INTVECTOR2 &inb,const int &numOfNode,INTARRAY1D &iCP,INTARRAY1D &CP,const int &numOfCP,const int &dim)
{
  ptr = (MKL_INT *)malloc((numOfNode*dim+numOfCP*dim+6+1)*sizeof(MKL_INT));

  CSR_ptr_initialize(inb,numOfNode,iCP,numOfCP,dim);
  index = (MKL_INT *)malloc( nnz*sizeof(MKL_INT) );
  value = (double *)malloc( nnz*sizeof(double));

  CSR_index_initialize(inb,numOfNode,iCP,CP,numOfCP,dim);

  // for(int ic=0;ic<numOfNode*dim+numOfCP*dim+6;ic++){
  //   printf("%d %d %d:",ic,ptr[ic],ptr[ic+1]);
  //   for(int i=ptr[ic];i<ptr[ic+1];i++){
  //     printf("%d ",index[i]);
  //   }
  //   printf("\n");
  // }

}

// #################################################################
/**
 * @brief set ptr (CSR matrix) format
 * @param [in] inb          numbers of adjacent nodes of each node
 * @param [in] numOfNode    number of nodes
 * @param [in] dim          total DOF
 */
void PARDISO_solver::CSR_ptr_initialize(const INTVECTOR2 &inb,const int &numOfNode,INTARRAY1D &iCP,const int &numOfCP,const int &dim)
{
  nnz = 0;

  //delta u
  for(int i=0;i<dim;i++){
    for(int ic=0;ic<numOfNode;ic++){
      ptr[ic+i*numOfNode] = nnz;
      nnz += inb[ic].size()*dim;
      if(iCP(ic)!=-1) nnz += 1;
    }
  }
  //delta lambda
  for(int i=0;i<dim;i++){
    for(int ic=0;ic<numOfCP;ic++){
      ptr[numOfNode*dim+ic+i*numOfCP] = nnz;
      nnz += 1+1+3;
    }
  }
  //delta U
  for(int i=0;i<dim;i++){
    ptr[numOfNode*dim+numOfCP*dim+i] = nnz;
    nnz += numOfCP;
  }
  //delta w
  for(int i=0;i<dim;i++){
    ptr[numOfNode*dim+numOfCP*dim+3+i] = nnz;
    nnz += numOfCP*dim+3;
  }

  ptr[numOfNode*dim+numOfCP*dim+3+3] = nnz;
}

// #################################################################
/**
 * @brief set index (CSR matrix) format
 * @param [in] inb          numbers of adjacent nodes of each node
 * @param [in] numOfNode    number of nodes
 * @param [in] dim          total DOF
 */
void PARDISO_solver::CSR_index_initialize(const INTVECTOR2 &inb,const int &numOfNode,INTARRAY1D &iCP,INTARRAY1D &CP,const int &numOfCP,const int &dim)
{
  int tmp = 0;

  //delta u
  for(int dim2=0;dim2<dim;dim2++){
    for(int ic=0;ic<numOfNode;ic++){

      for(int k=0;k<dim;k++){
        for(int i=0;i<inb[ic].size();i++){
          index[tmp] = inb[ic][i]+k*numOfNode;
          tmp++;
        }
      }
      if(iCP(ic)!=-1){
        index[tmp] = dim*numOfNode+iCP(ic)+dim2*numOfCP;
        tmp++;
      }

    }
  }

  //delta lambda
  for(int dim2=0;dim2<dim;dim2++){
    for(int ic=0;ic<numOfCP;ic++){
      index[tmp] = CP(ic)+dim2*numOfNode;
      tmp++;
      index[tmp] = dim*numOfNode+dim*numOfCP+dim2;
      tmp++;
      for(int k=0;k<3;k++){
        index[tmp] = dim*numOfNode+dim*numOfCP+3+k;
        tmp++;
      }
    }
  }

  //delta U
  for(int k=0;k<dim;k++){
    for(int ic=0;ic<numOfCP;ic++){
      index[tmp] = dim*numOfNode+iCP(CP(ic))+k*numOfCP;
      tmp++;
    }
  }

  //delta w
  for(int dim2=0;dim2<dim;dim2++){
    for(int k=0;k<dim;k++){
      for(int ic=0;ic<numOfCP;ic++){
        index[tmp] = dim*numOfNode+iCP(CP(ic))+k*numOfCP;
        tmp++;
      }
    }
    for(int k=0;k<dim;k++){
        index[tmp] = dim*numOfNode+dim*numOfCP+3+k;
        tmp++;
    }
  }
}

// #################################################################
/**
 * @brief calc set CSR matrix (Vector)
 * @param [in] K element stiffness matrix
 * @param [in] ie elements
 * @param [in] numOfNode number of nodes
 * @param [in] numOfElm number of elements
 * @param [in] numOfNodeInElm nubmer of node in each elements
 * @param [in] inb nodes around each node
 */
void PARDISO_solver::set_CSR_value_rigidBodyInteraction(const int &numOfNode,INTARRAY1D &iCP,DOUBLEARRAY3D &Rb,const double (&Kqq)[3][3],const int numOfCP)
{
  int tmp1,tmp2,tmp3;

  for(int k=0;k<3;k++){
    for(int ic=0;ic<numOfNode;ic++){
      if(iCP(ic)==-1) continue;
      for(int i=ptr[ic+k*numOfNode];i<ptr[ic+k*numOfNode+1];i++){
        if(index[i]>=numOfNode*3) value[i]=1e0;
      }
    }
  }

  for(int k=0;k<3;k++){
    for(int ic=0;ic<numOfCP;ic++){
      for(int i=ptr[3*numOfNode+ic+k*numOfCP];i<ptr[3*numOfNode+ic+k*numOfCP+1];i++){
        if(index[i]<numOfNode*3){
          value[i]=1e0;
        }else if(index[i]<numOfNode*3+numOfCP*3+3){
          value[i]=-1e0;
        }else if(index[i]==numOfNode*3+numOfCP*3+3){
          value[i]=Rb(ic,k,0);
        }else if(index[i]==numOfNode*3+numOfCP*3+3+1){
          value[i]=Rb(ic,k,1);
        }else if(index[i]==numOfNode*3+numOfCP*3+3+2){
          value[i]=Rb(ic,k,2);
        }
      }
    }
  }

  for(int ic=3*numOfNode+3*numOfCP;ic<3*numOfNode+3*numOfCP+3;ic++){
    for(int i=ptr[ic];i<ptr[ic+1];i++){
      value[i]=-1e0;
    }
  }

  for(int k=0;k<3;k++){
    tmp1=ptr[3*numOfNode+3*numOfCP+3+k];
    for(int ic=0;ic<numOfCP;ic++){
      value[tmp1+ic]          =-Rb(ic,k,0);
      value[tmp1+numOfCP+ic]  =-Rb(ic,k,1);
      value[tmp1+numOfCP*2+ic]=-Rb(ic,k,2);
    }
  }

  for(int k=0;k<3;k++){
    for(int i=ptr[3*numOfNode+3*numOfCP+3];i<ptr[3*numOfNode+3*numOfCP+6];i++){
      if(index[i]==numOfNode*3+numOfCP*3+3){
        value[i]=Kqq[k][0];
      }else if(index[i]==numOfNode*3+numOfCP*3+4){
        value[i]=Kqq[k][1];
      }else if(index[i]==numOfNode*3+numOfCP*3+5){
        value[i]=Kqq[k][2];
      }
    }
  }

  // for(int ic=0;ic<numOfNode*3+numOfCP*3+6;ic++){
  //   printf("%d %d %d:",ic,ptr[ic],ptr[ic+1]);
  //   for(int i=ptr[ic];i<ptr[ic+1];i++){
  //     printf("%e ",value[i]);
  //   }
  //   printf("\n");
  // }
  // exit(1);
}
