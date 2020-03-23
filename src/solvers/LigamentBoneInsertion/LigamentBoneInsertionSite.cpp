/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "LigamentBoneInsertion.h"
#include <ostream>
#include <fstream>

using namespace std;

// #################################################################
/**
 * @brief fem solid analysis routine
 * @param [inout] PARDISO   PARDISO class
 */
void InsertionSite::femSolidAnalysis()
{
  int output_iter=1;
  int increment=20;
  string output;

  if(Restart==0){
    output = outputDir + "/test_" + to_string(0) + ".vtu";
    fileIO::export_vtu(x,element,numOfNode,numOfElm,U,volumeChangeRatio,output);
  }

  // for(int ic=0;ic<numOfElm;ic++) calcVolume_hexa(ic,volume0,8,2,0);

  for(int loop=Restart+1;loop<=increment;loop++){

    for(int i=0;i<numOfNode;i++){
      for(int j=0;j<3;j++){
        if(ibd(i,j)==0) U(i,j)+=bd(i,j)/(double)increment;
      }
    }

    if(NRscheme()==true){
      printf("NaN detected in loop=%d. Exit...\n",loop);
      exit(1);
    }

    #pragma omp parallel for
    for(int i=0;i<numOfNode;i++){
      for(int j=0;j<3;j++) x(i,j) = x0(i,j) + U(i,j);
    }

    FILE *fp;
    output = "Restart/U_" + to_string(loop) + ".dat";
    if ((fp = fopen(output.c_str(), "w")) == NULL) {
      cout << "file open error" << endl;
      exit(1);
    }
    for(int i=0;i<numOfNode;i++){
      fprintf(fp,"%e %e %e\n",U(i,0),U(i,1),U(i,2));
    }
    fclose(fp);
    if(loop%output_iter==0){
      output = outputDir + "/test_" + to_string(loop/output_iter) + ".vtu";
      fileIO::export_vtu(x,element,numOfNode,numOfElm,U,volumeChangeRatio,lambda_ave,output);
    }
  }
}

// #################################################################
/**
 * @brief fem solid analysis routine
 * @param [inout] PARDISO   PARDISO class
 */
bool InsertionSite::NRscheme()
{
    double residual,residual0,norm,norm0;
    string output;

    for(int ic=1;ic<=NRiteration;ic++){

      calcStressTensor();  //calc K and Q
      set_rhs_statics();

      PARDISO.set_CSR_value3D(Ku,element,numOfNode,numOfElm,inb);
      PARDISO.set_CSR_dirichlet_boundary_condition3D(numOfNode,ibd);

      for(int i=0;i<numOfNode;i++){
        for(int j=0;j<3;j++) PARDISO.b[i+j*numOfNode]=RHS(i,j);
      }

      if(ic==1) residual0 = PARDISO.vector_norm(numOfNode*3,PARDISO.b);

      PARDISO.main(3*numOfNode,OMPnumThreads);
      norm = PARDISO.vector_norm(numOfNode*3,PARDISO.x);
      if(ic==1) norm0 = norm;
      // corrector_statistics(PARDISO.x,1e0);
      // residual=line_search(PARDISO.x);

      // if(isnan(residual)){
      //   cout << "residual is nan. Exit..." << endl;
      //   return true;
      // }

      // for(int i=0;i<numOfElm;i++){
      //   calcVolume_hexa(i,volume,8,2,1);
      //   volumeChangeRatio[i]=volume[i]/volume0[i];
      // }
      for(int i=0;i<numOfNode;i++){
        for(int j=0;j<3;j++) x(i,j) = x0(i,j) + U(i,j);
      }
      // output = outputDir + "/test_NR_" + to_string(ic) + ".vtu";
      // fileIO::export_vtu(x,element,numOfNode,numOfElm,U,volumeChangeRatio,lambda_ave,output);

      cout << "NewtonRaphson_iteration = " << ic << endl;

      if(residual/residual0<NRtolerance) break;
      // if(test!=1 && ic>50) break;
    }
    return false;
}

// #################################################################
/**
 * @brief calc stress tensor
 */
void InsertionSite::calcStressTensor()
{
  double elementVolume,volume;

  stress_tensor_initialize();

  volume = 0e0;

  #pragma omp parallel for
  for(int ic=0;ic<numOfElm;ic++){
    calcStressTensor_SantVenant_element_spatialForm(ic,U,true);

  //   switch(element[ic].meshType){
  //     case VTK_HEXAHEDRON:
  //       //calcStressTensor_ACL_element_spatialForm_hexa_Fbar(ic,U,8,2,true);
  //       // calcStressTensor_NeoHookean_element_spatialForm_hexa_Fbar(ic,U,8,2,true);
  //       // calcStressTensor_NeoHookean_element_spatialForm_hexa(ic,U,8,2,true);
  //       calcStressTensor_SantVenant_element_spatialForm_hexa(ic,U,8,2,true);
  //       break;
  //     case VTK_QUADRATIC_HEXAHEDRON:
  //       break;
  //     case VTK_WEDGE:
  //       calcStressTensor_SantVenant_element_spatialForm_prism(ic,U,6,2,true);
  //       break;
  //     default:
  //       cout << "undefined mesh type" << endl;
  //       exit(1);
  //     }
  //   volume += elementVolume;
  }

  totalVolume = volume;

  for(int ic=0;ic<numOfElm;ic++){
    for(int p=0;p<element[ic].node.size();p++){
      for(int i=0;i<3;i++){
        innerForce(element[ic].node[p],i) += Qu[ic](p,i);
      }
    }
  }
}

// #################################################################
/**
 * @brief assembly right-hand side
 */
void InsertionSite::set_rhs_statics()
{
 #pragma omp parallel for
  for(int ic=0;ic<numOfNode;ic++){
    for(int j=0;j<3;j++){
      RHS(ic,j) = (double)ibd(ic,j)*( - innerForce(ic,j));
    }
  }
}
