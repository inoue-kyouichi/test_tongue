/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "SFEM.h"
#include <ostream>
#include <fstream>

using namespace std;


// #################################################################
/**
 * @brief fem solid analysis routine
 * @param [inout] PARDISO   PARDISO class
 */
void SmoothedFEM::SFEM::femSolidAnalysis_linear()
{
  string output;

  calcStressTensor();  //calc K and Q
  calcBoundaryForce();
  set_rhs_statics();

  PARDISO.set_CSR_value3D(Ku,element,numOfNode,numOfElm,inb);
  PARDISO.set_CSR_dirichlet_boundary_condition3D(numOfNode,ibd);

  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) PARDISO.b[i+j*numOfNode]=RHS(i,j);
  }

  PARDISO.main(3*numOfNode,OMPnumThreads);
  corrector_statics(PARDISO.x,1e0);

  #pragma omp parallel for
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) x(i,j) = x0(i,j) + U(i,j);
  }

  output = outputDir + "/test.vtu";
  export_vtu(output);

}

// #################################################################
/**
 * @brief fem solid analysis routine
 * @param [inout] PARDISO   PARDISO class
 */
void SmoothedFEM::SFEM::femSolidAnalysis()
{
  int output_iter=1;
  string output;

  if(Restart==0){
    output = outputDir + "/test_" + to_string(0) + ".vtu";
    fileIO::export_vtu(x,element,numOfNode,numOfElm,U,volumeChangeRatio,output);
  }

  // for(int ic=0;ic<numOfElm;ic++) calcVolume_hexa(ic,volume0,8,2,0);

  for(int loop=Restart+1;loop<=maxIteration;loop++){

    if(NRscheme()==true){
      printf("NaN detected in loop=%d. Exit...\n",loop);
      exit(1);
    }

    #pragma omp parallel for
    for(int i=0;i<numOfNode;i++){
      for(int j=0;j<3;j++) x(i,j) = x0(i,j) + U(i,j);
    }

    // FILE *fp;
    // output = "Restart/U_" + to_string(loop) + ".dat";
    // if ((fp = fopen(output.c_str(), "w")) == NULL) {
    //   cout << "file open error" << endl;
    //   exit(1);
    // }
    // for(int i=0;i<numOfNode;i++){
    //   fprintf(fp,"%e %e %e\n",U(i,0),U(i,1),U(i,2));
    // }
    // fclose(fp);

    // if(loop%output_iter==0){
    //   output = outputDir + "/test_" + to_string(loop/output_iter) + ".vtu";
    //   fileIO::export_vtu(x,element,numOfNode,numOfElm,U,volumeChangeRatio,lambda_ave,output);
    // }
  }

  output = outputDir + "/test.vtu";
  export_vtu(output);

}

// #################################################################
/**
 * @brief fem solid analysis routine
 * @param [inout] PARDISO   PARDISO class
 */
bool SmoothedFEM::SFEM::NRscheme()
{
    double residual,residual0,norm,norm0;
    string output;

    for(int ic=1;ic<=NRiteration;ic++){

      calcStressTensor();  //calc K and Q
      calcBoundaryForce();
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
      corrector_statics(PARDISO.x,relaxation);
      // residual=line_search(PARDISO.x);

      if(isnan(norm)){
        cout << "residual is nan. Exit..." << endl;
        return true;
      }

      // for(int i=0;i<numOfElm;i++){
      //   calcVolume_hexa(i,volume,8,2,1);
      //   volumeChangeRatio[i]=volume[i]/volume0[i];
      // }
      for(int i=0;i<numOfNode;i++){
        for(int j=0;j<3;j++) x(i,j) = x0(i,j) + U(i,j);
      }
      // output = outputDir + "/test_NR_" + to_string(ic) + ".vtu";
      // fileIO::export_vtu(x,element,numOfNode,numOfElm,U,volumeChangeRatio,lambda_ave,output);

      printf("NR iter.=%d norm/norm0=%e\n",ic,norm/norm0);

      if(norm/norm0<NRtolerance) break;
      // if(test!=1 && ic>50) break;
    }
    return false;
}

// #################################################################
/**
 * @brief calc stress tensor
 */
void SmoothedFEM::SFEM::calcStressTensor()
{
  double elementVolume,volume;

  stress_tensor_initialize();

  volume = 0e0;

  #pragma omp parallel for
  for(int ic=0;ic<numOfElm;ic++){
        calcStressTensor_linearElasticMaterial_element(ic,true);
        // calcStressTensor_SantVenant_element_spatialForm(ic,U,true);
        // calcStressTensor_NeoHookean_element_spatialForm(ic,100e3,0.49e0,U,true);
  }

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
void SmoothedFEM::SFEM::set_rhs_statics()
{
 #pragma omp parallel for
  for(int ic=0;ic<numOfNode;ic++){
    for(int j=0;j<3;j++){
      RHS(ic,j) = (double)ibd(ic,j)*( - innerForce(ic,j)+boundaryForce(ic,j));
    }
  }
}


// #################################################################
/**
 * @brief calc boundary conditions
 * @param [in] stress
 */
void SmoothedFEM::SFEM::export_vtu(const string &file)
{
  FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL) {
    cout << file << " open error" << endl;
    exit(1); 
  }

  fprintf(fp,"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<UnstructuredGrid>\n");
  fprintf(fp,"<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n",numOfNode,numOfElm);
  fprintf(fp,"<Points>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",x(i,0),x(i,1),x(i,2));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Points>\n");
  fprintf(fp,"<Cells>\n");
  fprintf(fp,"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    for(int j=0;j<element[i].node.size();j++) fprintf(fp,"%d ",element[i].node[j]);
    fprintf(fp,"\n");
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num=0;
  for(int i=0;i<numOfElm;i++){
    num += element[i].node.size();
    fprintf(fp,"%d\n",num);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++) fprintf(fp,"%d\n",element[i].meshType);
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Cells>\n");

  fprintf(fp,"<PointData Vectors=\"displacement[m/s]\">\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"displacement[m/s]\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",U(i,0),U(i,1),U(i,2));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");

  fprintf(fp,"<CellData>");
  // fprintf(fp,"<DataArray type=\"Int64\" Name=\"Material\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  // for(int i=0;i<numOfElm;i++){
  //   fprintf(fp,"%d\n",element[i].materialType);
  // }
  // fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");
  fprintf(fp,"</Piece>");
  fprintf(fp,"</UnstructuredGrid>");
  fprintf(fp,"</VTKFile>");
  fclose(fp);
}

// #################################################################
/**
 * @brief set element face information
 */
void SmoothedFEM::SFEM::setFace()
{
  auto facesInElement=[](ARRAY3D<int> &tmp,std::vector<ElementType> &element,const int ic){
    tmp(ic,0,0) = element[ic].node[0]; tmp(ic,0,1) = element[ic].node[2]; tmp(ic,0,2) = element[ic].node[1];
    tmp(ic,1,0) = element[ic].node[0]; tmp(ic,1,1) = element[ic].node[1]; tmp(ic,1,2) = element[ic].node[3];
    tmp(ic,2,0) = element[ic].node[1]; tmp(ic,2,1) = element[ic].node[2]; tmp(ic,2,2) = element[ic].node[3];
    tmp(ic,3,0) = element[ic].node[0]; tmp(ic,3,1) = element[ic].node[3]; tmp(ic,3,2) = element[ic].node[2];
  };

  auto detectCorrespondingElement=[](VECTOR2D<int> ieb,const int node1,const int node2,const int node3,const int ic)->int{
    for(int ic=0;ic<ieb[node1].size();ic++){

      int element1 = ieb[node1][ic];
      if(element1==ic) continue;

      for(int ic2=0;ic2<ieb[node2].size();ic2++){

        int element2 = ieb[node2][ic2];
        if(element1!=element2) continue;

        for(int ic3=0;ic3<ieb[node3].size();ic3++){
          if(element1==ieb[node3][ic3]) return element1;
        }
      }
    }
    return -1;
  };

  auto searchCorrespondingFaceNumber=[](ARRAY3D<int> &tmp,const int node1,const int node2,const int node3,const int ic)->int{
    for(int i=0;i<4;i++){
      if(node1!=tmp(ic,i,0) && node1!=tmp(ic,i,1) && node1!=tmp(ic,i,2)) continue;
      if(node2!=tmp(ic,i,0) && node2!=tmp(ic,i,1) && node2!=tmp(ic,i,2)) continue;
      if(node3!=tmp(ic,i,0) && node3!=tmp(ic,i,1) && node3!=tmp(ic,i,2)) continue;
      return i;
    }
    return -1;
  };

  ARRAY3D<int> face_tmp(numOfElm,4,3);
  ARRAY2D<bool> selected(numOfElm,4);
    
  //initialize
  for(int ic=0;ic<numOfElm;ic++){
    element[ic].face.resize(4);
    facesInElement(face_tmp,element,ic);
    for(int j=0;j<4;j++) selected(ic,j) = false;
  }

  for(int ic=0;ic<numOfElm;ic++){
    for(int i=0;i<4;i++){
      if(selected(ic,i)==true) continue;

      int detectedElement = detectCorrespondingElement(ieb,face_tmp(ic,i,0),face_tmp(ic,i,1),face_tmp(ic,i,2),ic);

      FaceType face_dummy;
      face_dummy.node.resize(3);
      for(int j=0;j<3;j++) face_dummy.node[j]=face_tmp(ic,i,j);
      face_dummy.element[0] = ic;
      face_dummy.element[1] = detectedElement;
      face.push_back(face_dummy);
      selected(ic,i) = true;

      element[ic].face[i] = face.size()-1;

      if(detectedElement==-1) continue;
      int faceNumber = searchCorrespondingFaceNumber(face_tmp,face_tmp(ic,i,0),face_tmp(ic,i,1),face_tmp(ic,i,2),detectedElement);
      if(faceNumber==-1){
        cout <<"Error: wrong face detected. Exit..." << endl;
        exit(1);
      }
      selected(detectedElement,faceNumber) = true;
      element[detectedElement].face[faceNumber] = face.size()-1;
    }
  }

}

// #################################################################
/**
 * @brief calc boundary conditions
 * @param [in] stress
 */
void SmoothedFEM::SFEM::export_vtu_face(const string &file)
{
  FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL) {
    cout << file << " open error" << endl;
    exit(1); 
  }

  fprintf(fp,"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<UnstructuredGrid>\n");
  fprintf(fp,"<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n",numOfNode,(int)face.size());
  fprintf(fp,"<Points>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",x(i,0),x(i,1),x(i,2));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Points>\n");
  fprintf(fp,"<Cells>\n");
  fprintf(fp,"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for(int i=0;i<face.size();i++){
    for(int j=0;j<face[i].node.size();j++) fprintf(fp,"%d ",face[i].node[j]);
    fprintf(fp,"\n");
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num=0;
  for(int i=0;i<face.size();i++){
    num += face[i].node.size();
    fprintf(fp,"%d\n",num);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for(int i=0;i<face.size();i++) fprintf(fp,"%d\n",face[i].meshType);
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Cells>\n");

  fprintf(fp,"<PointData>\n");
  // fprintf(fp,"<DataArray type=\"Float64\" Name=\"displacement[m/s]\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  // for(int i=0;i<numOfNode;i++){
  //   fprintf(fp,"%e %e %e\n",U(i,0),U(i,1),U(i,2));
  // }
  // fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");

  fprintf(fp,"<CellData>");
  // fprintf(fp,"<DataArray type=\"Int64\" Name=\"Material\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  // for(int i=0;i<numOfElm;i++){
  //   fprintf(fp,"%d\n",element[i].materialType);
  // }
  // fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");
  fprintf(fp,"</Piece>");
  fprintf(fp,"</UnstructuredGrid>");
  fprintf(fp,"</VTKFile>");
  fclose(fp);
}
