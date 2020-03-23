//##################################################################################
//
// test: calc surface curvature
//
// Copyright (c) 2020  Mechanical and Bioengineering Systems Lab.,
//                     Department of Mechanical Science and Bioengineering,
//                     Graduate School of Engineering Science,
//                     Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   main_curvature.cpp
 * @author T. Otani
 * @detail Fem solid analysis with displacement control.
 */

#include "triangleSurfaceCurvature.h"
#include "glog/logging.h"

void read_plyFile(const std::string &file,DOUBLEARRAY2D &x,INTARRAY2D &ie,int &numOfNode,int &numOfElm);
void export_vtu(DOUBLEARRAY2D &x,INTARRAY2D &ie,const int &numOfNode,const int &numOfElm,DOUBLEARRAY2D &normal,DOUBLEARRAY1D &meanCurvature,const std::string &file);

using namespace std;

int main(int argc,char *argv[])
{
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();

  int numOfNode,numOfElm;
  DOUBLEARRAY2D x;
  INTARRAY2D ie;
  DOUBLEARRAY2D normal;
  DOUBLEARRAY1D meanCurvature;

  //read triangle surface
  if(argc!=2){
    cout << "please set input file name. Exit..." << endl;
    exit(1);
  }
  std::string input_file = argv[1];
  read_plyFile(input_file,x,ie,numOfNode,numOfElm);
  triangleSurfaceCurvature surface(x,ie,numOfNode,numOfElm);

  //calc triangle surface curvature
  surface.calcSurfaceMeanCurvature();

  //output mean curvature
  surface.exportSurfaceNormal(normal);
  surface.exportSurfaceMeanCurvature(meanCurvature);

  //export triangle surface data with curvature
  std::string output_file = "test.vtu";
  export_vtu(x,ie,numOfNode,numOfElm,normal,meanCurvature,output_file);
  return 0;
}

// #################################################################
/**
 * @brief read ply file
 */
void read_plyFile(const std::string &file,DOUBLEARRAY2D &x,INTARRAY2D &ie,int &numOfNode,int &numOfElm)
{
  int i,j,ret,dummy;
  char tmp[100],tmp2[100];
  FILE *fp;

  if((fp=fopen(file.c_str(),"r"))==NULL){
    cout << file << " open error" << endl;
    exit(1);
  }
  fgets(tmp,100,fp);
  fgets(tmp,100,fp);
  fgets(tmp,100,fp);
  fgets(tmp,100,fp);

  fscanf(fp,"%s %s %d\n",tmp,tmp2,&numOfNode);

  fgets(tmp,100,fp);
  fgets(tmp,100,fp);
  fgets(tmp,100,fp);

  fscanf(fp,"%s %s %d\n",tmp,tmp2,&numOfElm);

  fgets(tmp,100,fp);
  fgets(tmp,100,fp);

  x.allocate(numOfNode,3);
  ie.allocate(numOfElm,3);

  for(i=0;i<numOfNode;i++){
    fscanf(fp,"%lf %lf %lf\n",&x(i,0), &x(i,1), &x(i,2));
  }
  for(i=0;i<numOfElm;i++){
    fscanf(fp,"%d %d %d %d\n",&dummy,&ie(i,0), &ie(i,1), &ie(i,2));
  }

  fclose(fp);
}

// #################################################################
/**
 * @brief export surfaces with normal & curvature
 */
void export_vtu(DOUBLEARRAY2D &x,INTARRAY2D &ie,const int &numOfNode,const int &numOfElm,DOUBLEARRAY2D &normal,DOUBLEARRAY1D &meanCurvature,const std::string &file)
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
    fprintf(fp,"%d %d %d",ie(i,0),ie(i,1),ie(i,2));
    fprintf(fp,"\n");
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num=0;
  for(int i=0;i<numOfElm;i++){
    num += 3;
    fprintf(fp,"%d\n",num);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++) fprintf(fp,"5\n");
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Cells>\n");

  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"normal\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",normal(i,0),normal(i,1),normal(i,2));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"meanCurvature\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e\n",meanCurvature(i));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"<CellData>\n");
  fprintf(fp,"</CellData>\n");

  fprintf(fp,"</Piece>");
  fprintf(fp,"</UnstructuredGrid>");
  fprintf(fp,"</VTKFile>");
  fclose(fp);
}

