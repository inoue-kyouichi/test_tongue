/**
 * @file fileIO.cpp
 * @brief meshInfo class
 * @brief solverInfo class
 * @brief outputInfo class
 * @author T. Otani
 */

#include "fileIO.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>

using namespace std;

// #################################################################
/**
 * @brief counting number of column
 */
int fileIO::CountNumbersOfTextLines(const string &filePath )
{
  long i = 0;

  ifstream ifs( filePath );

  if( ifs )
  {
    string line;

    while( true )
    {
      getline( ifs, line );
      i++;
      if( ifs.eof() )
        break;
    }
  }
  return i-1;
}

// #################################################################
/**
 * @brief import triangle node
 */
void fileIO::read_geometry_node(DOUBLEARRAY2 &x,int &numOfNode,const string &file)
{
  FILE *fp;
  string str,tmp;
  numOfNode = CountNumbersOfTextLines(file);

  x=Allocation::allocate2dDOUBLE(numOfNode,3);

  if ((fp = fopen(file.c_str(), "r")) == NULL) {
    cout << "file open error" << endl;
    exit(1); 
  }
  for(int i=0;i<numOfNode;i++){
    fscanf(fp,"%lf %lf %lf\n",&x[i][0],&x[i][1],&x[i][2]);
  }
  fclose(fp);

  // ifstream node_file(file);
  // if(!node_file){
  //   cout << "Error:Input "<< file << " not found" << endl;
  //   exit(1);
  // }

  // for(int i=0;i<numOfNode;i++){

  //   getline(node_file,str);
  //   istringstream stream(str);

  //   for(int j=0;j<3;j++){
  //     getline(stream,tmp,' ');
  //     x[i][j] = stof(tmp);
  //   }
  // }
}
// #################################################################
/**
 * @brief import triangle element
 */
void fileIO::read_geometry_meshType(elementType &element,int &numOfElm,const string &file)
{
  string str,tmp;
  numOfElm = CountNumbersOfTextLines(file);
  element.resize(numOfElm);

  ifstream meshType_file(file);
  if(!meshType_file){
    cout << "Error:Input "<< file << " not found" << endl;
    exit(1);
  }
  for(int i=0;i<numOfElm;i++){
    getline(meshType_file,str);
    element[i].meshType = VTKCellType(stoi(str));
    switch(element[i].meshType){
      case VTK_TRIANGLE:
        element[i].node.resize(3);
        break;
      case VTK_QUAD:
        element[i].node.resize(4);
        break;
      case VTK_TETRA:
        element[i].node.resize(4);
        break;
      case VTK_HEXAHEDRON:
        element[i].node.resize(8);
        break;
      case VTK_WEDGE:
        element[i].node.resize(6);
        break;
      default:
        cout << "invalid mesh type" << endl;
        exit(1);
      }
  }
}

// #################################################################
/**
 * @brief import triangle element
 */
void fileIO::read_geometry_element(elementType &element,const int &numOfElm,const string &file)
{
  string str,tmp;
  ifstream element_file(file);
  if(!element_file){
    cout << "Error:Input "<< file << " not found" << endl;
    exit(1);
  }

  for(int i=0;i<numOfElm;i++){

    getline(element_file,str);
    istringstream stream(str);
    for(int j=0,n=element[i].node.size();j<n;j++){
      getline(stream,tmp,' ');
      element[i].node[j] = stoi(tmp);
    }
  }
}

// #################################################################
/**
 * @brief calc boundary conditions
 * @param [in] stress
 */
void fileIO::export_vtu(const DOUBLEARRAY2 &x,const elementType &element,
            const int &numOfNode,const int &numOfElm,
            const DOUBLEARRAY2 &U,const DOUBLEARRAY1 &volumeChangeRatio,const DOUBLEARRAY2 &lambda_ave,
            const DOUBLEARRAY2 &sigmaEigen_ave,const DOUBLEARRAY2 &AEigen_ave,
            const DOUBLEARRAY3 &sigmaEigenVector_ave,const DOUBLEARRAY3 &AEigenVector_ave,const DOUBLEARRAY2 &innerForce,const string &file)
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
    fprintf(fp,"%e %e %e\n",x[i][0],x[i][1],x[i][2]);
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
    fprintf(fp,"%e %e %e\n",U[i][0],U[i][1],U[i][2]);
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float64\" Name=\"innerForce[N]\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",innerForce[i][0],innerForce[i][1],innerForce[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");

  fprintf(fp,"<CellData>");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"volumeChangeRatio\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e\n",volumeChangeRatio[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"averageLambda\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e %e %e\n",lambda_ave[i][0],lambda_ave[i][1],lambda_ave[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"PrincipalStress\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e %e %e\n",sigmaEigen_ave[i][0],sigmaEigen_ave[i][1],sigmaEigen_ave[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"firstPrincipalStress\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e %e %e\n",sigmaEigenVector_ave[i][0][0],sigmaEigenVector_ave[i][0][1],sigmaEigenVector_ave[i][0][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"PrincipalStrain\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e %e %e\n",AEigen_ave[i][0],AEigen_ave[i][1],AEigen_ave[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"firstPrincipalStrain\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e %e %e\n",AEigenVector_ave[i][0][0],AEigenVector_ave[i][0][1],AEigenVector_ave[i][0][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");
  fprintf(fp,"</Piece>");
  fprintf(fp,"</UnstructuredGrid>");
  fprintf(fp,"</VTKFile>");
  fclose(fp);
}

// #################################################################
/**
 * @brief calc boundary conditions
 * @param [in] stress
 */
void fileIO::export_vtu(const DOUBLEARRAY2 &x,const elementType &element,
            const int &numOfNode,const int &numOfElm,
            const DOUBLEARRAY2 &U,const DOUBLEARRAY1 &volumeChangeRatio,const DOUBLEARRAY2 &lambda_ave,
            const DOUBLEARRAY1 &bundle,const INTARRAY1 &bundleElement,
            const DOUBLEARRAY2 &sigmaEigen_ave,const DOUBLEARRAY2 &AEigen_ave,
            const DOUBLEARRAY3 &sigmaEigenVector_ave,const DOUBLEARRAY3 &AEigenVector_ave,const DOUBLEARRAY2 &innerForce,const string &file)
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
    fprintf(fp,"%e %e %e\n",x[i][0],x[i][1],x[i][2]);
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
    fprintf(fp,"%e %e %e\n",U[i][0],U[i][1],U[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  // fprintf(fp,"<DataArray type=\"Float64\" Name=\"bundle\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  // for(int i=0;i<numOfNode;i++){
  //   fprintf(fp,"%e\n",bundle[i]);
  // }
  // fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float64\" Name=\"innerForce[N]\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",-innerForce[i][0],-innerForce[i][1],-innerForce[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");

  fprintf(fp,"<CellData>");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"volumeChangeRatio\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e\n",volumeChangeRatio[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"averageLambda\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e %e %e\n",lambda_ave[i][0],lambda_ave[i][1],lambda_ave[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"PrincipalStress\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e %e %e\n",sigmaEigen_ave[i][0],sigmaEigen_ave[i][1],sigmaEigen_ave[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"firstPrincipalStress\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e %e %e\n",sigmaEigenVector_ave[i][0][0],sigmaEigenVector_ave[i][0][1],sigmaEigenVector_ave[i][0][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"PrincipalStrain\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e %e %e\n",AEigen_ave[i][0],AEigen_ave[i][1],AEigen_ave[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"firstPrincipalStrain\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e %e %e\n",AEigenVector_ave[i][0][0],AEigenVector_ave[i][0][1],AEigenVector_ave[i][0][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Int64\" Name=\"bundle\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%d\n",bundleElement[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");
  fprintf(fp,"</Piece>");
  fprintf(fp,"</UnstructuredGrid>");
  fprintf(fp,"</VTKFile>");
  fclose(fp);
}


// #################################################################
/**
 * @brief calc boundary conditions
 * @param [in] stress
 */
void fileIO::export_vtu(const DOUBLEARRAY2 &x,const elementType &element,
            const int &numOfNode,const int &numOfElm,
            const DOUBLEARRAY2 &U,const DOUBLEARRAY1 &volumeChangeRatio,const DOUBLEARRAY2 &lambda_ave,const string &file)
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
    fprintf(fp,"%e %e %e\n",x[i][0],x[i][1],x[i][2]);
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
    fprintf(fp,"%e %e %e\n",U[i][0],U[i][1],U[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");

  fprintf(fp,"<CellData>");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"volumeChangeRatio\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e\n",volumeChangeRatio[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"averageLambda\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e %e %e\n",lambda_ave[i][0],lambda_ave[i][1],lambda_ave[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");
  fprintf(fp,"</Piece>");
  fprintf(fp,"</UnstructuredGrid>");
  fprintf(fp,"</VTKFile>");
  fclose(fp);
}

// #################################################################
/**
 * @brief calc boundary conditions
 * @param [in] stress
 */
void fileIO::export_vtu(const DOUBLEARRAY2 &x,const elementType &element,
            const int &numOfNode,const int &numOfElm,
            const DOUBLEARRAY2 &U,const DOUBLEARRAY1 &volumeChangeRatio,const string &file)
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
    fprintf(fp,"%e %e %e\n",x[i][0],x[i][1],x[i][2]);
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
    fprintf(fp,"%e %e %e\n",U[i][0],U[i][1],U[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");

  fprintf(fp,"<CellData Scalars=\"volumeChangeRatio\">\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"volumeChangeRatio\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e\n",volumeChangeRatio[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");
  fprintf(fp,"</Piece>");
  fprintf(fp,"</UnstructuredGrid>");
  fprintf(fp,"</VTKFile>");
  fclose(fp);
}

// #################################################################
/**
 * @brief calc boundary conditions
 * @param [in] stress
 */
void fileIO::export_vtu_boundary(const DOUBLEARRAY2 &x,const elementType &element,
            const int &numOfNode,const int &numOfElm,
            const INTARRAY2 &ibd,const DOUBLEARRAY2 &bd,const DOUBLEARRAY2 &fiberDirection_elm,const string &file)
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
    fprintf(fp,"%e %e %e\n",x[i][0],x[i][1],x[i][2]);
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

  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"bd\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",bd[i][0],bd[i][1],bd[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Int64\" Name=\"ibd\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%d %d %d\n",ibd[i][0],ibd[i][1],ibd[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"<CellData>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"fiberDirection\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e %e %e\n",fiberDirection_elm[i][0],fiberDirection_elm[i][1],fiberDirection_elm[i][2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");

  fprintf(fp,"</Piece>");
  fprintf(fp,"</UnstructuredGrid>");
  fprintf(fp,"</VTKFile>");
  fclose(fp);
}

