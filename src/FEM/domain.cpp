/**
 * @file domain.cpp
 * @brief domain class
 * @author T. Otani
 */

#include "domain.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include "fileIO.h"

using namespace std;

// #################################################################
/**
 * @brief read boundary information
 * @param [in] D_file filePath of Dirichlet condition
 */
void Domain::set_dirichlet(const string &D_file)
{
  int node;
  string str,tmp,dtmp;

  ifstream file(D_file);
  if(!file){
    cout << "Error:Input "<< D_file << " not found" << endl;
    exit(1);
  }

  ibd.allocate(numOfNode,3);
  bd.allocate(numOfNode,3);

  for(int i=0;i<numOfNode;i++){
    getline(file,str);
    istringstream stream(str);
    for(int j=0;j<3;j++){
      getline(stream,tmp,' ');
      ibd(i,j) = stoi(tmp);
    }
  }
}

// #################################################################
/**
 * @brief read boundary information
 * @param [in] D_file filePath of Dirichlet condition
 */
void Domain::set_dirichlet(const string &D_file,const string &Dvalue_file)
{
  int node;
  string str,tmp,dtmp;

  ifstream file(D_file);
  if(!file){
    cout << "Error:Input "<< D_file << " not found" << endl;
    exit(1);
  }

  ibd.allocate(numOfNode,3);
  bd.allocate(numOfNode,3);

  for(int i=0;i<numOfNode;i++){
    getline(file,str);
    istringstream stream(str);
    for(int j=0;j<3;j++){
      getline(stream,tmp,' ');
      ibd(i,j) = stoi(tmp);
    }
  }

  ifstream file2(Dvalue_file);
  if(!file){
    cout << "Error:Input "<< Dvalue_file << " not found" << endl;
    exit(1);
  }

  for(int i=0;i<numOfNode;i++){
    getline(file2,str);
    istringstream stream(str);
    for(int j=0;j<3;j++){
      getline(stream,tmp,' ');
      bd(i,j) = stof(tmp);
    }
  }

}
// #################################################################
/**
 * @brief read boundary information
 * @param [in] N_file filePath of Neumann condition
 */
void Domain::set_neumann(const string &N_file,const string &meshTypeFile,const string &NvalueFile)
{
  int node;
  string str,tmp;
  fileIO::read_geometry_meshType(belement,numOfNeumann,meshTypeFile);

  ifstream file2(N_file);
  if(!file2){
    cout << "Error:Input "<< N_file << " not found" << endl;
    exit(1);
  }
  for(int i=0;i<numOfNeumann;i++){
    getline(file2,str);
    istringstream stream(str);

    for(int j=0;j<belement[i].node.size();j++){
      getline(stream,tmp,' ');
      belement[i].node[j] = stoi(tmp);
    }
  }

  ifstream file3(NvalueFile);
  if(!file3){
    cout << "Error:Input "<< NvalueFile << " not found" << endl;
    exit(1);
  }
  bn.allocate(numOfNeumann,3);
  for(int i=0;i<numOfNeumann;i++){
    getline(file3,str);
    istringstream stream(str);
    for(int j=0;j<3;j++){
      getline(stream,tmp,' ');
      bn(i,j) = stof(tmp);
    }
  }

}
// #################################################################
/**
 * @brief read boundary information
 * @param [in] file1 filePath of nodes
 * @param [in] file2 filePath of elements
 */
void Domain::set_geometry(const string &file1,const string &file2,const string &file3)
{
  fileIO::read_geometry_node(x,numOfNode,file1);
  fileIO::read_geometry_meshType(element,numOfElm,file3);
  fileIO::read_geometry_element(element,numOfElm,file2);

  ieb.resize(numOfNode,vector<int>(0));
  inb.resize(numOfNode,vector<int>(0));

  x0.allocate(numOfNode,3);

  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) x0(i,j) = x(i,j);
  }

  calc_adjacent_elements();
  calc_adjacent_nodes();

}

// #################################################################
/**
 * @brief calc boundary conditions
 * @param [in] stress
 */
void Domain::calc_adjacent_nodes()
{
  int tmp;
  bool b1;

  for(int ic=0;ic<numOfNode;ic++){
    for(int i=0,n=ieb[ic].size();i<n;i++){
      for(int j=0;j<element[ieb[ic][i]].node.size();j++){
        tmp = element[ieb[ic][i]].node[j];
        b1 = true;
        for(int k=0,n=inb[ic].size();k<n;k++){
          if(inb[ic][k]==tmp){
            b1 = false;
            break;
          }
        }
        if(b1==false) continue;
        inb[ic].push_back(tmp);
      }
    }
    sort(inb[ic].begin(),inb[ic].end());
  }
}

// #################################################################
/**
 * @brief calc boundary conditions
 * @param [in] stress
 */
void Domain::calc_adjacent_elements()
{
  int tmp;
  for(int ic=0;ic<numOfElm;ic++){
    for(int j=0,n=element[ic].node.size();j<n;j++){
      tmp = element[ic].node[j];
      ieb[tmp].push_back(ic);
    }
  }
}

// #################################################################
/**
 * @brief calc boundary conditions
 * @param [in] stress
 */
void Domain::export_vtk(const string &file)
{
  ofstream ofs(file);

        ofs << "# vtk DataFile Version 2.0" << endl;
        ofs << "test" << endl;
        ofs << "ASCII" << endl;
        ofs << "DATASET UNSTRUCTURED_GRID" << endl;
        ofs << "POINTS "<< numOfNode << " double" << endl;

        for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) ofs << scientific << x(i,j) << " ";
    ofs << endl;
  }

  int totalNode=0;
  for(int i=0;i<numOfElm;i++){
    totalNode+=element[i].node.size()+1;
  }

  ofs << "CELLS " << numOfElm << " " << totalNode << endl;
  for(int i=0;i<numOfElm;i++){
    ofs << element[i].node.size() << " ";
    for(int j=0;j<element[i].node.size();j++) ofs << element[i].node[j] << " ";
    ofs << endl;
  }

  ofs << "CELL_TYPES " <<numOfElm << endl;
  for(int i=0;i<numOfElm;i++) ofs << element[i].meshType << endl;
}
