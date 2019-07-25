/**
 * @file TriangleSet.cpp
 * @brief rigid body class
 * @author T. Otani
 */
#include "rigidBody.h"

using namespace std;

// #########################################################
/**
 * @brief read ply file (ascii)
 */
void TriangleSet::translation(const double (&center)[3])
{
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) x[i][j]-=center[j];
    for(int j=0;j<3;j++) x0[i][j]=x[i][j];
  }
}

// #########################################################
/**
 * @brief read ply file (ascii)
 */
void TriangleSet::readPLY(const std::string &file)
{
  int ret,dummy;
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

  x=Allocation::allocate2dDOUBLE(numOfNode,3);
  x0=Allocation::allocate2dDOUBLE(numOfNode,3);
  elm=Allocation::allocate2dINT(numOfElm,3);

  for(int i=0;i<numOfNode;i++){
    fscanf(fp,"%lf %lf %lf\n",&x0[i][0], &x0[i][1], &x0[i][2]);
  }
  for(int i=0;i<numOfElm;i++){
    fscanf(fp,"%d %d %d %d\n",&dummy,&elm[i][0], &elm[i][1], &elm[i][2]);
  }
  fclose(fp);

  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) x[i][j]=x0[i][j];
  }

}