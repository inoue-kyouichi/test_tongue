/**
 * @file rigidBody.cpp
 * @brief rigid body class
 * @author T. Otani
 */
#include "rigidBody.h"

using namespace std;

void RigidBody::initialize(TextParser &tp)
{
  string label,fileName;
  label = "/RigidBody/shape";
  if ( !tp.getInspectedValue(label, fileName)){
    cout << fileName << "  is not set" << endl;
    exit(0);
  }
  obj.readPLY(fileName);
  rho=1e0;
  calcMassProperties();

  for(int j=0;j<3;j++) U[j]=0e0;

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      R[i][j]=0e0;
    }
    R[i][i]=1e0;
  }

}

// #########################################################
/**
 * @brief calc mass Properties by the method of Mirtich, Journal of Graphic Tools, 1, 31-50 (1996).
 */
void RigidBody::calcMassProperties()
{
  Mass=0e0;
  for(int i=0;i<3;i++){
    xg[i]=0e0;
    for(int j=0;j<3;j++) J[i][j]=0e0;
  }

  double a[3],b[3],n[3],center[3],Area,Tx2,Ty2,Tz2;

  for(int ic=0;ic<obj.numOfElm;ic++){
      for(int i=0;i<3;i++){
        center[i]=(obj.x(obj.elm(ic,0),i)+obj.x(obj.elm(ic,1),i)+obj.x(obj.elm(ic,2),i))/3e0;
        a[i]=obj.x(obj.elm(ic,1),i)-obj.x(obj.elm(ic,0),i);
        b[i]=obj.x(obj.elm(ic,2),i)-obj.x(obj.elm(ic,0),i);
      }
      mathTool::crossProduct(a,b,n,Area);
      for(int i=0;i<3;i++) n[i]/=Area;
      Area=5e-1*Area;
      Mass+=rho*(n[0]*center[0])*Area;
      xg[0]+=rho*5e-1*n[0]*center[0]*center[0]*Area;
      xg[1]+=rho*5e-1*n[1]*center[1]*center[1]*Area;
      xg[2]+=rho*5e-1*n[2]*center[2]*center[2]*Area;
      Tx2=rho*n[0]*center[0]*center[0]*center[0]*Area/3e0;
      Ty2=rho*n[1]*center[1]*center[1]*center[1]*Area/3e0;
      Tz2=rho*n[2]*center[2]*center[2]*center[2]*Area/3e0;
      J[0][0]+=Ty2+Tz2;
      J[1][1]+=Tz2+Tx2;
      J[2][2]+=Tx2+Ty2;
      J[0][1]+=rho*5e-1*n[0]*center[0]*center[0]*center[1]*Area;
      J[1][2]+=rho*5e-1*n[1]*center[1]*center[1]*center[2]*Area;
      J[2][0]+=rho*5e-1*n[2]*center[2]*center[2]*center[0]*Area;
  }
    for(int i=0;i<3;i++) xg[i]/=Mass;
    J[0][0]-=Mass*(xg[1]*xg[1]+xg[2]*xg[2]);
    J[1][1]-=Mass*(xg[2]*xg[2]+xg[0]*xg[0]);
    J[2][2]-=Mass*(xg[0]*xg[0]+xg[1]*xg[1]);
    J[0][1]-=Mass*(xg[0]*xg[1]);
    J[1][2]-=Mass*(xg[1]*xg[2]);
    J[2][0]-=Mass*(xg[2]*xg[0]);
    J[1][0]=J[0][1];
    J[2][1]=J[1][2];
    J[0][2]=J[2][0];

    mathTool::Jacobi3x3(10000,1e-8,J,J,coordinates_ref);
}

// #########################################################
/**
 * @brief export ply file (ascii format)
 */
void RigidBody::exportPLY(const std::string &file)
{
  FILE *fp;
  double R[3][3];
  // mathTool::quaternion2rotation(R,q);
  // for(int i=0;i<obj.numOfNode;i++){
  //   for(int j=0;j<3;j++){
  //     obj.x[i][j]=x[j];
  //     for(int k=0;k<3;k++) obj.x[i][j]+=R[j][k]*obj.x0[i][k];
  //   }
  // }

  if((fp=fopen(file.c_str(),"w"))==NULL){
    printf("file open error\n");
    exit(1);
  }
  fprintf(fp,"ply\n");
  fprintf(fp,"format ascii 1.0\n");
  fprintf(fp,"test\n");
  fprintf(fp,"obj_info testData\n");
  fprintf(fp,"element vertex %d\n",obj.numOfNode);
  fprintf(fp,"property float x\n");
  fprintf(fp,"property float y\n");
  fprintf(fp,"property float z\n");
  fprintf(fp,"element face %d\n",obj.numOfElm);
  fprintf(fp,"property list uchar int vertex_indices\n");
  fprintf(fp,"end_header\n");

  for(int i=0;i<obj.numOfNode;i++){
    fprintf(fp,"%e %e %e\n",obj.x(i,0),obj.x(i,1),obj.x(i,2));
  }
  for(int i=0;i<obj.numOfElm;i++){
    fprintf(fp,"3 %d %d %d\n",obj.elm(i,0),obj.elm(i,1),obj.elm(i,2));
  }
  fclose(fp);
}

// #########################################################
/**
 * @brief export ply file (ascii format)
 */
void RigidBody::updateRotationMatrix_spatialForm(const double (&w)[3])
{
  double R1[3][3],Rtmp[3][3],n[3],tilde[3][3];
  double theta=sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
  // printf("theta=%e\n",theta*180e0/PI);
  if(theta<1e-15) return;

  for(int j=0;j<3;j++) n[j]=w[j]/theta;
  mathTool::skewSymmetricTensor(tilde,n);

  double tilde2[3][3];
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      tilde2[i][j]=0e0;
      for(int k=0;k<3;k++){
        tilde2[i][j]+=tilde[i][k]*tilde[k][j];
      }
    }
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) R1[i][j]=sin(theta)*tilde[i][j]+(1e0-cos(theta))*tilde2[i][j];
    R1[i][i]+=1e0;
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      Rtmp[i][j]=0e0;
      for(int k=0;k<3;k++) Rtmp[i][j]+=R1[i][k]*R[k][j];
    }
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) R[i][j]=Rtmp[i][j];
  }
}

// #########################################################
/**
 * @brief export ply file (ascii format)
 */
void RigidBody::updateShape()
{
  double tmp[3],tmp2[3];
  for(int ic=0;ic<obj.numOfNode;ic++){
    for(int j=0;j<3;j++) tmp[j] = obj.x0(ic,j)-xg[j];
    for(int i=0;i<3;i++){
      tmp2[i]=0e0;
      for(int j=0;j<3;j++) tmp2[i]+=R[i][j]*tmp[j];
    }
   for(int j=0;j<3;j++) obj.x(ic,j)=xg[j]+tmp2[j]+U[j];
  }
}
