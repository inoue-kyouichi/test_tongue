/**
 * @file fem_preprocessing.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "fem.h"

using namespace std;

// #################################################################
/**
 * @brief initialize FEM class
 */
void Fem::initialize()
{
  inputDomainInfo();
  inputSolverInfo();
  inputOutputInfo();

  string restartDirName="Restart_"+to_string(dataNumber);
  mkdir(restartDirName.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);

  allocate();
  restart_setting();

  inputDirichletBoundaryInfo();
  // setFiberDirection_KogaModel();
  // inputFiberInfo();
}

// #################################################################
/**
 * @brief output information from TP file
 */
void Fem::restart_setting()
{
  FILE *fp;
  string output;

  string str,base_label,label;

  base_label = "/Domain";
  if(Restart==0){
    output = "Restart_"+to_string(dataNumber)+"/x0.dat";
    if ((fp = fopen(output.c_str(), "w")) == NULL) {
      cout << output << " open error" << endl;
      exit(1);
    }
    for(int i=0;i<numOfNode;i++){
      fprintf(fp,"%e %e %e\n",x0(i,0),x0(i,1),x0(i,2));
    }
    fclose(fp);
  }else{
    label = base_label + "/referenceNodeFile";
    if ( !tp.getInspectedValue(label, output)){
      cout << output <<" is not set" << endl;
      exit(0);
    }
    if ((fp = fopen(output.c_str(), "r")) == NULL) {
      cout << output << " open error" << endl;
      exit(1);
    }
    for(int i=0;i<numOfNode;i++){
      fscanf(fp,"%lf %lf %lf\n",&x0(i,0),&x0(i,1),&x0(i,2));
      for(int j=0;j<3;j++) x(i,j)=x0(i,j);
    }
    fclose(fp);

    label = base_label + "/restartNodeFile";
    if ( !tp.getInspectedValue(label, output)){
      cout << "outputDir is not set" << endl;
      exit(0);
    }
    if ((fp = fopen(output.c_str(), "r")) == NULL) {
      cout << output << " open error" << endl;
      exit(1);
    }
    for(int i=0;i<numOfNode;i++){
      fscanf(fp,"%lf %lf %lf\n",&U(i,0),&U(i,1),&U(i,2));
    }
    fclose(fp);
  }
}

// #################################################################
/**
 * @brief output information from TP file
 */
void Fem::inputOutputInfo()
{
  string str,base_label,label;

  base_label = "/Output";
  label = base_label + "/outputDir";
  if ( !tp.getInspectedValue(label, outputDir)){
    cout << "outputDir is not set" << endl;
    exit(0);
  }

  label = base_label + "/fileName";
  if ( !tp.getInspectedValue(label, fileName)){
    cout << "fileName is not set" << endl;
    exit(0);
  }
}

// #################################################################
/**
 * @brief solver information from TP file
 */
void Fem::inputSolverInfo()
{
  string str,base_label,label;
  int tmp;

  base_label = "/Solver";

  label = base_label + "/dataNumber";
  if ( !tp.getInspectedValue(label, dataNumber)){
    cout << "maxiteration is not set" << endl;
    exit(0);
  }

  label = base_label + "/maxIteration";
  if ( !tp.getInspectedValue(label, maxIteration)){
    cout << "maxiteration is not set" << endl;
    exit(0);
  }

  label = base_label + "/NR_iteration";
  if ( !tp.getInspectedValue(label, NRiteration)){
    cout << "NLiteration is not set" << endl;
    exit(0);
  }

  label = base_label + "/NR_tolerance";
  if ( !tp.getInspectedValue(label, tmp)){
    cout << "NRtolerance is not set" << endl;
    exit(0);
  }
  NRtolerance = pow(1e-1,tmp);

  label = base_label + "/Restart";
  if ( !tp.getInspectedValue(label, Restart)){
    cout << label <<" is not set" << endl;
    exit(0);
  }

  label = base_label + "/OMPnumThreads";
  if ( !tp.getInspectedValue(label, OMPnumThreads)){
    cout << "OMPnumThreads is not set" << endl;
    exit(0);
  }

  label = base_label + "/relaxation";
  if ( !tp.getInspectedValue(label,relaxation)){
    cout << "NLiteration is not set" << endl;
    exit(0);
  }
}



// #################################################################
/**
 * @brief domain information from tp file
 */
void Fem::inputDomainInfo()
{
  string str,base_label,label,inputDir;

  base_label = "/Domain";

  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }
  string file1,file2,file3,file4,file5;
  label = base_label + "/nodeFile";
  if ( !tp.getInspectedValue(label, file1)){
    cout << label << " is not found" << endl;
    exit(0);
  }
  label = base_label + "/elementFile";
  if ( !tp.getInspectedValue(label,file2)){
    cout << label << " is not found" << endl;
    exit(0);
  }
  label = base_label + "/meshTypeFile";
  if ( !tp.getInspectedValue(label,file3)){
    cout << label << " is not found" << endl;
    exit(0);
  }

  label = base_label + "/materialTypeFile";
  if ( !tp.getInspectedValue(label,file4)){
    cout << label << " is not found" << endl;
    exit(0);
  }
  file1=inputDir+"/"+file1;
  file2=inputDir+"/"+file2;
  file3=inputDir+"/"+file3;
  set_geometry(file1,file2,file3);

  fileIO::read_geometry_materialType(element,numOfElm,file4);
}

// #################################################################
/**
 * @brief Dirichlet information from TP file
 */
void Fem::inputDirichletBoundaryInfo()
{
  string str,base_label,label,inputDir;
  base_label = "/Domain";

  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }
  string D_file,Dvalue_file;
  label = base_label + "/dirichletFile";
  if ( !tp.getInspectedValue(label, D_file)){
    cout << label << " is not found" << endl;
  }
  label = base_label + "/dirichletValueFile";
  if ( !tp.getInspectedValue(label, Dvalue_file)){
    cout << label << " is not found" << endl;
  }
  D_file=inputDir+"/"+D_file;
  Dvalue_file=inputDir+"/"+Dvalue_file;
  set_dirichlet(D_file,Dvalue_file);
}

// #################################################################
/**
 * @brief Neumann information from TP file
 */
void Fem::inputNeumannBoundaryInfo()
{
  string str,base_label,label,inputDir;
  string N_file,Nvalue_file,NmeshType_file;

  base_label = "/Domain";
  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }

  label = base_label + "/neumannFile";
  if ( !tp.getInspectedValue(label, N_file)){
    cout << label << " is not found" << endl;
  }
  label = base_label + "/neumannValueFile";
  if ( !tp.getInspectedValue(label, Nvalue_file)){
    cout << label << " is not found" << endl;
  }
  label = base_label + "/neumannMeshTypeFile";
  if ( !tp.getInspectedValue(label, NmeshType_file)){
    cout << label << " is not found" << endl;
  }
  set_neumann(N_file,NmeshType_file,Nvalue_file);
}

// #################################################################
/**
 * @brief fiber information from TP file
 * @TODO linear hexa element only
 */
void Fem::setFiberDirection()
{
  double X_elm[4][3],normal[3];
  for(int ic=0;ic<numOfElm;ic++){
    for(int p=0;p<4;p++){
      for(int j=0;j<3;j++) X_elm[p][j]=x(element[ic].node[p],j);
    }
    calc_normal_quad(normal,X_elm);

    for(int j=0;j<3;j++) fiberDirection_elm(ic,j)=normal[j];
  }
}

// #################################################################
/**
 * @brief fiber information from TP file
 * @TODO linear hexa element only
 */
void Fem::calc_normal_quad(double (&normal)[3],double (&X)[4][3])
{
  double dXdr1[3],dXdr2[3],dXdr3[3],Jacobian,normalVector[3];
  double tmp;
  DOUBLEARRAY2D dNdr(4,2);
  DOUBLEARRAY2D dNdX(4,3);

  Gauss g1(2);

  for(int j=0;j<3;j++) normal[j] = 0e0;

  for(int iy=0;iy<2;iy++){
    for(int ix=0;ix<2;ix++){

      ShapeFunction2D::C2D4_dNdr(dNdr,g1.point[ix],g1.point[iy]);
      for(int i=0;i<3;i++){
        dXdr1[i] = 0e0;
        dXdr2[i] = 0e0;
        for(int p=0;p<4;p++){
          dXdr1[i] += dNdr(p,0) * X[p][i];
          dXdr2[i] += dNdr(p,1) * X[p][i];
        }
      }

      mathTool::crossProduct(dXdr1,dXdr2,dXdr3,Jacobian);

      for(int i=0;i<3;i++) normalVector[i] = dXdr3[i] / Jacobian;
      for(int i=0;i<3;i++) normal[i] += normalVector[i];
    }
  }

  tmp = sqrt(pow(normal[0],2.0e0)+pow(normal[1],2.0e0)+pow(normal[2],2.0e0));
  for(int i=0;i<3;i++) normal[i] /= tmp;
}

// #################################################################
/**
 * @brief fiber information from TP file
 */
void Fem::setFiberDirection_KogaModel()
{
  string str,base_label,label,inputDir;
  double thetaMax,horizontalFiberPosition;

  base_label = "/Domain";

  label = base_label + "/thetaMax";
  if ( !tp.getInspectedValue(label,thetaMax)){
    cout << "data format is not set" << endl;
    exit(0);
  }
  label = base_label + "/horizontalFiberPosition";
  if ( !tp.getInspectedValue(label,horizontalFiberPosition)){
    cout << "data format is not set" << endl;
    exit(0);
  }

  double A1=4e0*sqrt(3)*thetaMax/9e0;
  double A2=2e0*sqrt(3)*thetaMax/9e0;
  double omega=1e0/horizontalFiberPosition;
  double s_coordinates,zAxis[3]={0e0,0e0,1e0};
  double angle;

  double X_elm[4][3],normal[3],center[3];
  double R[3][3],rotAxis[3],absoluteValue,direction[3];

  double length,s_max = arcLength(0e0,3e-3);

  for(int ic=0;ic<numOfElm;ic++){

    for(int p=0;p<4;p++){
      for(int j=0;j<3;j++) X_elm[p][j]=x(element[ic].node[p],j);
    }

    calc_normal_quad(normal,X_elm);
    for(int j=0;j<3;j++) normal[j]=-1e0*normal[j];

    mathTool::crossProduct(normal,zAxis,rotAxis,absoluteValue);
    for(int j=0;j<3;j++) rotAxis[j]/=absoluteValue;
    for(int j=0;j<3;j++) center[j] = (X_elm[0][j]+X_elm[1][j]+X_elm[2][j]+X_elm[3][j])/4e0;

    length = sqrt(center[0]*center[0]+center[1]*center[1]);
    s_coordinates = PI*arcLength(0e0,length)/s_max;

    angle = A1*sin(omega*s_coordinates) + A2*sin(2e0*omega*s_coordinates);
    angle = angle * PI/180e0;

    mathTool::calcRotationMatrix(R,rotAxis,angle);

    for(int i=0;i<3;i++){
      direction[i]=0e0;
      for(int j=0;j<3;j++) direction[i] += R[i][j] * normal[j];
    }

    for(int j=0;j<3;j++) fiberDirection_elm(ic,j)=direction[j];
  }
}

// #################################################################
/**
 * @brief arcLength
 */
double Fem::arcLength(const double xMin,const double xMax)
{
  double s = (xMax*sqrt(4e6*xMax*xMax+1e0)/2e0+asinh(2e3*xMax)/4e3)-(xMin*sqrt(4e6*xMin*xMin+1e0)/2e0+asinh(2e3*xMin)/4e3);
  return s;
}

// #################################################################
/**
 * @brief fiber information from TP file
 */
void Fem::inputFiberInfo()
{
  string str,base_label,label,inputDir;

  base_label = "/Domain";

  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }

  string file1;
  fiberDirection.allocate(numOfElm,8,3); //gauss point
  label = base_label + "/fiberDirectionFile";
  if ( !tp.getInspectedValue(label,file1)){
    cout << label << " is not found" << endl;
    exit(0);
  }
  file1=inputDir+"/"+file1;
  FILE *fp;
  if ((fp = fopen(file1.c_str(), "r")) == NULL) {
    cout << file1 << " open error" << endl;
    exit(1);
  }
  for(int i=0;i<numOfElm;i++){
    for(int j=0;j<8;j++) fscanf(fp,"%lf %lf %lf\n",&fiberDirection(i,j,0),&fiberDirection(i,j,1),&fiberDirection(i,j,2));
  }
  fclose(fp);
}

// #################################################################
/**
 * @brief allocation. 高次要素を用いる場合は修正が必要。
 */
void Fem::allocate()
{
  Mass.allocate(numOfElm,20,20);
  U.allocate(numOfNode,3);
  RHS.allocate(numOfNode,3);
  innerForce.allocate(numOfNode,3);
  externalForce.allocate(numOfNode,3);
  K.allocate(numOfElm,20,20,3,3);
  Qu.allocate(numOfElm,20,3);
  BFe.allocate(numOfElm,9,3);
  volume.allocate(numOfElm);
  volume0.allocate(numOfElm);
  volumeChangeRatio.allocate(numOfElm);

  Mises.allocate(numOfElm);
  lambda_ave.allocate(numOfElm,3);
  AEigen_Ave.allocate(numOfElm,3);
  sigmaEigen_Ave.allocate(numOfElm,3);
  AEigenVector_Ave.allocate(numOfElm,3,3);
  sigmaEigenVector_Ave.allocate(numOfElm,3,3);

  fiberDirection_elm.allocate(numOfElm,3); //gauss point

  for(int ic=0;ic<numOfElm;ic++) volumeChangeRatio(ic)=1e0;

  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) U(i,j) = 0e0;
  }

  //dirichlet boundary conditions
  ibd.allocate(numOfNode,3);
  bd.allocate(numOfNode,3);
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++){
      ibd(i,j)=1;
      bd(i,j)=0e0;
    }
  }

}
