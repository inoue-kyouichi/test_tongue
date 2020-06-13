/**
 * @file fem_preprocessing.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "PDL_analysis.h"

using namespace std;

// #################################################################
/**
 * @brief calc stress tensor
 */
void humanPDL::PeriodontalLigament::calcStressTensor()
{
  stress_tensor_initialize();

  // #pragma omp parallel for
  //for(int ic=0;ic<numOfElm;ic++) calcStressTensor_PDL_element_spatialForm_hexa_SRI(ic,U,8,2,true);
  //for(int ic=0;ic<numOfElm;ic++) calcStressTensor_hyperFoam_element_spatialForm_hexa(ic,U,8,2,true);
  // for(int ic=0;ic<numOfElm;ic++) calcStressTensor_PDL_element_spatialForm_hexa_SRI_2018(ic,U,8,2,true);
  // for(int ic=0;ic<numOfElm;ic++) calcStressTensor_PDL_element_spatialForm_hexa_2018(ic,U,8,2,true);

  #pragma omp parallel for
  for(int ic=0;ic<numOfElm;ic++){
    switch(element[ic].meshType){
      case VTK_HEXAHEDRON:
        calcStressTensor_PDL_element_2018(ic,U);
        break;
      default:
        cout << "undefined meshType. Exit..." << endl;
        exit(1);
    }
  }

  setInnerForce();
}

// #################################################################
/**
 * @brief fiber information from TP file
 * @TODO linear hexa element only
 */
void humanPDL::PeriodontalLigament::setFiberDirection()
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
void humanPDL::PeriodontalLigament::calc_normal_quad(double (&normal)[3],double (&X)[4][3])
{
  double dXdr1[3],dXdr2[3],dXdr3[3],Jacobian,normalVector[3];
  double tmp;
  ARRAY2D<double> dNdr(4,2);
  ARRAY2D<double> dNdX(4,3);

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
void humanPDL::PeriodontalLigament::setFiberDirection_KogaModel(TextParser &tp)
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

    fibreAngle(ic) = angle * 180e0 / PI;
    for(int j=0;j<3;j++) fiberDirection_elm(ic,j)=direction[j];
  }
}

// #################################################################
/**
 * @brief arcLength
 */
double humanPDL::PeriodontalLigament::arcLength(const double xMin,const double xMax)
{
  double s = (xMax*sqrt(4e6*xMax*xMax+1e0)/2e0+asinh(2e3*xMax)/4e3)-(xMin*sqrt(4e6*xMin*xMin+1e0)/2e0+asinh(2e3*xMin)/4e3);
  return s;
}

// #################################################################
/**
 * @brief fiber information from TP file
 */
void humanPDL::PeriodontalLigament::inputFiberInfo(TextParser &tp)
{
    string str, base_label, label, inputDir;

    base_label = "/Domain";

    label = base_label + "/inputDir";
    if (!tp.getInspectedValue(label, inputDir))
    {
        cout << "data format is not set" << endl;
        exit(0);
    }

    string file1;
    fiberDirection.allocate(numOfElm, 8, 3); //gauss point
    label = base_label + "/fiberDirectionFile";
    if (!tp.getInspectedValue(label, file1))
    {
        cout << label << " is not found" << endl;
        exit(0);
    }
    file1 = inputDir + "/" + file1;
    FILE *fp;
    if ((fp = fopen(file1.c_str(), "r")) == NULL)
    {
        cout << file1 << " open error" << endl;
        exit(1);
    }
    for (int i = 0; i < numOfElm; i++)
    {
        for (int j = 0; j < 8; j++)
            fscanf(fp, "%lf %lf %lf\n", &fiberDirection(i, j, 0), &fiberDirection(i, j, 1), &fiberDirection(i, j, 2));
    }
    fclose(fp);
}

void humanPDL::PeriodontalLigament::allocatePDLvariables()
{
  lambda_ave.allocate(numOfElm,3);
  fibreStress.allocate(numOfElm);
  fiberDirection_elm.allocate(numOfElm,3); //gauss point
  angleVariation.allocate(numOfElm);
  fibreAngle.allocate(numOfElm);
}

// #################################################################
/**
 * @brief calc boundary conditions
 * @param [in] stress
 */
void humanPDL::PeriodontalLigament::export_vtu(const string &file)
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
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"volumeChangeRatio\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e\n",volumeChangeRatio(i));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"averageLambda\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e %e %e\n",lambda_ave(i,0),lambda_ave(i,1),lambda_ave(i,2));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"PrincipalStress\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e %e %e\n",sigmaEigen_Ave(i,0),sigmaEigen_Ave(i,1),sigmaEigen_Ave(i,2));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"firstPrincipalStress\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e %e %e\n",sigmaEigenVector_Ave(i,0,0),sigmaEigenVector_Ave(i,0,1),sigmaEigenVector_Ave(i,0,2));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"vonMisesStress\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  double vonMisesStress;
  for(int i=0;i<numOfElm;i++){
    vonMisesStress=pow(sigmaEigen_Ave(i,0)-sigmaEigen_Ave(i,1),2e0)
                  +pow(sigmaEigen_Ave(i,1)-sigmaEigen_Ave(i,2),2e0)
                  +pow(sigmaEigen_Ave(i,2)-sigmaEigen_Ave(i,0),2e0);
    vonMisesStress=sqrt(5e-1*vonMisesStress);
    fprintf(fp,"%e\n",vonMisesStress);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"fibreStress\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e\n",fibreStress(i));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"angleVariation\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e\n",angleVariation(i));
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
void humanPDL::PeriodontalLigament::export_vtu_Mises(const string &file)
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

  fprintf(fp,"<CellData>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Mises\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e\n",Mises(i));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"UInt8\" Name=\"Material\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%d\n",element[i].materialType);
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
void humanPDL::PeriodontalLigament::export_vtu_boundary(const string &file)
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

  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"bd\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",bd(i,0),bd(i,1),bd(i,2));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Int64\" Name=\"ibd\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%d %d %d\n",ibd(i,0),ibd(i,1),ibd(i,2));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  fprintf(fp,"<CellData>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"fiberDirection\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e %e %e\n",fiberDirection_elm(i,0),fiberDirection_elm(i,1),fiberDirection_elm(i,2));
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"Float64\" Name=\"fiberAngle\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%e\n",fibreAngle(i));
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"</CellData>\n");

  fprintf(fp,"</Piece>");
  fprintf(fp,"</UnstructuredGrid>");
  fprintf(fp,"</VTKFile>");
  fclose(fp);
}
