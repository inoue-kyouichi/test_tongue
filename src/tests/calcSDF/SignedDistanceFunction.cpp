  #include "SignedDistanceFunction.h"
  
using namespace std;


// #################################################################
/**
 * @brief tentative
 */
void SignedDistanceFunction::calcSDF()
{
  int number=0;
  int unknown;
  int loop=0;

  mask.allocate(FEM.numOfNode);
  for(int ic=0;ic<FEM.numOfNode;ic++) mask(ic) = 1;

  for(int ic=0;ic<FEM.numOfElm;ic++){
    for(int i=0;i<knownNodes[ic].size();i++){
      mask(knownNodes[ic][i]) = 0;
    }
  }

  for(int ic=0;ic<elementToSolve.size();ic++){
    int elm = elementToSolve[ic];
    double sdf = calcSDFofUnknownNode(elm,unknown);
    if(sdf==-1){
      cout << "test" << endl;
      number++;
      continue;
    }
    SDF(unknown) = sdf;
    mask(unknown) = 0;
  }

  calc_dt();

  do{
    calcEikonalEquation();
    loop++;
  }while(loop!=1000);

  // do{
  //   int elm = elementToSolve[number];
  //   double sdf = calcSDFofUnknownNode(elm,unknown);
  //   if(sdf==-1){
  //     cout << "test" << endl;
  //     number++;
  //     continue;
  //   }

  //   //if calclation works well,
  //   SDF(unknown) = sdf;
  //   for(int i=0;i<elementsPerNodes[unknown].size();i++){
  //     knownNodes[elementsPerNodes[unknown][i]].push_back(unknown);
  //   }

  //   elementToSolve.clear();
  //   for(int ic=0;ic<FEM.numOfElm;ic++){
  //     if(knownNodes[ic].size()==3) elementToSolve.push_back(ic);
  //   }
  //   number = 0;
  //   loop++;
  //   cout << loop <<  " " <<  elementToSolve.size()<< endl;
  //   if(loop>FEM.numOfElm) break;
  // }while(elementToSolve.size()!=0);

}


// #################################################################
/**
 * @brief tentative
 */
double SignedDistanceFunction::calcSDFofUnknownNode(const int elm,int &unknown)
{
  int e1=knownNodes[elm][0];
  int e2=knownNodes[elm][1];
  int e3=knownNodes[elm][2];

  if(FEM.element[elm].node[0]!=e1 && FEM.element[elm].node[0]!=e2 && FEM.element[elm].node[0]!=e3){
    unknown = 0;
  }
  else if(FEM.element[elm].node[1]!=e1 && FEM.element[elm].node[1]!=e2 && FEM.element[elm].node[1]!=e3){
    unknown = 1;
  }
  else if(FEM.element[elm].node[2]!=e1 && FEM.element[elm].node[2]!=e2 && FEM.element[elm].node[2]!=e3){
    unknown = 2;
  }
  else if(FEM.element[elm].node[3]!=e1 && FEM.element[elm].node[3]!=e2 && FEM.element[elm].node[3]!=e3){
    unknown = 3;
  }else{
    cout << "invalid value." << endl;
    exit(1);
  }

  ARRAY2D<double> dNdr(4,3),dNdx(4,3),x_ref(4,3);
  GaussTetra gTet(1);
  ShapeFunction3D::C3D4_dNdr(dNdr,gTet.point[0][0],gTet.point[0][1],gTet.point[0][2],gTet.point[0][3]);

  double dXdr[3][3];
  for(int i=0;i<4;i++){
    for(int j=0;j<3;j++) x_ref(i,j) = FEM.x(FEM.element[elm].node[i],j);
  }
  FEM_MathTool::calc_dXdr(dXdr,dNdr,x_ref,4);
  FEM_MathTool::calc_dNdx(dNdx,dNdr,dXdr,4);
  double detJ = mathTool::calcDeterminant_3x3(dXdr);
  double volume = detJ/6e0 ;

  double dx=0e0,dy=0e0,dz=0e0;

  for(int i=0;i<4;i++){
    if(unknown==i) continue;
    dx += dNdx(i,0)*SDF(FEM.element[elm].node[i]);
    dy += dNdx(i,1)*SDF(FEM.element[elm].node[i]);
    dz += dNdx(i,2)*SDF(FEM.element[elm].node[i]);
  }

  double a = dNdx(unknown,0)*dNdx(unknown,0)+dNdx(unknown,1)*dNdx(unknown,1)+dNdx(unknown,2)*dNdx(unknown,2);
  double b = 2e0*dNdx(unknown,0)*dx + 2e0*dNdx(unknown,1)*dy + 2e0*dNdx(unknown,2)*dz;
  double c = dx*dx+dy*dy+dz*dz-1e0;
  


  if(b*b-4*a*c < 0e0) return -1e0;

  double d4_1 = (-b + sqrt(b*b-4e0*a*c)) / (2e0*a);
  double d4_2 = (-b - sqrt(b*b-4e0*a*c)) / (2e0*a);

  double sdf;
  if(d4_1>d4_2){
    sdf= d4_1;
  }else{
    sdf= d4_2;
  }

  // dx+=dNdx(unknown,0)*sdf;
  // dy+=dNdx(unknown,1)*sdf;
  // dz+=dNdx(unknown,2)*sdf;
  // printf("%e %e\n",sdf,dx*dx + dy*dy + dz*dz);


  unknown = FEM.element[elm].node[unknown];
  return sdf;
}

// #################################################################
/**
 * @brief tentative
 */
void SignedDistanceFunction::calcKnownNodes()
{
  knownNodes.resize(FEM.numOfElm);

  for(int ic=0;ic<FEM.numOfElm;ic++){
    if(FEM.element[ic].materialType==0){
      VECTOR1D<int> tmp(4);
      tmp[0] = FEM.element[ic].node[0];
      tmp[1] = FEM.element[ic].node[1];
      tmp[2] = FEM.element[ic].node[2];
      tmp[3] = FEM.element[ic].node[3];
      knownNodes[ic].push_back(tmp[0]);
      knownNodes[ic].push_back(tmp[1]);
      knownNodes[ic].push_back(tmp[2]);
      knownNodes[ic].push_back(tmp[3]);
    }
  }

  for(int ic=0;ic<FEM.numOfElm;ic++){
    if(knownNodes[ic].size()==4){
      for(int j=0;j<4;j++){
        int tmp = knownNodes[ic][j];
        for(int i=0;i<elementsPerNodes[tmp].size();i++){
         if(knownNodes[elementsPerNodes[tmp][i]].size()==4) continue;

        bool tmp_bool=true;
        for(int k=0;k<knownNodes[elementsPerNodes[tmp][i]].size();k++){
          if(knownNodes[elementsPerNodes[tmp][i]][k]==tmp) tmp_bool = false;
        }
         if(tmp_bool==true) knownNodes[elementsPerNodes[tmp][i]].push_back(tmp);
        }
      }
    }
  }

  for(int ic=0;ic<FEM.numOfElm;ic++){
    if(knownNodes[ic].size()==3) elementToSolve.push_back(ic);
  }

}

// #################################################################
/**
 * @brief tentative
 */
void SignedDistanceFunction::calcElementsPerNodes()
{
  elementsPerNodes.resize(FEM.numOfNode);

  for(int ic=0;ic<FEM.numOfElm;ic++){
    elementsPerNodes[FEM.element[ic].node[0]].push_back(ic);
    elementsPerNodes[FEM.element[ic].node[1]].push_back(ic);
    elementsPerNodes[FEM.element[ic].node[2]].push_back(ic);
    elementsPerNodes[FEM.element[ic].node[3]].push_back(ic);
  }

}

// #################################################################
/**
 * @brief tentative
 */
void SignedDistanceFunction::set1stTetra()
{
  FEM.numOfElm = FEMorg.numOfElm;
  FEM.element.resize(FEM.numOfElm);

  for(int ic=0;ic<FEM.numOfElm;ic++){
    FEM.element[ic].node.resize(4);
    FEM.element[ic].meshType = VTK_TETRA;
    for(int j=0;j<4;j++) FEM.element[ic].node[j] = FEMorg.element[ic].node[j];
  }

  ARRAY1D<int> mask_nodes(FEMorg.numOfNode);
  for(int ic=0;ic<FEMorg.numOfNode;ic++) mask_nodes(ic) = 0;
  for(int ic=0;ic<FEM.numOfElm;ic++){
    for(int j=0;j<4;j++) mask_nodes(FEM.element[ic].node[j]) = 1;
  }

  int tmp=0;
  for(int ic=0;ic<FEMorg.numOfNode;ic++) tmp+=mask_nodes(ic);
  FEM.numOfNode = tmp;
  FEM.x.allocate(FEM.numOfNode,3);

  ARRAY1D<int> tag_nodes(FEMorg.numOfNode);
  tmp=0;
  for(int ic=0;ic<FEMorg.numOfNode;ic++){
    if(mask_nodes(ic)==1){
      for(int j=0;j<3;j++) FEM.x(tmp,j) = FEMorg.x0(ic,j);
      tag_nodes(ic) = tmp;
      tmp++;
    }
  }

  for(int ic=0;ic<FEM.numOfElm;ic++){
    for(int j=0;j<4;j++){
      FEM.element[ic].node[j] = tag_nodes(FEM.element[ic].node[j]);
    }
  }

}

// #################################################################
/**
 * @brief tentative
 */
void SignedDistanceFunction::exportVTU(std::string file)
{
  FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL) {
    cout << file << " open error" << endl;
    exit(1); 
  }

  fprintf(fp,"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<UnstructuredGrid>\n");
  fprintf(fp,"<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n",FEM.numOfNode,FEM.numOfElm);
  fprintf(fp,"<Points>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<FEM.numOfNode;i++){
    fprintf(fp,"%e %e %e\n",FEM.x(i,0),FEM.x(i,1),FEM.x(i,2));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Points>\n");
  fprintf(fp,"<Cells>\n");
  fprintf(fp,"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for(int i=0;i<FEM.numOfElm;i++){
    for(int j=0;j<FEM.element[i].node.size();j++) fprintf(fp,"%d ",FEM.element[i].node[j]);
    fprintf(fp,"\n");
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num=0;
  for(int i=0;i<FEM.numOfElm;i++){
    num += FEM.element[i].node.size();
    fprintf(fp,"%d\n",num);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for(int i=0;i<FEM.numOfElm;i++) fprintf(fp,"%d\n",FEM.element[i].meshType);
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Cells>\n");

  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"SDF\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<FEM.numOfNode;i++){
    fprintf(fp,"%e\n",SDF(i));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");

  fprintf(fp,"<CellData>");
  fprintf(fp,"<DataArray type=\"UInt8\" Name=\"Material\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<FEM.numOfElm;i++){
    fprintf(fp,"%d\n",FEM.element[i].materialType);
  }
  fprintf(fp,"</DataArray>\n");

  fprintf(fp,"<DataArray type=\"UInt8\" Name=\"knownNodes\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<FEM.numOfElm;i++){
    fprintf(fp,"%d\n",knownNodes[i].size());
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");
  fprintf(fp,"</Piece>");
  fprintf(fp,"</UnstructuredGrid>");
  fprintf(fp,"</VTKFile>");
  fclose(fp);

}