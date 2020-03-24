/**
 * @file triangleSurfaceCurvature.cpp
 * @brief fem class
 * @author t. otani
 */

#include "triangleSurfaceCurvature.h"

// #################################################################
/**
 * @brief calc surface mean curvature
 */
void triangleSurfaceCurvature::modifyNormalVectorDirection()
{
  auto calcTriangleNormal=[](DOUBLEARRAY2D &x,const int i0,const int i1,const int i2,double (&normal)[3]){
    double x1[3],x2[3],area;
    for(int j=0;j<3;j++){
      x1[j] = x(i1,j)-x(i0,j);
      x2[j] = x(i2,j)-x(i0,j);
    }
    mathTool::crossProduct(x1,x2,normal,area);
    for(int j=0;j<3;j++) normal[j] /= area;
  };

  for(int ic=0;ic<numOfNode;ic++){

    double surfaceNormal[3]={0e0,0e0,0e0};
    for(int i=0;i<neighborElements[ic].size();i++){

      int elementNumber = neighborElements[ic][i];
      int i1,i2,element[3]={ie(elementNumber,0),ie(elementNumber,1),ie(elementNumber,2)};

      if(element[0]==ic){
        i1 = element[1]; i2 = element[2];
      }else if(element[1]==ic){
        i1 = element[2]; i2 = element[0];
      }else if(element[2]==ic){
        i1 = element[0]; i2 = element[1];
      }
      double normal_tmp[3];
      calcTriangleNormal(x,ic,i1,i2,normal_tmp);
      for(int j=0;j<3;j++) surfaceNormal[j] += normal_tmp[j];
    } 
    for(int j=0;j<3;j++) surfaceNormal[j] /= neighborElements[ic].size();

    double tmp = surfaceNormal[0]*normal(ic,0)+surfaceNormal[1]*normal(ic,1)+surfaceNormal[2]*normal(ic,2);
    if(tmp<0e0){
      for(int j=0;j<3;j++) normal(ic,j) = -normal(ic,j);
      meanCurvature(ic) = -1e0*meanCurvature(ic);
    } 
  }

}

// #################################################################
/**
 * @brief calc surface mean curvature
 */
void triangleSurfaceCurvature::calcSurfaceMeanCurvature()
{
  auto directionVector=[](DOUBLEARRAY2D &x,const int i1,const int i2,double (&x1)[3]){
    for(int j=0;j<3;j++) x1[j] = x(i2,j)-x(i1,j);
  };

  auto innerProduct=[](const double (&x1)[3],const double (&x2)[3]){
    return x1[0]*x2[0]+x1[1]*x2[1]+x1[2]*x2[2];
  };

  auto calcTriangleArea=[](const double (&x1)[3],const double (&x2)[3]){
    double x3[3],area=0e0;
    mathTool::crossProduct(x1,x2,x3,area);
    return 5e-1*area;
  };

  for(int ic=0;ic<numOfNode;ic++){

    double area = 0e0;
    for(int j=0;j<3;j++) normal(ic,j) = 0e0;

    for(int i=0;i<neighborElements[ic].size();i++){

      int elementNumber = neighborElements[ic][i];
      int i1,i2,element[3]={ie(elementNumber,0),ie(elementNumber,1),ie(elementNumber,2)};

      if(element[0]==ic){
        i1 = element[1]; i2 = element[2];
      }else if(element[1]==ic){
        i1 = element[2]; i2 = element[0];
      }else if(element[2]==ic){
        i1 = element[0]; i2 = element[1];
      }

      double x10[3],x12[3];
      directionVector(x,i1,ic,x10);
      directionVector(x,i1,i2,x12);
      double cot012 = cot(x10,x12);

      double x20[3],x21[3];
      directionVector(x,i2,ic,x20);
      directionVector(x,i2,i1,x21);
      double cot120 = cot(x20,x21);

      for(int j=0;j<3;j++){
        normal(ic,j) += 5e-1 * (cot012 * x20[j] + cot120 * x10[j]);
      }

      if(acuteTriangle(x10,x20)==true){
        for(int j=0;j<3;j++) area += (innerProduct(x10,x10)*cot120+innerProduct(x20,x20)*cot012)/8e0;
      }else{
        double x0c[3],x1[3],x2[3];
        for(int j=0;j<3;j++){
          x0c[j] = -5e-1 * (x10[j] + x20[j]);
          x1[j] = -5e-1 * x10[j];
          x2[j] = -5e-1 * x20[j];
        }
        area += calcTriangleArea(x1,x0c);
        area += calcTriangleArea(x2,x0c);
      }
    }

    for(int j=0;j<3;j++) normal(ic,j) /= area;
    double tmp = sqrt(normal(ic,0)*normal(ic,0)+normal(ic,1)*normal(ic,1)+normal(ic,2)*normal(ic,2));
    meanCurvature(ic) = tmp;
    for(int j=0;j<3;j++) normal(ic,j) = normal(ic,j)/tmp;
  }

  modifyNormalVectorDirection();
}

// #################################################################
/**
 * @brief calc triangle area in 1-ring
 */
bool triangleSurfaceCurvature::acuteTriangle(double (&x1)[3],double (&x2)[3])
{
  double x01[3],x02[3],x12[3];
  for(int j=0;j<3;j++){
    x01[j] = -x1[j];
    x02[j] = -x2[j];
    x12[j] = x01[j] - x02[j];
  }

  auto vectorLength=[](const double (&x)[3]){
    return std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  };

  double edgeLength[3]={vectorLength(x01),vectorLength(x02),vectorLength(x12)};
  std::sort(edgeLength,edgeLength+3);

  double a2 = edgeLength[0]*edgeLength[0];
  double b2 = edgeLength[1]*edgeLength[1];
  double c2 = edgeLength[2]*edgeLength[2];

  double c[3];

  if(a2+b2>c2) return true;

  return false;
}

// #################################################################
/**
 * @brief calc cotangent of angle between edge i0-i1 and i0-i2
 */
double triangleSurfaceCurvature::cot(const double (&x1)[3],const double (&x2)[3])
{
  auto innerProduct=[](const double (&x1)[3],const double (&x2)[3]){
    return x1[0]*x2[0]+x1[1]*x2[1]+x1[2]*x2[2];
  };

  auto crossProductMag=[](const double (&x1)[3],const double (&x2)[3]){
    double x3[3],tmp;
    mathTool::crossProduct(x1,x2,x3,tmp);
    return tmp;
  };

  double tmp = sqrt(innerProduct(x1,x1))*sqrt(innerProduct(x2,x2));
  double cosA = innerProduct(x1,x2)/tmp;
  double sinA = crossProductMag(x1,x2)/tmp;

  double cotA;
  if(fabs(sinA)<1e-15){
    std::cout << "triangle surface may be collapsed (sin<1e-15). Please check. Exit..." << std::endl;
     exit(1);
  }else{
    cotA = cosA/sinA;
  }
  return cotA;
}

// #################################################################
/**
 * @brief export surface normal
 */
void triangleSurfaceCurvature::exportSurfaceNormal(DOUBLEARRAY2D &normal_export)
{
  normal_export.allocate(numOfNode,3);

  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) normal_export(i,j) = normal(i,j);
  }
}

// #################################################################
/**
 * @brief export surface mean curvature
 */
void triangleSurfaceCurvature::exportSurfaceMeanCurvature(DOUBLEARRAY1D &meanCurvature_export)
{
  meanCurvature_export.allocate(numOfNode);

  for(int i=0;i<numOfNode;i++){
    meanCurvature_export(i) = meanCurvature(i);
  }
}

// #################################################################
/**
 * @brief calc adjacent nodes
 */
void triangleSurfaceCurvature::calc_adjacent_nodes()
{
  int tmp;
  bool b1;

  for(int ic=0;ic<numOfNode;ic++){
    for(int i=0,n=neighborElements[ic].size();i<n;i++){
      for(int j=0;j<3;j++){
        tmp = ie(neighborElements[ic][i],j);
        b1 = true;
        for(int k=0,n=neighborNodes[ic].size();k<n;k++){
          if(neighborNodes[ic][k]==tmp){
            b1 = false;
            break;
          }
        }
        if(b1==false) continue;
        neighborNodes[ic].push_back(tmp);
      }
    }
    sort(neighborNodes[ic].begin(),neighborNodes[ic].end());
  }
}

// #################################################################
/**
 * @brief calc adjacent elements
 */
void triangleSurfaceCurvature::calc_adjacent_elements()
{
  int tmp;
  for(int ic=0;ic<numOfElm;ic++){
    for(int j=0,n=3;j<n;j++){
      tmp = ie(ic,j);
      neighborElements[tmp].push_back(ic);
    }
  }
}
