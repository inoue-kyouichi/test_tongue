#ifndef _ALLOCATION_H_
#define _ALLOCATION_H_

//##################################################################################
//
// allocation.h
//
// Copyright (c) 2016 Biomechanics Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
// Copyright (c) 2016 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   allocation.h
 * @brief  bdim definition Header
 * @author T.Otani, S.Ii
 */
#include <string>
#include <vector>

typedef std::vector<int> INTVECTOR1;
typedef std::vector<std::vector<int>> INTVECTOR2;
typedef std::vector<std::vector<std::vector<int>>> INTVECTOR3;
typedef std::vector<std::vector<std::vector<std::vector<int>>>> INTVECTOR4;
typedef std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> INTVECTOR5;

typedef std::vector<double> DOUBLEVECTOR1;
typedef std::vector<std::vector<double>> DOUBLEVECTOR2;
typedef std::vector<std::vector<std::vector<double>>> DOUBLEVECTOR3;
typedef std::vector<std::vector<std::vector<std::vector<double>>>> DOUBLEVECTOR4;
typedef std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> DOUBLEVECTOR5;
typedef std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> DOUBLEVECTOR6;
typedef std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>> DOUBLEVECTOR7;

typedef int*        INTARRAY1;
typedef int**       INTARRAY2;
typedef int***      INTARRAY3;
typedef int****     INTARRAY4;
typedef int*****    INTARRAY5;
typedef double*     DOUBLEARRAY1;
typedef double**    DOUBLEARRAY2;
typedef double***   DOUBLEARRAY3;
typedef double****  DOUBLEARRAY4;
typedef double***** DOUBLEARRAY5;


class Allocation {

 public:

    static INTARRAY1 allocate1dINT(const int &size,const int &guideCell=0);
    static INTARRAY2 allocate2dINT(const int &size1,const int &size2,const int &guideCell=0);
    static INTARRAY3 allocate3dINT(const int &size1,const int &size2,const int &size3,const int &guideCell=0);
    static INTARRAY4 allocate4dINT(const int &size1,const int &size2,const int &size3,const int &size4,const int &guideCell=0);
    static INTARRAY5 allocate5dINT(const int &size1,const int &size2,const int &size3,const int &size4,const int &size5,const int &guideCell=0);
    static DOUBLEARRAY1 allocate1dDOUBLE(const int &size,const int &guideCell=0);
    static DOUBLEARRAY2 allocate2dDOUBLE(const int &size1,const int &size2,const int &guideCell=0);
    static DOUBLEARRAY3 allocate3dDOUBLE(const int &size1,const int &size2,const int &size3,const int &guideCell=0);
    static DOUBLEARRAY4 allocate4dDOUBLE(const int &size1,const int &size2,const int &size3,const int &size4,const int &guideCell=0);
    static DOUBLEARRAY5 allocate5dDOUBLE(const int &size1,const int &size2,const int &size3,const int &size4,const int &size5,const int &guideCell=0);

    template <typename T>
    static void free2d(T var){
        free(var[0]);
        free(var);
    }
    template <typename T>
    static void free3d(T var){
        free(var[0][0]);
        free(var[0]);
        free(var);
    }
    template <typename T>
    static void free4d(T var){
        free(var[0][0][0]);
        free(var[0][0]);
        free(var[0]);
        free(var);
    }
    template <typename T>
    static void free5d(T var){
        free(var[0][0][0][0]);
        free(var[0][0][0]);
        free(var[0][0]);
        free(var[0]);
        free(var);
    }

};


#endif // _ALLOCATION_H_
