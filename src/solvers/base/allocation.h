#ifndef _ALLOCATION_H_
#define _ALLOCATION_H_

//##################################################################################
//
// allocation.h
//
// Copyright (c) 2019 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   allocation.h
 * @brief  bdim definition Header
 * @author T.Otani
 */
#include <string>
#include <iostream>

template <typename T>
class ARRAY1D
{
public:
    ARRAY1D()
    {
        data = new T[1];
    }
    ARRAY1D(const int X)
    {
        nx = X;
        data = new T[nx];
    }
    ~ARRAY1D()
    {
        delete[] data;
    }
    T &operator()(const int i)
    {
        return data[i];
    }
    void allocate(const int X)
    {
        delete[] data;
        nx = X;
        data = new T[nx];
    }

private:
    T *data;
    int nx;
};

template <typename T>
class ARRAY2D
{
public:
    ARRAY2D()
    {
        data = new T[1];
    }
    ARRAY2D(const int X, const int Y)
    {
        nx = X;
        ny = Y;
        data = new T[nx*ny];
    }
    ~ARRAY2D()
    {
        delete[] data;
    }
    T &operator()(const int i, const int j)
    {
        return data[i*ny+j];
    }
    void allocate(const int X, const int Y)
    {
        delete[] data;
        nx = X;
        ny = Y;
        data = new T[nx*ny];
    }
private:
    T *data;
    int nx;
    int ny;
};

template <typename T>
class ARRAY3D
{
public:
    ARRAY3D()
    {
        data = new T[1];
    }
    ARRAY3D(const int X, const int Y, const int Z)
    {
        nz = Z;
        nx = X;
        ny = Y;
        data = new T[nx*ny*nz];
    }
    ~ARRAY3D()
    {
        delete[] data;
    }
    T &operator()(const int i, const int j, const int k)
    {
        return data[i*ny*nz+j*nz+k];
    }
    void allocate(const int X, const int Y, const int Z)
    {
        delete[] data;
        nz = Z;
        nx = X;
        ny = Y;
        data = new T[nx*ny*nz];
    }

private:
    T *data;
    int nx;
    int ny;
    int nz;
};

template <typename T>
class ARRAY4D
{
public:
    ARRAY4D()
    {
        data = new T[1];
    }
    ARRAY4D(const int H, const int X, const int Y, const int Z)
    {
        nh = H;
        nz = Z;
        nx = X;
        ny = Y;
        data = new T[nx*ny*nz*nh];
    }
    ~ARRAY4D()
    {
        delete[] data;
    }
    T &operator()(const int h,const int i, const int j, const int k)
    {
        return data[h*nx*ny*nz+i*ny*nz+j*nz+k];
    }
    void allocate(const int H, const int X, const int Y, const int Z)
    {
        delete[] data;
        nh = H;
        nz = Z;
        nx = X;
        ny = Y;
        data = new T[nh*nx*ny*nz];
    }

private:
    T *data;
    int nx;
    int ny;
    int nz;
    int nh;
};

template <typename T>
class ARRAY5D
{
public:
    ARRAY5D()
    {
        data = new T[1];
    }
    ARRAY5D(const int L,const int H, const int X, const int Y, const int Z)
    {
        nl = L;
        nh = H;
        nz = Z;
        nx = X;
        ny = Y;
        data = new T[nl*nx*ny*nz*nh];
    }
    ~ARRAY5D()
    {
        delete[] data;
    }
    T &operator()(const int l,const int h,const int i, const int j, const int k)
    {
        return data[l*nh*nx*ny*nz+h*nx*ny*nz+i*ny*nz+j*nz+k];
    }
    void allocate(const int L,const int H, const int X, const int Y, const int Z)
    {
        delete[] data;
        nl = L;
        nh = H;
        nz = Z;
        nx = X;
        ny = Y;
        data = new T[nl*nh*nx*ny*nz];
    }

private:
    T *data;
    int nx;
    int ny;
    int nz;
    int nh;
    int nl;
};


class INTARRAY1D
{
public:
    INTARRAY1D()
    {
        data=new int[1];
    }
    INTARRAY1D(const int X)
    {
        nx = X;
        data = new int[nx];
    }
    ~INTARRAY1D()
    {
        delete[] data;
    }
    int &operator()(const int i)
    {
        return data[i];
    }
    void allocate(const int X)
    {
        delete[] data;
        nx = X;
        data = new int[nx];
    }
private:
    int *data;
    int nx;
};

class INTARRAY2D
{
public:
    INTARRAY2D()
    {
        data=new int[1];
    }
    INTARRAY2D(const int X, const int Y)
    {
        nx = X;
        ny = Y;
        data = new int[nx*ny];
    }
    ~INTARRAY2D()
    {
        delete[] data;
    }
    int &operator()(const int i, const int j)
    {
        return data[i*ny+j];
    }
    void allocate(const int X, const int Y)
    {
        delete[] data;
        nx = X;
        ny = Y;
        data = new int[nx*ny];
    }
    void importData(const std::string &file);
    void exportData(const std::string &file);

private:
    int *data;
    int nx;
    int ny;
};

class INTARRAY3D
{
public:
    INTARRAY3D()
    {
        data = new int[1];
    }
    INTARRAY3D(const int X, const int Y, const int Z)
    {
        nz = Z;
        nx = X;
        ny = Y;
        data = new int[nx*ny*nz];
    }
    ~INTARRAY3D()
    {
        delete[] data;
    }
    int &operator()(const int i, const int j, const int k)
    {
        return data[i*ny*nz+j*nz+k];
    }
    void allocate(const int X, const int Y, const int Z)
    {
        delete[] data;
        nz = Z;
        nx = X;
        ny = Y;
        data = new int[nx*ny*nz];
    }
private:
    int *data;
    int nx;
    int ny;
    int nz;
};

class DOUBLEARRAY1D
{
public:
    DOUBLEARRAY1D()
    {
        data = new double[1];
    }
    DOUBLEARRAY1D(const int X)
    {
        nx = X;
        data = new double[nx];
    }
    ~DOUBLEARRAY1D()
    {
        delete[] data;
    }
    double &operator()(const int i)
    {
        return data[i];
    }
    void allocate(const int X)
    {
        delete[] data;
        nx = X;
        data = new double[nx];
    }

private:
    double *data;
    int nx;
};

class DOUBLEARRAY2D
{
public:
    DOUBLEARRAY2D()
    {
        data = new double[1];
    }
    DOUBLEARRAY2D(const int X, const int Y)
    {
        nx = X;
        ny = Y;
        data = new double[nx*ny];
    }
    ~DOUBLEARRAY2D()
    {
        delete[] data;
    }
    double &operator()(const int i, const int j)
    {
        return data[i*ny+j];
    }
    void allocate(const int X, const int Y)
    {
        delete[] data;
        nx = X;
        ny = Y;
        data = new double[nx*ny];
    }
private:
    double *data;
    int nx;
    int ny;
};

class DOUBLEARRAY3D
{
public:
    DOUBLEARRAY3D()
    {
        data = new double[1];
    }
    DOUBLEARRAY3D(const int X, const int Y, const int Z)
    {
        nz = Z;
        nx = X;
        ny = Y;
        data = new double[nx*ny*nz];
    }
    ~DOUBLEARRAY3D()
    {
        delete[] data;
    }
    double &operator()(const int i, const int j, const int k)
    {
        return data[i*ny*nz+j*nz+k];
    }
    void allocate(const int X, const int Y, const int Z)
    {
        delete[] data;
        nz = Z;
        nx = X;
        ny = Y;
        data = new double[nx*ny*nz];
    }

private:
    double *data;
    int nx;
    int ny;
    int nz;
};

class DOUBLEARRAY4D
{
public:
    DOUBLEARRAY4D()
    {
        data = new double[1];
    }
    DOUBLEARRAY4D(const int H, const int X, const int Y, const int Z)
    {
        nh = H;
        nz = Z;
        nx = X;
        ny = Y;
        data = new double[nx*ny*nz*nh];
    }
    ~DOUBLEARRAY4D()
    {
        delete[] data;
    }
    double &operator()(const int h,const int i, const int j, const int k)
    {
        return data[h*nx*ny*nz+i*ny*nz+j*nz+k];
    }
    void allocate(const int H, const int X, const int Y, const int Z)
    {
        delete[] data;
        nh = H;
        nz = Z;
        nx = X;
        ny = Y;
        data = new double[nh*nx*ny*nz];
    }

private:
    double *data;
    int nx;
    int ny;
    int nz;
    int nh;
};

class DOUBLEARRAY5D
{
public:
    DOUBLEARRAY5D()
    {
        data = new double[1];
    }
    DOUBLEARRAY5D(const int L,const int H, const int X, const int Y, const int Z)
    {
        nl = L;
        nh = H;
        nz = Z;
        nx = X;
        ny = Y;
        data = new double[nl*nx*ny*nz*nh];
    }
    ~DOUBLEARRAY5D()
    {
        delete[] data;
    }
    double &operator()(const int l,const int h,const int i, const int j, const int k)
    {
        return data[l*nh*nx*ny*nz+h*nx*ny*nz+i*ny*nz+j*nz+k];
    }
    void allocate(const int L,const int H, const int X, const int Y, const int Z)
    {
        delete[] data;
        nl = L;
        nh = H;
        nz = Z;
        nx = X;
        ny = Y;
        data = new double[nl*nh*nx*ny*nz];
    }

private:
    double *data;
    int nx;
    int ny;
    int nz;
    int nh;
    int nl;
};


#endif // _ALLOCATION_H_
