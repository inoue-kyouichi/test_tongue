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




#endif // _ALLOCATION_H_
