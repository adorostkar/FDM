//
//  myData.h
//  PDE
//
//  Created by Ali Dorostkar on 4/8/12.
//  Copyright (c) 2012 Uppsala. All rights reserved.
//
#include <iostream>
#include <stdlib.h>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <cmath>

#ifndef PDE_myData_h
#define PDE_myData_h

enum direction {
    X_DIR = 0,
    Y_DIR = 1,
    Z_DIR = 2
    };

class myData {
private:
    ///// MEMBERS /////
    unsigned int nx;
    unsigned int ny;
    unsigned int nz;
    unsigned int d_size;
    
    double D1P(int i, int j, int k, const double h, const int order, const int n1, const int n2, const int n3, const int c1, const int c2, const int c3) const;
    double D2P(int i, int j, int k, const double h, const int order, const int n1, const int n2, const int n3, const int c1, const int c2, const int c3) const;
public:
    int * dim;
    double *data;
    ///// CONSTRUCTORS - DESTRUCTOR /////
    myData(const unsigned int _nx = 1,const unsigned int _ny = 1, const unsigned int _nz = 1);
    myData(bool init,const unsigned int _nx,const unsigned int _ny, const unsigned int _nz);
    myData(const myData & md);
    virtual ~myData();
    
    ///// OPERATORS /////
    //index checking and bound correctness left to the user due to performance issue
    const double & operator()(const unsigned int i,const unsigned int j /*= 0*/, const unsigned int k = 0) const
    { return data[i + j*nx + k*nx*ny];};
    double & operator()(const unsigned int i,const unsigned int j /*= 0*/, const unsigned int k = 0)
    { return data[i + j*nx + k*nx*ny];};

    myData & operator=(const myData &);
    
    myData & operator+=(const myData &);
    myData & operator-=(const myData &);
    myData & operator*=(const myData &);
    myData & operator/=(const myData &);
    
    myData operator+(const myData &) const;
    myData operator-(const myData &) const;
    myData operator*(const myData &) const;
    myData operator/(const myData &) const;

    myData & operator+=(const double);
    myData & operator-=(const double);
    myData & operator*=(const double);
    myData & operator/=(const double);
    myData operator+(const double) const;
    myData operator-(const double) const;
    myData operator*(const double) const;
    myData operator/(const double) const;
    
    ///// METHODS /////
    const unsigned int Nx() const{return nx;};
    const unsigned int Ny() const{return ny;};
    const unsigned int Nz() const{return nz;};
    const int Data_Size() const{return d_size;};
    const double Norm();
    const double NormSquare();
    void Set_Data(const double * _data,const unsigned int _nx,const unsigned int _ny,const unsigned int _nz);

    // First and second derivatives.
    myData D1(double h,int order, direction d) const;
    myData D2(double h, int order, direction d) const;
    
    myData B_D1(double h,int order, direction d) const;
    
    myData S1(double h,int order, direction d) const;
    
    ///// FRIENDS /////
    friend myData operator*(const double a, const myData& m);
    friend myData operator/(const double a, const myData& m);
    friend myData operator+(const double a, const myData& m);
    friend myData operator-(const double a, const myData& m);
};

#endif
