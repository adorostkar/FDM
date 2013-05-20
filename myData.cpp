//
//  myData.cpp
//  PDE
//
//  Created by Ali Dorostkar on 4/8/12.
//  Copyright (c) 2012 Uppsala. All rights reserved.
//

#include "myData.h"
#define MG(i,j,k,c1,c2,c3) data[c1*(i) + c2*(j) + c3*(k)]

///// CONSTRUCTORS - DESTRUCTOR /////
myData::myData(const unsigned int _nx,const unsigned int _ny, const unsigned int _nz){
    nx = _nx; ny = _ny; nz = _nz; d_size = nx*ny*nz;
    dim = (int*)calloc(3, sizeof(int));
    dim[0] = nx; dim[1] = ny; dim[2] = nz;
    data = (double*)calloc(d_size, sizeof(double));
}
myData::myData(bool init,const unsigned int _nx,const unsigned int _ny, const unsigned int _nz){
    nx = _nx; ny = _ny; nz = _nz; d_size = nx*ny*nz;
    dim = (int*)calloc(3, sizeof(int));
    dim[0] = nx; dim[1] = ny; dim[2] = nz;
    data = (double*)malloc(d_size*sizeof(double));
}

myData::myData(const myData & md){
    nx = md.nx;
    ny = md.ny;
    nz = md.nz;
    dim = (int*)calloc(3, sizeof(int));
    dim[0] = nx; dim[1] = ny; dim[2] = nz;
    d_size = md.d_size;
    
    data = (double*)malloc(d_size*sizeof(double));
    memcpy(data, md.data, d_size*sizeof(double));
}

myData::~myData(){
    free(data);
    free(dim);
}

///// OPERATORS /////

myData & myData::operator=(const myData & md){
    if (this == &md)
        return *this;
    
    nx = md.nx;
    ny = md.ny;
    nz = md.nz;
    dim[0] = nx; dim[1] = ny; dim[2] = nz;
    d_size = md.d_size;
    
    free(data);
    data = (double*)malloc(d_size*sizeof(double));
    memcpy(data, md.data, d_size*sizeof(double));
    return *this;
}

myData & myData::operator+=(const myData & md){
    if (nx != md.nx || ny != md.ny || nz != md.nz) {
        throw std::invalid_argument("+= not applicable. Sizes of the two operand don't match");
    }
#pragma omp parallel for
    for (int i = 0; i < d_size; i++) {
        data[i] += md.data[i];
    }
    return *this;
}
myData & myData::operator-=(const myData & md){
    if (nx != md.nx || ny != md.ny || nz != md.nz) {
        throw std::invalid_argument("-= not applicable. Sizes of the two operand don't match");
    }
#pragma omp parallel for
    for (int i = 0; i < d_size; i++) {
        data[i] -= md.data[i];
    }
    return *this;
}
myData & myData::operator*=(const myData & md){
    if (nx != md.nx || ny != md.ny || nz != md.nz) {
        throw std::invalid_argument("*= not applicable. Sizes of the two operand don't match");
    }
#pragma omp parallel for
    for (int i = 0; i < d_size; i++) {
        data[i] *= md.data[i];
    }
    return *this;
}
myData & myData::operator/=(const myData & md){
    if (nx != md.nx || ny != md.ny || nz != md.nz) {
        throw std::invalid_argument("/= not applicable. Sizes of the two operand don't match");
    }
#pragma omp parallel for
    for (int i = 0; i < d_size; i++) {
        data[i] /= md.data[i];
    }
    return *this;
}

myData myData::operator+(const myData & md) const{
    if (nx != md.nx || ny != md.ny || nz != md.nz) {
        throw std::invalid_argument("+ not applicable. Sizes of the two operand don't match");
    }
    myData res(false,nx,ny,nz);
#pragma omp parallel for
    for (int i = 0; i < d_size; i++) {
        res.data[i] = data[i] + md.data[i];
    }
    return res;
}
myData myData::operator-(const myData & md) const{
    if (nx != md.nx || ny != md.ny || nz != md.nz) {
        throw std::invalid_argument("- not applicable. Sizes of the two operand don't match");
    }
    myData res(false,nx,ny,nz);
#pragma omp parallel for
    for (int i = 0; i < d_size; i++) {
        res.data[i] = data[i] - md.data[i];
    }
    return res;
}
myData myData::operator*(const myData & md) const{
    if (nx != md.nx || ny != md.ny || nz != md.nz) {
        throw std::invalid_argument("* not applicable. Sizes of the two operand don't match");
    }
    myData res(false,nx,ny,nz);
#pragma omp parallel for
    for (int i = 0; i < d_size; i++) {
        res.data[i] = data[i] * md.data[i];
    }
    return res;
}
myData myData::operator/(const myData & md) const{
    if (nx != md.nx || ny != md.ny || nz != md.nz) {
        throw std::invalid_argument("/ not applicable. Sizes of the two operand don't match");
    }
    myData res(false,nx,ny,nz);
#pragma omp parallel for
    for (int i = 0; i < d_size; i++) {
        res.data[i] = data[i] / md.data[i];
    }
    return res;
}

myData & myData::operator+=(const double a){
#pragma omp parallel for
    for (int i = 0; i < d_size; i++) {
        data[i] += a;
    }
    return *this;
}
myData & myData::operator-=(const double a){
#pragma omp parallel for
    for (int i = 0; i < d_size; i++) {
        data[i] -= a;
    }
    return *this;
}
myData & myData::operator*=(const double a){
#pragma omp parallel for
    for (int i = 0; i < d_size; i++) {
        data[i] *= a;
    }
    return *this;
}
myData & myData::operator/=(const double a){
    for (int i = 0; i < d_size; i++) {
        data[i] /= a;
    }
    return *this;
}
myData myData::operator+(const double a) const{
    myData res(false,nx,ny,nz);
#pragma omp parallel for
    for (int i = 0; i < d_size; i++) {
        res.data[i] = data[i] + a;
    }
    return res;
}
myData myData::operator-(const double a) const{
    myData res(false,nx,ny,nz);
#pragma omp parallel for
    for (int i = 0; i < d_size; i++) {
        res.data[i] = data[i] - a;
    }
    return res;
}
myData myData::operator*(const double a) const{
    myData res(false,nx,ny,nz);
#pragma omp parallel for
    for (int i = 0; i < d_size; i++) {
        res.data[i] = data[i] * a;
    }
    return res;
}
myData myData::operator/(const double a) const{
    myData res(nx,ny,nz);
#pragma omp parallel for
    for (int i = 0; i < d_size; i++) {
        res.data[i] = data[i] / a;
    }
    return res;
}

///// METHODS /////
void myData::Set_Data(const double * _data,const unsigned int _nx,const unsigned int _ny,const unsigned int _nz){
    nx = _nx;
    ny = _ny;
    nz = _nz;
    dim[0] = nx; dim[1] = ny; dim[2] = nz;
    d_size = nx*ny*nz;
    free(data);
    
    data = (double*)malloc(d_size*sizeof(double));
    memcpy(data, _data, d_size*sizeof(double));
}

const double myData::NormSquare(){
    double res = 0;
#pragma omp parallel for
    for (int i = 0; i < nx*ny*nz; i++) {
        res += data[i]*data[i];
    }
    return res;
}

const double myData::Norm(){
    double res = 0;
#pragma omp parallel for
    for (int i = 0; i < nx*ny*nz; i++) {
        res += data[i]*data[i];
    }
    return sqrt(res);
}

///// FRIENDS /////
myData operator*(const double a, const myData& md){
    return md*a;
}
myData operator/(const double a, const myData& md){
    myData res(false,md.nx,md.ny,md.nz);
#pragma omp parallel for
    for (int i = 0; i < md.nx*md.ny*md.nz; i++) {
        res.data[i] = a / md.data[i];
    }
    return res;
}
myData operator+(const double a, const myData& md){
    return md+a;
}
myData operator-(const double a, const myData& md){
    myData res(false,md.nx,md.ny,md.nz);
#pragma omp parallel for
    for (int i = 0; i < md.nx*md.ny*md.nz; i++) {
        res.data[i] = a - md.data[i];
    }
    return res;
}


myData myData::D1(double h,int order, direction d) const{
    // The number of elements that is passed in the array if index in the direction is added by one
    int c1,c2,c3;
    // number of points in each directio but regardless of name.
    int n1,n2,n3;
    int i,j,k;
    myData u_1(nx,ny,nz);
    
    switch (d) {
        case X_DIR:
            c1 = 1;     c2 = nx; c3 = nx*ny;
            n1 = nx;    n2 = ny; n3 = nz;
            break;
        case Y_DIR:
            c1 = nx; c2 = nx*ny;    c3 = 1;
            n1 = ny; n2 = nz;       n3 = nx;
            break;
        case Z_DIR:
            c1 = nx*ny; c2 = 1;     c3 = nx;
            n1 = nz;    n2 = nx;    n3 = ny;
            break;
    }
    
    if(( order == 2 && n1 < 2 ) || ( order == 4 && n1 < 12 ) || ( order == 6 && n1 < 16) )
        return u_1;
    
    for (k = 0; k < n3; k++) {
        for (j = 0; j < n2; j++) {
            for (i = 0; i < n1; i++) {
                u_1.MG(i,j,k,c1,c2,c3) = D1P(i, j, k, h, order, n1, n2, n3, c1, c2, c3);
            }
        }
    }
    return u_1;
}

myData myData::D2(double h, int order, direction d) const{
    // The number of elements that is passed in the array if index in the direction is added by one
    int c1,c2,c3;
    // number of points in each directio but regardless of name.
    int n1,n2,n3;
    int i,j,k;
    
    myData u_2(nx,ny,nz);
    
    switch (d) {
        case X_DIR:
            c1 = 1;     c2 = nx; c3 = nx*ny;
            n1 = nx;    n2 = ny; n3 = nz;
            break;
        case Y_DIR:
            c1 = nx; c2 = nx*ny;    c3 = 1;
            n1 = ny; n2 = nz;       n3 = nx;
            break;
        case Z_DIR:
            c1 = nx*ny; c2 = 1;     c3 = nx;
            n1 = nz;    n2 = nx;    n3 = ny;
            break;
    }
    
    if(( order == 2 && n1 < 3 ) || ( order == 4 && n1 < 12 ) || ( order == 6 && n1 < 16) )
        return u_2;
    
    for (k = 0; k < n3; k++) {
        for (j = 0; j < n2; j++) {
            for (i = 0; i < n1; i++) {
                u_2.MG(i,j,k,c1,c2,c3) = D2P(i, j, k, h, order, n1, n2, n3, c1, c2, c3);
            }
        }
    }
    return u_2;
}

myData myData::B_D1(double h,int order, direction d) const{
    // The number of elements that is passed in the array if index in the direction is added by one
    int c1,c2,c3;
    // number of points in each directio but regardless of name.
    int n1,n2,n3;
    int j,k;
    myData u_1(nx,ny,nz);
    
    switch (d) {
        case X_DIR:
            c1 = 1;     c2 = nx; c3 = nx*ny;
            n1 = nx;    n2 = ny; n3 = nz;
            break;
        case Y_DIR:
            c1 = nx; c2 = nx*ny;    c3 = 1;
            n1 = ny; n2 = nz;       n3 = nx;
            break;
        case Z_DIR:
            c1 = nx*ny; c2 = 1;     c3 = nx;
            n1 = nz;    n2 = nx;    n3 = ny;
            break;
    }
    
    if(( order == 2 && n1 < 2 ) || ( order == 4 && n1 < 12 ) || ( order == 6 && n1 < 16) )
        return u_1;
    
    for (k = 0; k < n3; k++) {
        for (j = 0; j < n2; j++) {
            u_1.MG(0,j,k,c1,c2,c3) = - D1P(0, j, k, h, order, n1, n2, n3, c1, c2, c3);
            u_1.MG(n1-1,j,k,c1,c2,c3) = D1P(n1-1, j, k, h, order, n1, n2, n3, c1, c2, c3);
        }
    }
    
    return u_1;
}

myData myData::S1(double h,int order, direction d) const{
    // The number of elements that is passed in the array if index in the direction is added by one
    int c1,c2,c3;
    // number of points in each directio but regardless of name.
    int n1,n2,n3;
    int j,k;
    myData u_1(nx,ny,nz);
    
    switch (d) {
        case X_DIR:
            c1 = 1;     c2 = nx; c3 = nx*ny;
            n1 = nx;    n2 = ny; n3 = nz;
            break;
        case Y_DIR:
            c1 = nx; c2 = nx*ny;    c3 = 1;
            n1 = ny; n2 = nz;       n3 = nx;
            break;
        case Z_DIR:
            c1 = nx*ny; c2 = 1;     c3 = nx;
            n1 = nz;    n2 = nx;    n3 = ny;
            break;
    }
    
    if(( order == 2 && n1 < 2 ) || ( order == 4 && n1 < 12 ) || ( order == 6 && n1 < 16) )
        return u_1;
    
    for (k = 0; k < n3; k++) {
        for (j = 0; j < n2; j++) {
            u_1.MG(0,j,k,c1,c2,c3) = D1P(0, j, k, h, order, n1, n2, n3, c1, c2, c3);
            u_1.MG(n1-1,j,k,c1,c2,c3) = D1P(n1-1, j, k, h, order, n1, n2, n3, c1, c2, c3);
        }
    }
    return u_1;
}

// Computes the second derivative of a function in each point.
double myData::D2P(int i, int j, int k, const double h, const int order, const int n1, const int n2, const int n3, const int c1, const int c2, const int c3) const{
    int i0 = i, i1 = i, i2 = i, i3 = i, i4 = i, i5 = i;
    int i6 = i, i7 = i, i8 = i, i9 = i, i10 = i;
    double a0 = 0, a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0;
    double a6 = 0, a7 = 0, a8 = 0, a9 = 0, a10 = 0;
	if(order == 2){
		if(i == 0 || i == n1 - 1)
			return 0;
		else{
			a0 = 1;                     i0 = i - 1;
			a1 = -2;                    i1 = i;
			a2 = 1;                     i2 = i + 1;
		}
	}
    else if (order == 4){
		if(i == 0){
			a0 = 1.319753349375470;		i0 = 0;
			a1 =-2.279013397501878;     i1 = 1;
			a2 =-0.081479903747182;     i2 = 2;
			a3 = 1.720986602498121;		i3 = 3;
			a4 =-0.680246650624530;     i4 = 4;
		}else if( i == 1){
			a0 = 0.901112191978006;		i0 = 0;
			a1 =-1.593787843034701;     i1 = 1;
			a2 = 0.364029452358750;		i2 = 2;
			a3 = 0.459516781351903;		i3 = 3;
			a4 =-0.141531507531278;     i4 = 4;
			a5 = 0.010660924877321;		i5 = 5;
		}else if( i == 2){
			a0 = 0.023650158582586;		i0 = 0;
			a1 = 0.796439026381821;		i1 = 1;
			a2 =-1.422257691353143;     i2 = 2;
			a3 = 0.251637329942644;		i3 = 3;
			a4 = 0.459491515733928;		i4 = 4;
			a5 =-0.108960339287835;     i5 = 5;
		}else if(i == 3){
			a0 =-0.120536120137894;     i0 = 0;
			a1 = 0.518872114613017;		i1 = 1;
			a2 = 0.129872742926870;		i2 = 2;
			a3 =-1.297489715079774;     i3 = 3;
			a4 = 0.732553343616339;		i4 = 4;
			a5 = 0.036727634061442;		i5 = 5;
		}else if (i == 4){
			a0 = 0.044476945566742;		i0 = 0;
			a1 =-0.212643239032514;     i1 = 1;
			a2 = 0.315544004411490;		i2 = 2;
			a3 = 0.974717461344343;		i3 = 3;
			a4 =-2.358138203677955;     i4 = 4;
			a5 = 1.326302527439042;		i5 = 5;
			a6 =-0.090259496051147;     i6 = 6;
		}else if(i == 5){
            a0 = 0;                     i0 = 0;
			a1 = 0.014642919234304;		i1 = 1;
			a2 =-0.068404586240661;     i2 = 2;
			a3 = 0.044675228394891;		i3 = 3;
			a4 = 1.212486564140966;		i4 = 4;
			a5 =-2.441108988900194;     i5 = 5;
			a6 = 1.320222787595407;		i6 = 6;
			a7 =-0.082513924224713;     i7 = 7;
		}else if(i == n1-6){
			a0 =-0.082513924224713;     i0 = n1-8;
			a1 = 1.320222787595407;		i1 = n1-7;
			a2 =-2.441108988900194;     i2 = n1-6;
			a3 = 1.212486564140966;		i3 = n1-5;
			a4 = 0.044675228394891;		i4 = n1-4;
			a5 =-0.068404586240661;     i5 = n1-3;
			a6 = 0.014642919234304;		i6 = n1-2;
		}else if(i == n1-5){
			a0 =-0.090259496051147;     i0 = n1-7;
			a1 = 1.326302527439042;		i1 = n1-6;
			a2 =-2.358138203677955;     i2 = n1-5;
			a3 = 0.974717461344343;		i3 = n1-4;
			a4 = 0.315544004411490;		i4 = n1-3;
			a5 =-0.212643239032514;     i5 = n1-2;
			a6 = 0.044476945566742;		i6 = n1-1;
		}else if(i == n1-4){
			a0 = 0.036727634061442;		i0 = n1-6;
			a1 = 0.732553343616339;		i1 = n1-5;
			a2 =-1.297489715079774;     i2 = n1-4;
			a3 = 0.129872742926870;		i3 = n1-3;
			a4 = 0.518872114613017;		i4 = n1-2;
			a5 =-0.120536120137894;     i5 = n1-1;
		}else if(i == n1-3){
			a0 =-0.108960339287835;     i0 = n1-6;
			a1 = 0.459491515733928;		i1 = n1-5;
			a2 = 0.251637329942644;		i2 = n1-4;
			a3 =-1.422257691353143;     i3 = n1-3;
			a4 = 0.796439026381821;		i4 = n1-2;
			a5 = 0.023650158582586;		i5 = n1-1;
		}else if(i == n1-2){
			a0 = 0.010660924877321;		i0 = n1-6;
			a1 =-0.141531507531278;     i1 = n1-5;
			a2 = 0.459516781351903;		i2 = n1-4;
			a3 = 0.364029452358750;		i3 = n1-3;
			a4 =-1.593787843034701;     i4 = n1-2;
			a5 = 0.901112191978006;		i5 = n1-1;
		}else if(i == n1-1){
			a0 =-0.680246650624530;     i0 = n1-5;
			a1 = 1.720986602498121;		i1 = n1-4;
			a2 =-0.081479903747182;     i2 = n1-3;
			a3 =-2.279013397501878;     i3 = n1-2;
			a4 = 1.319753349375470;		i4 = n1-1;
		}else{
			a0 =-0.083333333333333;     i0 = i - 2;
			a1 = 1.333333333333333;		i1 = i - 1;
			a2 =-2.5;                   i2 = i;
			a3 = 1.333333333333333;		i3 = i + 1;
			a4 =-0.083333333333333;     i4 = i + 2;
		}
	}
    else{
        if(i == 0){
            a0 = 1.677132509479325;     i0 = 0;
            a1 =-3.038214359671524;     i1 = 1;
            a2 =-0.049249177165591;     i2 = 2;
            a3 = 2.036490115791102;     i3 = 3;
            a4 = 0.411180670145605;     i4 = 4;
            a5 =-1.606558237520481;     i5 = 5;
            a6 = 0.569218478941565;     i6 = 6;
        }else if(i == 1){
            a0 = 0.884778483741509;     i0 = 0;
            a1 = -1.558734408191758;    i1 = 1;
            a2 = 0.443635669435212;     i2 = 2;
            a3 = 0.112257511414158;     i3 = 3;
            a4 = 0.322054492873589;     i4 = 4;
            a5 =-0.275397277153155;     i5 = 5;
            a6 = 0.076379745812837;     i6 = 6;
            a7 =-0.004974217932391;     i7 = 7;
        }else if(i == 2){
            a0 =-0.270618371836400;     i0 = 0;
            a1 = 2.608899664176471;     i1 = 1;
            a2 =-6.229396164399805;     i2 = 2;
            a3 = 7.401795540292633;     i3 = 3;
            a4 =-6.019570786427897;     i4 = 4;
            a5 = 3.491391607443278;     i5 = 5;
            a6 =-1.143341576978401;     i6 = 6;
            a7 = 0.160840087730119;     i7 = 7;
        }else if(i == 3){
            a0 =-0.038775768114076;     i0 = 0;
            a1 = 0.095390654740871;     i1 = 1;
            a2 = 1.069538898422949;     i2 = 2;
            a3 =-2.504735762568496;     i3 = 3;
            a4 = 1.772892855219715;     i4 = 4;
            a5 =-0.602221665195144;     i5 = 5;
            a6 = 0.256099771243681;     i6 = 6;
            a7 =-0.048188983749501;     i7 = 7;
        }else if(i == 4){
            a0 =-0.161791607686480;     i0 = 0;
            a1 = 1.177738635455464;     i1 = 1;
            a2 =-3.743295966085387;     i2 = 2;
            a3 = 7.629773234288446;     i3 = 3;
            a4 =-8.978199716382386;     i4 = 4;
            a5 = 5.320463633854483;     i5 = 5;
            a6 =-1.442831784216072;     i6 = 6;
            a7 = 0.198143570771935;     i7 = 7;
        }else if(i == 5){
            a0 = 0.061315880970547;     i0 = 0;
            a1 =-0.329403900039681;     i1 = 1;
            a2 = 0.710128166312073;     i2 = 2;
            a3 =-0.847686078265718;     i3 = 3;
            a4 = 1.740199396764717;     i4 = 4;
            a5 =-2.591052963604934;     i5 = 5;
            a6 = 1.357538031860576;     i6 = 6;
            a7 =-0.109753113995674;     i7 = 7;
            a8 = 0.008714579998094;     i8 = 8;
        }else if(i == 6){
            a0 =-0.021840052913065;     i0 = 0;
            a1 = 0.125943642664851;     i1 = 1;
            a2 =-0.320584941290064;     i2 = 2;
            a3 = 0.496954711803904;     i3 = 3;
            a4 =-0.650570224322142;     i4 = 4;
            a5 = 1.871461736877637;     i5 = 5;
            a6 =-2.917045539749188;     i6 = 6;
            a7 = 1.565851441490449;     i7 = 7;
            a8 =-0.162184436527372;     i8 = 8;
            a9 = 0.012013661964991;     i9 = 9;
        }else if(i == 7){
            a0 = 0;                     i0 = 0;
            a1 =-0.007517959778232;     i1 = 1;
            a2 = 0.041336963439990;     i2 = 2;
            a3 =-0.085710205677551;     i3 = 3;
            a4 = 0.081890894364834;     i4 = 4;
            a5 =-0.138682909879055;     i5 = 5;
            a6 = 1.435250490561686;     i6 = 6;
            a7 =-2.675494883306515;     i7 = 7;
            a8 = 1.486573284792685;     i8 = 8;
            a9 =-0.148657328479269;     i9 = 9;
            a10 = 0.011011653961427;    i10 = 10;
        }else if(i == n1 - 8){
            a0 = 0.011011653961427;     i0 = n1-11;
            a1 =-0.148657328479269;     i1 = n1-10;
            a2 = 1.486573284792685;     i2 = n1-9;
            a3 =-2.675494883306515;     i3 = n1-8;
            a4 = 1.435250490561686;     i4 = n1-7;
            a5 =-0.138682909879055;     i5 = n1-6;
            a6 = 0.081890894364834;     i6 = n1-5;
            a7 =-0.085710205677551;     i7 = n1-4;
            a8 = 0.041336963439990;     i8 = n1-3;
            a9 =-0.007517959778232;     i9 = n1-2;
        }else if(i == n1 - 7){
            a0 = 0.012013661964991;     i0 = n1-10;
            a1 =-0.162184436527372;     i1 = n1-9;
            a2 = 1.565851441490449;     i2 = n1-8;
            a3 =-2.917045539749188;     i3 = n1-7;
            a4 = 1.871461736877637;     i4 = n1-6;
            a5 =-0.650570224322142;     i5 = n1-5;
            a6 = 0.496954711803904;     i6 = n1-4;
            a7 =-0.320584941290064;     i7 = n1-3;
            a8 = 0.125943642664851;     i8 = n1-2;
            a9 =-0.021840052913065;     i9 = n1-1;
        }else if(i == n1 - 6){
            a0 = 0.008714579998094;     i0 = n1-9;
            a1 =-0.109753113995674;     i1 = n1-8;
            a2 = 1.357538031860576;     i2 = n1-7;
            a3 =-2.591052963604934;     i3 = n1-6;
            a4 = 1.740199396764717;     i4 = n1-5;
            a5 =-0.847686078265718;     i5 = n1-4;
            a6 = 0.710128166312073;     i6 = n1-3;
            a7 =-0.329403900039681;     i7 = n1-2;
            a8 = 0.061315880970547;     i8 = n1-1;
        }else if(i == n1 - 5){
            a0 = 0.198143570771935;     i0 = n1-8;
            a1 =-1.442831784216072;     i1 = n1-7;
            a2 = 5.320463633854483;     i2 = n1-6;
            a3 =-8.978199716382386;     i3 = n1-5;
            a4 = 7.629773234288446;     i4 = n1-4;
            a5 =-3.743295966085387;     i5 = n1-3;
            a6 = 1.177738635455464;     i6 = n1-2;
            a7 =-0.161791607686480;     i7 = n1-1;
        }else if(i == n1 - 4){
            a0 =-0.048188983749501;     i0 = n1-8;
            a1 = 0.256099771243681;     i1 = n1-7;
            a2 =-0.602221665195144;     i2 = n1-6;
            a3 = 1.772892855219715;     i3 = n1-5;
            a4 =-2.504735762568496;     i4 = n1-4;
            a5 = 1.069538898422949;     i5 = n1-3;
            a6 = 0.095390654740871;     i6 = n1-2;
            a7 =-0.038775768114076;     i7 = n1-1;
        }else if(i == n1 - 3){
            a0 = 0.160840087730119;     i0 = n1-8;
            a1 =-1.143341576978401;     i1 = n1-7;
            a2 = 3.491391607443278;     i2 = n1-6;
            a3 =-6.019570786427897;     i3 = n1-5;
            a4 = 7.401795540292633;     i4 = n1-4;
            a5 =-6.229396164399805;     i5 = n1-3;
            a6 = 2.608899664176471;     i6 = n1-2;
            a7 =-0.270618371836400;     i7 = n1-1;
        }else if(i == n1 - 2){
            a0 =-0.004974217932391;     i0 = n1-8;
            a1 = 0.076379745812837;     i1 = n1-7;
            a2 =-0.275397277153155;     i2 = n1-6;
            a3 = 0.322054492873589;     i3 = n1-5;
            a4 = 0.112257511414158;     i4 = n1-4;
            a5 = 0.443635669435212;     i5 = n1-3;
            a6 =-1.558734408191758;     i6 = n1-2;
            a7 = 0.884778483741509;     i7 = n1-1;
        }else if(i == n1 - 1){
            a0 = 0.569218478941565;     i0 = n1-7;
            a1 =-1.606558237520481;     i1 = n1-6;
            a2 = 0.411180670145605;     i2 = n1-5;
            a3 = 2.036490115791102;     i3 = n1-4;
            a4 =-0.049249177165591;     i4 = n1-3;
            a5 =-3.038214359671524;     i5 = n1-2;
            a6 = 1.677132509479325;     i6 = n1-1;
        }else{
            a0 = 0.011111111111111;     i0 = i - 3;
            a1 =-0.15;                  i1 = i - 2;
            a2 = 1.5;                   i2 = i - 1;
            a3 =-2.722222222222222;     i3 = i;
            a4 = 1.5;                   i4 = i + 1;
            a5 =-0.15;                  i5 = i + 2;
            a6 = 0.011111111111111;     i6 = i + 3;
        }
    }
	return (1./(h*h))*( a0*MG(i0,j,k,c1,c2,c3) + a1*MG(i1,j,k,c1,c2,c3) + a2*MG(i2,j,k,c1,c2,c3) + a3*MG(i3,j,k,c1,c2,c3) + a4*MG(i4,j,k,c1,c2,c3) + a5*MG(i5,j,k,c1,c2,c3) + a6*MG(i6,j,k,c1,c2,c3) + a7*MG(i7,j,k,c1,c2,c3) + a8*MG(i8,j,k,c1,c2,c3) + a9*MG(i9,j,k,c1,c2,c3) + a10*MG(i10,j,k,c1,c2,c3) );
}

double myData::D1P(int i, int j, int k, const double h, const int order, const int n1, const int n2, const int n3, const int c1, const int c2, const int c3) const{
    int i0 = i, i1 = i, i2 = i, i3 = i, i4 = i, i5 = i;
    int i6 = i, i7 = i, i8 = i, i9 = i, i10 = i;
    double a0 = 0, a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0;
    double a6 = 0, a7 = 0, a8 = 0, a9 = 0, a10 = 0;
    
    if(order == 2){
        if(i == 0){
            a0 =-1;                     i0 = 0;
            a1 = 1;                     i1 = 1;
        }else if(i == n1-1){
            a0 =-1;                     i0 = n1-2;
            a1 = 1;                     i1 = n1-1;
        }else{
            a0 =-0.5;                   i0 = i - 1;
            a1 = 0;                     i1 = i;
            a2 = 0.5;                   i2 = i + 1;
        }
    }
    else if(order == 4){
        if (i == 0) {
            a0 =-1.576527260783884;     i0 = 0;
            a1 = 1.972775709802204;     i1 = 1;
            a2 = 0.040836435296694;     i2 = 2;
            a3 =-0.693890956864462;     i3 = 3;
            a4 = 0.256806072549449;     i4 = 4;
        }else if(i == 1){
            a0 =-0.451047142261160;     i0 = 0;
            a1 = 0;                     i1 = 1;
            a2 = 0.185148101793909;     i2 = 2;
            a3 = 0.455027117229871;     i3 = 3;
            a4 =-0.210262828535670;     i4 = 4;
            a5 = 0.021134751773050;     i5 = 5;
        }else if(i == 2){
            a0 =-0.020427163198248;     i0 = 0;
            a1 =-0.405074844833881;     i1 = 1;
            a2 = 0;                     i2 = 2;
            a3 = 0.134720700985761;     i3 = 3;
            a4 = 0.452993793355239;     i4 = 4;
            a5 =-0.162212486308872;     i5 = 5;
        }else if(i == 3){
            a0 = 0.179140757490107;     i0 = 0;
            a1 =-0.513802524967025;     i1 = 1;
            a2 =-0.069530808366309;     i2 = 2;
            a3 = 0;                     i3 = 3;
            a4 = 0.368098737516488;     i4 = 4;
            a5 = 0.036093838326738;     i5 = 5;
        }else if(i == 4){
            a0 =-0.088216121348878;     i0 = 0;
            a1 = 0.315908236179015;     i1 = 1;
            a2 =-0.311081860348502;     i2 = 2;
            a3 =-0.489783126488655;     i3 = 3;
            a4 = 0;                     i4 = 0;
            a5 = 0.663432368058167;     i5 = 5;
            a6 =-0.090259496051147;     i6 = 6;
        }else if(i == 5){
            a0 = 0;                     i0 = 0;
            a1 =-0.029028856952944;     i1 = 1;
            a2 = 0.101835934814000;     i2 = 2;
            a3 =-0.043904283847899;     i3 = 3;
            a4 =-0.606500263586147;     i4 = 4;
            a5 = 0;                     i5 = 0;
            a6 = 0.660111393797703;     i6 = 6;
            a7 =-0.082513924224713;     i7 = 7;
        }else if(i == n1-6){
            a0 = 0.082513924224713;     i0 = n1-8;
            a1 =-0.660111393797703;     i1 = n1-7;
            a2 = 0;                     i2 = n1-6;
            a3 = 0.606500263586147;     i3 = n1-5;
            a4 = 0.043904283847899;     i4 = n1-4;
            a5 =-0.101835934814000;     i5 = n1-3;
            a6 = 0.029028856952944;     i6 = n1-2;
        }else if(i == n1-5){
            a0 = 0.090259496051147;     i0 = n1-7;
            a1 =-0.663432368058167;     i1 = n1-6;
            a2 = 0;                     i2 = n1-5;
            a3 = 0.489783126488655;     i3 = n1-4;
            a4 = 0.311081860348502;     i4 = n1-3;
            a5 =-0.315908236179015;     i5 = n1-2;
            a6 = 0.088216121348878;     i6 = n1-1;
        }else if(i == n1-4){
            a0 =-0.036093838326738;     i0 = n1-6;
            a1 =-0.368098737516488;     i1 = n1-5;
            a2 = 0;                     i2 = n1-4;
            a3 = 0.069530808366309;     i3 = n1-3;
            a4 = 0.513802524967025;     i4 = n1-2;
            a5 =-0.179140757490107;     i5 = n1-1;
        }else if(i == n1-3){
            a0 = 0.162212486308872;     i0 = n1-6;
            a1 =-0.452993793355239;     i1 = n1-5;
            a2 =-0.134720700985761;     i2 = n1-4;
            a3 = 0;                     i3 = n1-3;
            a4 = 0.405074844833881;     i4 = n1-2;
            a5 = 0.020427163198248;     i5 = n1-1;
        }else if(i == n1-2){
            a0 =-0.021134751773050;     i0 = n1-6;
            a1 = 0.210262828535670;     i1 = n1-5;
            a2 =-0.455027117229871;     i2 = n1-4;
            a3 =-0.185148101793909;     i3 = n1-3;
            a4 = 0;                     i4 = n1-2;
            a5 = 0.451047142261160;     i5 = n1-1;
        }else if(i == n1-1){
            a0 =-0.256806072549449;     i0 = n1-5;
            a1 = 0.693890956864462;     i1 = n1-4;
            a2 =-0.040836435296694;     i2 = n1-3;
            a3 =-1.972775709802204;     i3 = n1-2;
            a4 = 1.576527260783884;     i4 = n1-1;
        }else{
            a0 = 0.083333333333333;     i0 = i - 2;
            a1 =-0.666666666666667;     i1 = i - 1;
            a2 = 0;                     i2 = i;
            a3 = 0.666666666666667;     i3 = i + 1;
            a4 =-0.083333333333333;     i4 = i + 2;
        }
    }
    else{
        if(i == 0){
            a0 =-1.694834962162858;     i0 = 0;
            a1 = 2.245634824947698;     i1 = 1;
            a2 =-0.055649692295628;     i2 = 2;
            a3 =-0.670383570370653;     i3 = 3;
            a4 =-0.188774952148393;     i4 = 4;
            a5 = 0.552135032829910;     i5 = 5;
            a6 =-0.188126680800077;     i6 = 6;
        }else if(i == 1){
            a0 =-0.434411786832708;     i0 = 0;
            a1 = 0;                     i1 = 1;
            a2 = 0.107043134706685;     i2 = 2;
            a3 = 0.420172642668695;     i3 = 3;
            a4 = 0.119957288069806;     i4 = 4;
            a5 =-0.328691543801578;     i5 = 5;
            a6 = 0.122487487014485;     i6 = 6;
            a7 =-0.006557221825386;     i7 = 7;
        }else if(i == 2){
            a0 = 0.063307644169533;     i0 = 0;
            a1 =-0.629491308812471;     i1 = 1;
            a2 = 0;                     i2 = 2;
            a3 = 0.809935419586724;     i3 = 3;
            a4 =-0.699016381364484;     i4 = 4;
            a5 = 0.850345731199969;     i5 = 5;
            a6 =-0.509589652965290;     i6 = 6;
            a7 = 0.114508548186019;     i7 = 7;
        }else if(i == 3){
            a0 = 0.110198643174386;     i0 = 0;
            a1 =-0.357041083340051;     i1 = 1;
            a2 =-0.117033418681039;     i2 = 2;
            a3 = 0;                     i3 = 3;
            a4 = 0.120870009174558;     i4 = 4;
            a5 = 0.349168902725368;     i5 = 5;
            a6 =-0.104924741749615;     i6 = 6;
            a7 =-0.001238311303608;     i7 = 7;
        }else if(i == 4){
            a0 = 0.133544619364965;     i0 = 0;
            a1 =-0.438678347579289;     i1 = 1;
            a2 = 0.434686341173840;     i2 = 2;
            a3 =-0.520172867814934;     i3 = 3;
            a4 = 0;                     i4 = 4;
            a5 = 0.049912002176267;     i5 = 5;
            a6 = 0.504693510958978;     i6 = 6;
            a7 =-0.163985258279827;     i7 = 7;
        }else if(i == 5){
            a0 =-0.127754693486067;     i0 = 0;
            a1 = 0.393149407857401;     i1 = 1;
            a2 =-0.172955234680916;     i2 = 2;
            a3 =-0.491489487857764;     i3 = 3;
            a4 =-0.016325050231672;     i4 = 4;
            a5 = 0;                     i5 = 5;
            a6 = 0.428167552785852;     i6 = 6;
            a7 =-0.025864364383975;     i7 = 7;
            a8 = 0.013071869997141;     i8 = 8;
        }else if(i == 6){
            a0 = 0.060008241515128;     i0 = 0;
            a1 =-0.201971348965594;     i1 = 1;
            a2 = 0.142885356631256;     i2 = 2;
            a3 = 0.203603636754774;     i3 = 3;
            a4 =-0.227565385120003;     i4 = 4;
            a5 =-0.590259111130048;     i5 = 5;
            a6 = 0;                     i6 = 6;
            a7 = 0.757462553894374;     i7 = 7;
            a8 =-0.162184436527372;     i8 = 8;
            a9 = 0.018020492947486;     i9 = 9;
        }else if(i == 7){
            a0 = 0;                     i0 = 0;
            a1 = 0.009910488565285;     i1 = 1;
            a2 =-0.029429452176588;     i2 = 2;
            a3 = 0.002202493355677;     i3 = 3;
            a4 = 0.067773581604826;     i4 = 4;
            a5 = 0.032681945726690;     i5 = 5;
            a6 =-0.694285851935105;     i6 = 6;
            a7 = 0;                     i7 = 7;
            a8 = 0.743286642396343;     i8 = 8;
            a9 =-0.148657328479269;     i9 = 9;
            a10 =0.016517480942141;     i10 = 10;
        }else if(i == n1-8){
            a0 =-0.016517480942141;     i0 = n1-11;
            a1 = 0.148657328479269;     i1 = n1-10;
            a2 =-0.743286642396343;     i2 = n1-9;
            a3 = 0;                     i3 = n1-8;
            a4 = 0.694285851935105;     i4 = n1-7;
            a5 =-0.032681945726690;     i5 = n1-6;
            a6 =-0.067773581604826;     i6 = n1-5;
            a7 =-0.002202493355677;     i7 = n1-4;
            a8 = 0.029429452176588;     i8 = n1-3;
            a9 =-0.009910488565285;     i9 = n1-2;
        }else if(i == n1-7){
            a0 =-0.018020492947486;     i0 = n1-10;
            a1 = 0.162184436527372;     i1 = n1-9;
            a2 =-0.757462553894374;     i2 = n1-8;
            a3 = 0;                     i3 = n1-7;
            a4 = 0.590259111130048;     i4 = n1-6;
            a5 = 0.227565385120003;     i5 = n1-5;
            a6 =-0.203603636754774;     i6 = n1-4;
            a7 =-0.142885356631256;     i7 = n1-3;
            a8 = 0.201971348965594;     i8 = n1-2;
            a9 =-0.060008241515128;     i9 = n1-1;
        }else if(i == n1-6){
            a0 =-0.013071869997141;     i0 = n1-9;
            a1 = 0.025864364383975;     i1 = n1-8;
            a2 =-0.428167552785852;     i2 = n1-7;
            a3 = 0;                     i3 = n1-6;
            a4 = 0.016325050231672;     i4 = n1-5;
            a5 = 0.491489487857764;     i5 = n1-4;
            a6 = 0.172955234680916;     i6 = n1-3;
            a7 =-0.393149407857401;     i7 = n1-2;
            a8 = 0.127754693486067;     i8 = n1-1;
        }else if(i == n1-5){
            a0 = 0.163985258279827;     i0 = n1-8;
            a1 =-0.504693510958978;     i1 = n1-7;
            a2 =-0.049912002176267;     i2 = n1-6;
            a3 = 0;                     i3 = n1-5;
            a4 = 0.520172867814934;     i4 = n1-4;
            a5 =-0.434686341173840;     i5 = n1-3;
            a6 = 0.438678347579289;     i6 = n1-2;
            a7 =-0.133544619364965;     i7 = n1-1;
        }else if(i == n1-4){
            a0 = 0.001238311303608;     i0 = n1-8;
            a1 = 0.104924741749615;     i1 = n1-7;
            a2 =-0.349168902725368;     i2 = n1-6;
            a3 =-0.120870009174558;     i3 = n1-5;
            a4 = 0;                     i4 = n1-4;
            a5 = 0.117033418681039;     i5 = n1-3;
            a6 = 0.357041083340051;     i6 = n1-2;
            a7 =-0.110198643174386;     i7 = n1-1;
        }else if(i == n1-3){
            a0 =-0.114508548186019;     i0 = n1-8;
            a1 = 0.509589652965290;     i1 = n1-7;
            a2 =-0.850345731199969;     i2 = n1-6;
            a3 = 0.699016381364484;     i3 = n1-5;
            a4 =-0.809935419586724;     i4 = n1-4;
            a5 = 0;                     i5 = n1-3;
            a6 = 0.629491308812471;     i6 = n1-2;
            a7 =-0.063307644169533;     i7 = n1-1;
        }else if(i == n1-2){
            a0 = 0.006557221825386;     i0 = n1-8;
            a1 =-0.122487487014485;     i1 = n1-7;
            a2 = 0.328691543801578;     i2 = n1-6;
            a3 =-0.119957288069806;     i3 = n1-5;
            a4 =-0.420172642668695;     i4 = n1-4;
            a5 =-0.107043134706685;     i5 = n1-3;
            a6 = 0;                     i6 = n1-2;
            a7 = 0.434411786832708;     i7 = n1-1;
        }else if(i == n1-1){
            a0 = 0.188126680800077;     i0 = n1-7;
            a1 =-0.552135032829910;     i1 = n1-6;
            a2 = 0.188774952148393;     i2 = n1-5;
            a3 = 0.670383570370653;     i3 = n1-4;
            a4 = 0.055649692295628;     i4 = n1-3;
            a5 =-2.245634824947698;     i5 = n1-2;
            a6 = 1.694834962162858;     i6 = n1-1;
        }else{
            a0 =-0.016666666666667;     i0 = i - 3;
            a1 = 0.15;                  i1 = i - 2;
            a2 =-0.75;                  i2 = i - 1;
            a3 = 0;                     i3 = i;
            a4 = 0.75;                  i4 = i + 1;
            a5 =-0.15;                  i5 = i + 2;
            a6 = 0.016666666666667;     i6 = i + 3;
        }
    }
    
    return (1./h)*( a0*MG(i0,j,k,c1,c2,c3) + a1*MG(i1,j,k,c1,c2,c3) + a2*MG(i2,j,k,c1,c2,c3) + a3*MG(i3,j,k,c1,c2,c3) + a4*MG(i4,j,k,c1,c2,c3) + a5*MG(i5,j,k,c1,c2,c3) + a6*MG(i6,j,k,c1,c2,c3) + a7*MG(i7,j,k,c1,c2,c3) + a8*MG(i8,j,k,c1,c2,c3) + a9*MG(i9,j,k,c1,c2,c3) + a10*MG(i10,j,k,c1,c2,c3) );
}

