//
//  main.cpp
//  PDE
//
//  Created by Ali Dorostkar on 3/23/12.
//  Copyright (c) 2012 Uppsala. All rights reserved.
//

#include <iostream>
#include <cstdio>
#include <cstring>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <sstream>
#include "solver.h"
#include "myData.h"


#define idx(i,j,k) (i) + (j)*LN_X + (k)*LN_X*LN_Y

void RikerWavelet(Vdata &, Vdata &, FDM::baseSolver &, const Vdata&, int, double, void *);
void finale(vector<myData> & u, FDM::baseSolver & sol, double tEnd, void * args);
void start(vector<myData> & u, FDM::baseSolver & sol, void * args);
void MMS(vector<myData> & F_out, vector<myData> & Ftt_out, FDM::baseSolver& sol, const vector<myData>& u,int step,double t, void * args);
void plot(const Vdata &u, FDM::baseSolver& sol,int k, double t, void *args);

int main(int argc, const char * argv[])
{
//    cout << "Setting program parameters..."<<endl;
//
//    int order = 4;
//    int nx[4] = {51, 101, 201, 401};
//    double err1[4];
//    int ny = 1;
//    int nz = 1;
//    for (int n = 0; n < 4; n++) {
//        double h = 2*M_PI/(nx[n]-1);
//        myData u(nx[n],ny,nz),ux(nx[n],ny,nz),uxx(nx[n],ny,nz);
//        
//        for (int i = 0; i < u.Nx(); i++) {
//            double a = sin(i*h);
//            double b = cos(i*h);
//            for (int j = 0; j < u.Ny(); j++) {
//                for (int k = 0; k < u.Nz(); k++) {
//                    u(i,j,k) = a;
//                    ux(i,j,k) = b;
//                    uxx(i,j,k) = -a;
//                }
//            }
//        }
//
//        
//        
//        myData vx = u.D1(h, order, X_DIR);
////        myData vxx = u.D2(h, order, X_DIR);
//    
//        myData w = vx - ux;
////        w = vxx - uxx;
//        err1[n] = w.Norm()*sqrt(h);
////        err1[n] = w.Norm()*h;
//        
//        for (int i = 0; i < u.Nx(); i++)
//            for (int j = 0; j < u.Ny(); j++)
//                for (int k = 0; k < u.Nz(); k++)
//                    cout<<w(i,j,k)<<endl;
//    }
//
//    for (int i = 0; i < 4; i++)
//        cout<<"Error :"<<err1[i]<<endl;
//    
//    for (int i = 0; i < 4; i++)
//        cout<<"Value :"<<log2(err1[i])<<endl;
//    
//    for (int i = 0; i < 3; i++)
//        cout<<"D1 :"<<log2(err1[i]) - log2(err1[i+1])<<endl;
    

    double h = 0.04;
    int order = 2;
    double tEnd = 2;
    double dt_coef = 0.1;
    
    for (int i = 0; i < argc; i++) {
        if (strcmp(argv[i],"-o") == 0)
            order = atoi(argv[i+1]);
        if (strcmp(argv[i],"-h") == 0)
            h = atof(argv[i+1]);
        if (strcmp(argv[i],"-dt") == 0)
            dt_coef = atof(argv[i+1]);
        if (strcmp(argv[i],"-t") == 0)
            tEnd = atof(argv[i+1]);
        if (strcmp(argv[i],"--help") == 0){
            cout<<"Useage : ./fdm [options]"<<endl;
            cout<<"\t -o : order of accuracy (Default is 2)"<<endl;
            cout<<"\t -h : spatial step size (Default is 0.04)"<<endl;
            cout<<"\t -dt : time step size coefficient(Default is c = 0.1 in dt = c*h)"<<endl;
            cout<<"\t -t : Final calculation time (Default is 2)"<<endl;
            return 0;
        }
    }
    double dt = dt_coef*h;
    
    cout<<"Initialising..."<<endl;
    cout<<"Order of accuracy: "<<order<<endl;
    cout<<"Spatial step size: "<<h<<endl;
    cout<<"Time step: "<<dt<<endl;
    cout<<"Final time: "<<tEnd<<endl;
    
    FDM::pmlSolver pmlSol(h,dt,order);
    FDM::solver refSol(h,dt,order);
    
    pmlSol.Set_EndTime(tEnd);
    refSol.Set_EndTime(tEnd);
    

    int LN_X = 2;
    int LN_Y = 1;
    int LN_Z = 1;
    
    int ref_LN_Y = 2;
    
    vector<double> x_layers(LN_X + 1);
    x_layers[0] = 0;
    x_layers[1] = 0.5;
    x_layers[2] = 1;    

    vector<double> y_layers(ref_LN_Y + 1);
    y_layers[0] = 0;
    y_layers[1] = 1;
    y_layers[2] = 6;
    
    vector<double> z_layers(LN_Z + 1);
    z_layers[0] = 0;
    z_layers[1] = 1;
    
    vector<double> pml_y_layers(LN_Y + 1);
    pml_y_layers[0] = 0;
    pml_y_layers[1] = 1;
    
    vector<double> ref_spd(LN_X*ref_LN_Y*LN_Z);
    ref_spd[0] = 1;
    ref_spd[1] = 1;
    ref_spd[2] = 1;
    ref_spd[3] = 1;
    
    vector<double> spd(LN_X*LN_Y*LN_Z);
    spd[0] = 1;
    spd[1] = 1;

    refSol.Set_Domain(x_layers, y_layers, z_layers, ref_spd);
    pmlSol.Set_Domain(x_layers, pml_y_layers, z_layers, spd);
    
    refSol.Set_FourceFunction(RikerWavelet, NULL);
    pmlSol.Set_FourceFunction(RikerWavelet, NULL);

    pmlSol.Set_Plotter(plot, NULL);
//////////
//    pmlSol.Set_Initializer(start, NULL);
//    pmlSol.Set_Finilizer(finale, NULL);
//    pmlSol.Set_FourceFunction(MMS, NULL);
//    pmlSol.Set_DampingCoeff(false,0);
//
//    refSol.Set_Initializer(start, NULL);
//    refSol.Set_Finilizer(finale, NULL);
//    refSol.Set_FourceFunction(MMS, NULL);
    
    Vdata pml = pmlSol.run();
    Vdata ref = refSol.run();
    double err = 0;
    double val = 0;

    for (int i = 0; i < LN_X; i++) {
        for (int j = 0; j < ref_LN_Y - 1; j++) {
            for (int k = 0; k < LN_Z; k++) {
                int ind = idx(i, j, k);
                val = (pml[ind]-ref[ind]).NormSquare();
                err+=val;
            }
        }
    }
    
    err = sqrt(err)*sqrt(h*h*h);
    cout << "-----------------------"<<endl;
    cout<<"Error:"<<err<<endl;
    cout<<"log2(error):"<<log2(err)<<endl;
    cout << "-----------------------"<<endl;
    
    return 0;
}

void RikerWavelet(Vdata & F_out, Vdata & Ftt_out, FDM::baseSolver& sol, const Vdata& u, int step, double t, void * args){
    int layerNum_x = sol.X_LayerNumber();
    int layerNum_y = sol.Y_LayerNumber();
    int layerNum_z = sol.Z_LayerNumber();
    double f0 = 10;
    double a = 2*(M_PI*f0)*(M_PI*f0);
    double t0 = 1/f0;
    double g = (2*a*(t - t0))*(2*a*(t - t0)) - 2*a;
    double h = exp(-a*(t - t0)*(t - t0));
    double F = h*g;
    double x_0 = 0.5;
    double y_0 = 0.5;
    double z_0 = 0.5;
    double delta = 0.1;
    
    
    double hg__ = 8*a*a*h;
    double h__g = F*g;
    double h_g_ = -16*a*a*a*(t-t0)*(t-t0)*h;
    
    if (h == 0)
        return;
    
    double gx = 0;
    double x,y,z;
    int Nx,Ny,Nz;
    int index = 0;
    for (int ln_z = 0; ln_z < layerNum_z; ln_z++) {
        for (int ln_y = 0; ln_y < layerNum_y; ln_y++) {
            for (int ln = 0; ln < layerNum_x; ln++) {
                
                index = ln + ln_y*layerNum_x + ln_z*layerNum_x*layerNum_y;
                Nx = u[index].Nx();
                Ny = u[index].Ny();
                Nz = u[index].Nz();
#pragma omp parallel for private(x,y,z)
                for (int k = 0; k < Nz; k++) {
                    z = sol.Z(k,ln_z)-z_0;
                    for (int j = 0; j < Ny; j++) {
                        y = sol.Y(j,ln_y)-y_0;
                        for (int i = 0; i < Nx; i++) {
                            x = sol.X(i,ln)-x_0;
                            gx = exp(-(x*x + y*y + z*z)/(delta*delta));
                            F_out[index](i,j,k) = F*gx;
                            Ftt_out[index](i,j,k) = (hg__ + h__g + 2*h_g_)*gx;
                        }
                    }
                }
            }
        }
    }
}

void plot(const Vdata &u, FDM::baseSolver& sol,int k, double t, void *args){
    int layerNum_x = sol.X_LayerNumber();
    int layerNum_y = sol.Y_LayerNumber();
    int layerNum_z = sol.Z_LayerNumber();

    int Nx,Ny,Nz;
    double x,y,z;
    int index = 0;
    
    ofstream myfile;
    string s = "../../data";
    stringstream ss;
    ss << k;
    s += ss.str();
    s+= ".dat";
    myfile.open (s.data());

    
    for (int ln_z = 0; ln_z < layerNum_z; ln_z++) {
        for (int ln_y = 0; ln_y < layerNum_y; ln_y++) {
            for (int ln = 0; ln < layerNum_x; ln++) {
                
                index = ln + ln_y*layerNum_x + ln_z*layerNum_x*layerNum_y;
                Nx = u[index].Nx();
                Ny = u[index].Ny();
                Nz = u[index].Nz();

#pragma omp parallel for private(x,y,z)
                for (int k = 0; k < Nz; k++) {
                    z = sol.Z(k,ln_z);
                    for (int j = 0; j < Ny; j++) {
                        y = sol.Y(j,ln_y);
                        for (int i = 0; i < Nx; i++) {
                            x = sol.X(i,ln);
                            myfile <<x<<" "<<y<<" "<<z<<" "<<setprecision(12)<<u[index](i,j,k)<<"\n";
                            
                        }
                    }
                }
            }
        }
    }
    myfile.close();
}

void MMS(vector<myData> & F_out, vector<myData> & Ftt_out, FDM::baseSolver& sol, const vector<myData>& u,int step,double t, void * args){
    double x = 0;
    double y = 0;
    double z = 0;
    
    int LN_X = sol.X_LayerNumber();
    int LN_Y = sol.Y_LayerNumber();
    int LN_Z = sol.Z_LayerNumber();
    
    int Nx,Ny,Nz;
    int index = 0;
    for (int lz = 0; lz < LN_Z; lz++) {
        for (int ly = 0; ly < LN_Y; ly++) {
            for (int lx = 0; lx < LN_X; lx++) {
                
                index = idx(lx, ly, lz);
                Nx = u[index].Nx();
                Ny = u[index].Ny();
                Nz = u[index].Nz();
                
                myData ux(Nx,Ny,Nz);
                myData uy(Nx,Ny,Nz);
                myData uz(Nx,Ny,Nz);
                
                myData Bux(Nx,Ny,Nz);
                myData Buy(Nx,Ny,Nz);
                myData Buz(Nx,Ny,Nz);
                
                for (int k = 0; k < Nz; k++) {
                    z = sol.Z(k,lz);
                    for (int j = 0; j < Ny; j++) {
                        y = sol.Y(j,ly);
                        for (int i = 0; i < Nx; i++) {
                            x = sol.X(i,lx);
                            
                            
                            F_out[index](i,j,k)     =  8*M_PI*M_PI*cos(2*M_PI*t)*cos(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z);
                            Ftt_out[index](i,j,k)   = -4*M_PI*M_PI*F_out[index](i,j,k);
                            
                            ux(i,j,k)               = -2*M_PI*cos(2*M_PI*t)*sin(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z);//ux boundary
                            uy(i,j,k)               = -2*M_PI*cos(2*M_PI*t)*cos(2*M_PI*x)*sin(2*M_PI*y)*cos(2*M_PI*z);//uy
                            uz(i,j,k)               = -2*M_PI*cos(2*M_PI*t)*cos(2*M_PI*x)*cos(2*M_PI*y)*sin(2*M_PI*z);//uz
                        }
                    }
                }
                
                //Bux setup
                if (lx == 0) {
                    for (int k = 0; k < Nz; k++) {
                        for (int j = 0; j < Ny; j++) {
                            Bux(0,j,k) = -1;
                        }
                    }
                }
                
                if (lx == LN_X - 1) {
                    for (int k = 0; k < Nz; k++) {
                        for (int j = 0; j < Ny; j++) {
                            Bux(Nx - 1,j,k) = 1;
                        }
                    }
                }
                //
                
                ///Buy setup
                if (ly == 0) {
                    for (int k = 0; k < Nz; k++) {
                        for (int i = 0; i < Nx; i++) {
                            Buy(i,0,k) = -1;
                        }
                    }
                }
                
                if (ly == LN_Y - 1) {
                    for (int k = 0; k < Nz; k++) {
                        for (int i = 0; i < Nx; i++) {
                            Buy(i,Ny - 1,k) = 1;
                        }
                    }
                }
                //
                
                //Buz setup
                if (lz == 0) {
                    for (int j = 0; j < Ny; j++) {
                        for (int i = 0; i < Nx; i++) {
                            Buz(i,j,0) = -1;
                        }
                    }
                }
                
                if (lz == LN_Z - 1) {
                    for (int j = 0; j < Ny; j++) {
                        for (int i = 0; i < Nx; i++) {
                            Buz(i,j,Nz - 1) = 1;
                        }
                    }
                }
                //
                
                F_out[index] += sol.H_INV()*(Bux*ux + Buy*uy + Buz*uz);
            }
        }
    }
}

void start(vector<myData> & u, FDM::baseSolver & sol, void * args){
    double x = 0;
    double y = 0;
    double z = 0;
    
    int LN_X = sol.X_LayerNumber();
    int LN_Y = sol.Y_LayerNumber();
    int LN_Z = sol.Z_LayerNumber();
    
    int Nx,Ny,Nz;
    int index = 0;
    for (int lz = 0; lz < LN_Z; lz++) {
        for (int ly = 0; ly < LN_Y; ly++) {
            for (int lx = 0; lx < LN_X; lx++) {
                
                index = idx(lx, ly, lz);
                Nx = u[index].Nx();
                Ny = u[index].Ny();
                Nz = u[index].Nz();
                
                for (int k = 0; k < Nz; k++) {
                    z = sol.Z(k,lz);
                    for (int j = 0; j < Ny; j++) {
                        y = sol.Y(j,ly);
                        for (int i = 0; i < Nx; i++) {
                            x = sol.X(i,lx);
                            
                            u[index](i,j,k) = cos(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z);
                        }
                    }
                }
            }
        }
    }
}
void finale(vector<myData> & u, FDM::baseSolver & sol, double tEnd, void * args){
    double error = 0;
    
    double h = sol.dx();
    double dt = sol.Dt();
    double x = 0;
    double y = 0;
    double z = 0;
    double ex = 0;
    
    int LN_X = sol.X_LayerNumber();
    int LN_Y = sol.Y_LayerNumber();
    int LN_Z = sol.Z_LayerNumber();
    
    int Nx,Ny,Nz;
    int index = 0;
    for (int lz = 0; lz < LN_Z; lz++) {
        for (int ly = 0; ly < LN_Y; ly++) {
            for (int lx = 0; lx < LN_X; lx++) {
                
                index = idx(lx, ly, lz);
                Nx = u[index].Nx();
                Ny = u[index].Ny();
                Nz = u[index].Nz();
                for (int k = 0; k < Nz; k++) {
                    z = sol.Z(k,lz);
                    for (int j = 0; j < Ny; j++) {
                        y = sol.Y(j,ly);
                        for (int i = 0; i < Nx; i++) {
                            x = sol.X(i,lx);
                            
                            ex = cos(2*M_PI*tEnd)*cos(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z) - u[index](i,j,k);
                            error += ex*ex;
                        }
                    }
                }
            }
        }
    }
    
    error = sqrt(error)*sqrt(h*h*h);
//    std::system("clear");
//    cout << "\033[2J\033[1;1H";
    cout<<"--------------------------------------------------------------------"<<endl;
    cout
    << setw(9) << left << "Order"
    << setw(9) << left << "h"
    << setw(9) << left << "dt"
    << setw(15) << left << "final time"
    << setw(15) << left << "Error"
    << setw(15) << left << "log2(error)"
    << endl;
    cout<<"____________________________________________________________________"<<endl;
    cout
    << setw(9) << left << sol.Order()
    << setw(9) << left << h
    << setw(9) << left << dt
    << setw(15) << left << sol.EndTime()
    << setw(15) << left << error
    << setw(15) << left << log(error)/log(2)
    << endl;
    cout<<"===================================================================="<<endl;
}