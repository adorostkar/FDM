//
//  solver.cpp
//  PDE
//
//  Created by Ali Dorostkar on 3/24/12.
//  Copyright (c) 2012 Uppsala. All rights reserved.
//
//  SOLVER class solves 2D PDE
#include "solver.h"

// Compute index of the block with respect to placement
#define idx(i,j,k) (i) + (j)*LN_X + (k)*LN_X*LN_Y

// Compute index of the PML block with respect to placement
#define PML_idx(i,k) (i) + (k)*LN_X


/* ----======== Base class ========----- */
FDM::baseSolver::baseSolver(double _h, double _dt, int _order){
    h = _h;
    dt = _dt;
    order = _order;
    tEnd = 2;
    
    X_layers = vector<double>(2);
    Y_layers = vector<double>(2);
    Z_layers = vector<double>(2);
    speeds   = vector<double>(1);
    
    X_layers[0] = 0;
    X_layers[1] = 1;
    Y_layers[0] = 0;
    Y_layers[1] = 1;
    Z_layers[0] = 0;
    Z_layers[1] = 1;
    
    speeds[0] = 1;
    
    LN_X = 1;
    LN_Y = 1;
    LN_Z = 1;
    
    totalBlocks = LN_X*LN_Y*LN_Z;
    
    Ta_0 = 5;
    SPP = 10;
    
    init_args = NULL;
    initializer = NULL;
    fource_args = NULL;
    Fource = NULL;
    plot_args = NULL;
    plotter = NULL;
    finilizer = NULL;
    finilizer_args = NULL;
    
    Set_H_INV();
}

void FDM::baseSolver::Set_FourceFunction(void (*SF)(Vdata &, Vdata &, baseSolver &, const Vdata&, int, double, void *), void * args){
    if(SF != NULL){
        Fource = SF;
        fource_args = args;
    }
}

void FDM::baseSolver::Set_Initializer(void (*SI)(Vdata &, baseSolver &, void *), void * args){
    if(SI != NULL){
        initializer = SI;
        init_args = args;
    }
}

void FDM::baseSolver::Set_Plotter(void (*SP)(const Vdata &, baseSolver&,int, double , void *),void * args){
    plotter = SP;
    plot_args = args;
}

void FDM::baseSolver::Set_Finilizer(void (*FI)(Vdata &, baseSolver &, double, void *), void * args){
    finilizer = FI;
    finilizer_args = args;
}

void FDM::baseSolver::X_Interfc_L(const myData & u1, const myData & u1x,const double c1 ,const myData & u2, const myData & u2x,const double c2,myData & I2){
    int Nx1 = u1.Nx();
    int Ny = u2.Ny();
    int Nz = u2.Nz();
    double Ta_n = 0.5;
    
    myData I_u2_0(1,Ny,Nz);
    myData I_u2x_0(1,Ny,Nz);
    
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            //            I_u2_0(0,j,k) = Ta_0/h*(c1+c2)*(u2(0,j,k) - u1(Nx1-1,j,k));
            I_u2_0(0,j,k) = Ta_0/h*(c1+c2)*(u2(0,j,k) - u1(Nx1-1,j,k));
            I_u2x_0(0,j,k) = -Ta_n*(c2*u2x(0,j,k) - c1*u1x(Nx1-1,j,k));
        }
    }
    
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            I2(0,j,k) = I_u2_0(0,j,k) + I_u2x_0(0,j,k);
        }
    }
}

void FDM::baseSolver::X_Interfc_R(const myData & u1, const myData & u1x,const double c1,myData & I1 ,const myData & u2, const myData & u2x,const double c2){
    int Nx1 = u1.Nx();
    int Ny = u1.Ny();
    int Nz = u1.Nz();
    double Ta_n = 0.5;
    
    myData I_u1_n(1,Ny,Nz);
    myData I_u1x_n(1,Ny,Nz);
    
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            //            I_u1_n(0,j,k) = Ta_0/h*(c1+c2)*(u1(Nx1-1,j,k) - u2(0,j,k));
            I_u1_n(0,j,k) = Ta_0/h*(c1+c2)*(u1(Nx1-1,j,k) - u2(0,j,k));
            I_u1x_n(0,j,k) =  Ta_n*(c1*u1x(Nx1-1,j,k) - c2*u2x(0,j,k));
        }
    }
    
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            I1(Nx1-1,j,k) = I_u1_n(0,j,k) + I_u1x_n(0,j,k);
        }
    }
}

void FDM::baseSolver::Y_Interfc_L(const myData & u1, const myData & u1y,const double c1 ,const myData & u2, const myData & u2y,const double c2,myData & I2){
    int Nx = u1.Nx();
    int Ny1 = u1.Ny();
    int Nz = u1.Nz();
    double Ta_n = 0.5;
    
    myData I_u2_0(Nx,1,Nz);
    myData I_u2y_0(Nx,1,Nz);
    
    for (int k = 0; k < Nz; k++) {
        for (int i = 0; i < Nx; i++) {
            //            I_u2_0(i,0,k) = Ta_0/h*(c1+c2)*(u2(i,0,k) - u1(i,Ny1-1,k));
            I_u2_0(i,0,k) = Ta_0/h*(c1+c2)*(u2(i,0,k) - u1(i,Ny1-1,k));
            I_u2y_0(i,0,k) = -Ta_n*(c2*u2y(i,0,k) - c1*u1y(i,Ny1-1,k));
        }
    }
    
    for (int k = 0; k < Nz; k++) {
        for (int i = 0; i < Nx; i++) {
            I2(i,0,k) = I_u2_0(i,0,k) + I_u2y_0(i,0,k);
        }
    }
}

void FDM::baseSolver::Y_Interfc_R(const myData & u1, const myData & u1y,const double c1,myData & I1 ,const myData & u2, const myData & u2y,const double c2){
    int Nx = u1.Nx();
    int Ny1 = u1.Ny();
    int Nz = u1.Nz();
    double Ta_n = 0.5;
    
    myData I_u1_n(Nx,1,Nz);
    myData I_u1y_n(Nx,1,Nz);
    
    for (int k = 0; k < Nz; k++) {
        for (int i = 0; i < Nx; i++) {
            //            I_u1_n(i,0,k) = Ta_0/h*(c1+c2)*(u1(i,Ny1-1,k) - u2(i,0,k));
            I_u1_n(i,0,k) = Ta_0/h*(c1+c2)*(u1(i,Ny1-1,k) - u2(i,0,k));
            I_u1y_n(i,0,k) =  Ta_n*(c1*u1y(i,Ny1-1,k) - c2*u2y(i,0,k));
        }
    }
    
    for (int k = 0; k < Nz; k++) {
        for (int i = 0; i < Nx; i++) {
            I1(i,Ny1-1,k) = I_u1_n(i,0,k) + I_u1y_n(i,0,k);
        }
    }
}

void FDM::baseSolver::Z_Interfc_L(const myData & u1, const myData & u1z,const double c1 ,const myData & u2, const myData & u2z,const double c2,myData & I2){
    int Nx = u1.Nx();
    int Ny = u1.Ny();
    int Nz1 = u1.Nz();
    double Ta_n = 0.5;
    
    myData I_u2_0(Nx,Ny,1);
    myData I_u2z_0(Nx,Ny,1);
    
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            //            I_u2_0(i,j,0) = Ta_0/h*(c1+c2)*(u2(i,j,0) - u1(i,j,Nz1-1));
            I_u2_0(i,j,0) = Ta_0/h*(c1+c2)*(u2(i,j,0) - u1(i,j,Nz1-1));
            I_u2z_0(i,j,0) = -Ta_n*(c2*u2z(i,j,0) - c1*u1z(i,j,Nz1-1));
        }
    }
    
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i<Nx; i++) {
            I2(i,j,0) = I_u2_0(i,j,0) + I_u2z_0(i,j,0);
        }
    }
}

void FDM::baseSolver::Z_Interfc_R(const myData & u1, const myData & u1z,const double c1,myData & I1 ,const myData & u2, const myData & u2z,const double c2){
    int Nx = u1.Nx();
    int Ny = u1.Ny();
    int Nz1 = u1.Nz();
    double Ta_n = 0.5;
    
    myData I_u1_n(Nx,Ny,1);
    myData I_u1z_n(Nx,Ny,1);
    
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            //            I_u1_n(i,j,0) = Ta_0/h*(c1+c2)*(u1(i,j,Nz1-1) - u2(i,j,0));
            I_u1_n(i,j,0) = Ta_0/h*(c1+c2)*(u1(i,j,Nz1-1) - u2(i,j,0));
            I_u1z_n(i,j,0) =  Ta_n*(c1*u1z(i,j,Nz1-1) - c2*u2z(i,j,0));
        }
    }
    
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i<Nx; i++) {
            I1(i,j,Nz1-1) = I_u1_n(i,j,0) + I_u1z_n(i,j,0);
        }
    }
}

void FDM::baseSolver::Set_H_INV(){
    double a;
    switch (order) {
        case 2:
            a = 0.5;
            break;
        case 4:
            a = 4567./14400.;
            break;
        case 6:
            a = 7493827./25401600.;
            break;
        default:
            a = 0.5;
            break;
    }
    H1 = 1/(a*h);
}

/* ----======== Solver class ========----- */
FDM::solver::solver(double _h, double _dt, int _order) : FDM::baseSolver::baseSolver(_h, _dt, _order){}

Vdata FDM::solver::run(){
    ///// INITIALIZE /////
    
    Vdata u(totalBlocks);
    Vdata u0(totalBlocks);
    Vdata ut0(totalBlocks);
    Vdata u_out(totalBlocks);
    
    int Nz = 0;
    int Nx = 0;
    int Ny = 0;
    for (int lz = 0; lz < LN_Z; lz++) {
        Nz = static_cast<int>( (Z_layers[lz + 1] - Z_layers[lz])/h + 1 );
        
        for (int ly = 0; ly < LN_Y; ly++) {
            Ny = static_cast<int>( (Y_layers[ly + 1] - Y_layers[ly])/h + 1 );
            
            for (int lx = 0; lx < LN_X; lx++) {
                Nx = static_cast<int>( (X_layers[lx+1] - X_layers[lx])/h + 1 );
                int index = idx(lx,ly,lz);
                
                u[index] = myData(Nx,Ny,Nz);
                u0[index] = myData(Nx,Ny,Nz);
                ut0[index] = myData(Nx,Ny,Nz);
                u_out[index] = myData(Nx,Ny,Nz);
            }
        }
    }
    
    if (initializer != NULL)
        initializer(u, *this, init_args);
    if(plotter != NULL)
        plotter(u,*this,0,0, plot_args);
    
    ///// TIME MARCH /////
    double t = 0;
    
    int iterNr = static_cast<int>( tEnd/dt );
    for (int k = 0; k < iterNr; k++) {
        t = k*dt;
        
        Tstep(u_out, u, u0, ut0, k, t);
        u0 = u;
        u = u_out;
        
        ///// PLOT /////
        if(k%SPP == 0 && plotter != NULL)
            plotter(u,*this,k+1,t+dt, plot_args);
    }
    
    if(finilizer != NULL)
        finilizer(u,*this,tEnd, finilizer_args);
    
    return u;
    
}

void FDM::solver::Tstep(Vdata & u_out,const Vdata & u, const Vdata & u0,const Vdata & ut0, int k,double t){
    Vdata rhs4(totalBlocks);
    Vdata F(totalBlocks);
    Vdata Ftt(totalBlocks);
    Vdata F4(totalBlocks);
    
    for (int ln = 0; ln < totalBlocks; ln++){
        int Nx = u[ln].Nx();
        int Ny = u[ln].Ny();
        int Nz = u[ln].Nz();
        
        F[ln]    = myData(Nx,Ny,Nz);
        Ftt[ln]  = myData(Nx,Ny,Nz);
        F4[ln]   = myData(Nx,Ny,Nz);
        rhs4[ln] = myData(Nx,Ny,Nz);
    }
    
    Fource(F, Ftt, *this, u, k, t, fource_args);
    
    Vdata rhs = RHS(u);
    
    if(order > 2){
        rhs4 = RHS(rhs);
        F4   = RHS(F);
        
        for (int ln = 0; ln < totalBlocks; ln++)
            F4[ln] += Ftt[ln];
    }
    
    for (int ln = 0; ln < totalBlocks; ln++) {
        switch (k) {
            case 0:
                u_out[ln] = u[ln] + dt*ut0[ln] + (dt*dt/2)*(rhs[ln] + (dt*dt/12)*(rhs4[ln] + F4[ln]) + F[ln]);
                break;
            default:
                //                u_out[ln] = (2*u[ln] - u0[ln] + dt*dt*(rhs[ln] + (dt*dt/12)*(rhs4[ln] + F4[ln]) + F[ln]));
                u_out[ln] = rhs4[ln];
                u_out[ln] += F4[ln];
                u_out[ln] *= (dt*dt/12);
                u_out[ln] += F[ln];
                u_out[ln] += rhs[ln];
                u_out[ln] *= dt*dt;
                u_out[ln] -= u0[ln];
                u_out[ln] += 2*u[ln];
                break;
        }
    }
}

Vdata FDM::solver::RHS(const Vdata & u){
    Vdata rhs(totalBlocks);
    Vdata ux(totalBlocks);
    Vdata uy(totalBlocks);
    Vdata uz(totalBlocks);
    
    for (int ln = 0; ln < totalBlocks;ln++){
        int Nx = u[ln].Nx();
        int Ny = u[ln].Ny();
        int Nz = u[ln].Nz();
        
        rhs[ln] = myData(Nx,Ny,Nz);
        
        ux[ln] = u[ln].D1(h, order, X_DIR);
        uy[ln] = u[ln].D1(h, order, Y_DIR);
        uz[ln] = u[ln].D1(h, order, Z_DIR);
    }
    
    for (int lz = 0; lz < LN_Z;lz++){
        for (int ly = 0; ly < LN_Y; ly++) {
            for (int lx = 0; lx < LN_X; lx++) {
                int index = idx(lx, ly, lz);
                int Nx = u[index].Nx();
                int Ny = u[index].Ny();
                int Nz = u[index].Nz();
                
                myData Bux(Nx,Ny,Nz);
                myData Buy(Nx,Ny,Nz);
                myData Buz(Nx,Ny,Nz);
                
                myData X_Iu1(Nx,Ny,Nz);
                myData X_Iu2(Nx,Ny,Nz);
                myData Y_Iu1(Nx,Ny,Nz);
                myData Y_Iu2(Nx,Ny,Nz);
                myData Z_Iu1(Nx,Ny,Nz);
                myData Z_Iu2(Nx,Ny,Nz);
                
                //Bux setup
                if(lx == 0){
                    for (int k = 0; k < Nz; k++) {
                        for (int j = 0; j < Ny; j++) {
                            Bux(0,j,k) = -1;
                        }
                    }
                }
                
                if(lx == LN_X - 1){
                    for (int k = 0; k < Nz; k++) {
                        for (int j = 0; j < Ny; j++) {
                            Bux(Nx - 1,j,k) = 1;
                        }
                    }
                }
                //
                
                //Buy setup
                if(ly == 0){
                    for (int k = 0; k < Nz; k++) {
                        for (int i = 0; i < Nx; i++) {
                            Buy(i,0,k) = -1;
                        }
                    }
                }
                if(ly == LN_Y - 1){
                    for (int k = 0; k < Nz; k++) {
                        for (int i = 0; i < Nx; i++) {
                            Buy(i,Ny - 1,k) = 1;
                        }
                    }
                }
                //
                
                //Buz setup
                if(lz == 0){
                    for (int j = 0; j < Ny; j++) {
                        for (int i = 0; i < Nx; i++) {
                            Buz(i,j,0) = -1;
                        }
                    }
                }
                if(lz == LN_Z - 1){
                    for (int j = 0; j < Ny; j++) {
                        for (int i = 0; i < Nx; i++) {
                            Buz(i,j,Nz - 1) = 1;
                        }
                    }
                }
                //
                
                ///// interface
                // x direction
                if (lx != 0) {
                    int index2 = idx(lx - 1, ly, lz);
                    X_Interfc_L(u[index2], ux[index2], speeds[index2], u[index], ux[index], speeds[index], X_Iu1);
                }
                
                if (lx != LN_X - 1) {
                    int index2 = idx(lx + 1, ly, lz);
                    X_Interfc_R(u[index], ux[index], speeds[index], X_Iu2, u[index2], ux[index2], speeds[index2]);
                }
                
                // y direction
                if (ly != 0) {
                    int index2 = idx(lx, ly - 1, lz);
                    Y_Interfc_L(u[index2],uy[index2],speeds[index2],u[index], uy[index], speeds[index], Y_Iu1);
                }
                
                if (ly != LN_Y - 1) {
                    int index2 = idx(lx, ly + 1, lz);
                    Y_Interfc_R(u[index],uy[index],speeds[index],Y_Iu2,u[index2], uy[index2], speeds[index2]);
                }
                
                // z direction
                if (lz != 0) {
                    int index2 = idx(lx, ly, lz - 1);
                    Z_Interfc_L(u[index2],uz[index2],speeds[index2],u[index], uz[index], speeds[index], Z_Iu1);
                }
                
                if (lz != LN_Z - 1) {
                    int index2 = idx(lx, ly, lz + 1);
                    Z_Interfc_R(u[index],uz[index],speeds[index],Z_Iu2,u[index2], uz[index2], speeds[index2]);
                }
                
                myData Dxxu = u[index].D2(h, order, X_DIR);
                myData Dyyu = u[index].D2(h, order, Y_DIR);
                myData Dzzu = u[index].D2(h, order, Z_DIR);
                //                rhs[index] = speeds[index]*(Dxxu + Dyyu + Dzzu) - H1*(speeds[index]*Bux*ux[index] + speeds[index]*Buy*uy[index] + speeds[index]*Buz*uz[index] + X_Iu1 + X_Iu2 + Y_Iu1 + Y_Iu2+ Z_Iu1 + Z_Iu2);
                
                rhs[index] += Dxxu;
                
                rhs[index] += Dyyu;
                
                rhs[index] += Dzzu;
                
                rhs[index] *= speeds[index];
                
                Bux *= ux[index];
                Bux *= speeds[index];
                Buy *= uy[index];
                Buy *= speeds[index];
                Buz *= uz[index];
                Buz *= speeds[index];
                Bux += Buy;
                Bux += Buz;
                Bux += X_Iu1;
                Bux += X_Iu2;
                Bux += Y_Iu1;
                Bux += Y_Iu2;
                Bux += Z_Iu1;
                Bux += Z_Iu2;
                Bux *= H1;
                
                rhs[index] -= Bux;
            }
        }
    }
    
    return rhs;
}

void FDM::solver::Set_Domain(vector<double> & _X_layers, vector<double> & _Y_layers, vector<double> & _Z_layers,vector<double> & _speeds){
    LN_X = (int)_X_layers.size() - 1;
    LN_Y = (int)_Y_layers.size() - 1;
    LN_Z = (int)_Z_layers.size() - 1;
    
    totalBlocks = LN_X*LN_Y*LN_Z;
    
    if (_speeds.size() != totalBlocks) {
        exit(EXIT_FAILURE);
    }
    
    speeds = _speeds;
    X_layers = _X_layers;
    Y_layers = _Y_layers;
    Z_layers = _Z_layers;
}

/* ----======== PML class ========----- */
FDM::pmlSolver::pmlSolver(double _h, double _dt, int _order) : FDM::baseSolver::baseSolver(_h, _dt, _order){
    // Set options to defaults
    
    LN_X = 1;
    LN_Y = 2; // one layer and PML
    LN_Z = 1;
    
    Ta_0 = 5;
    SPP = 10;
    
    PML_Polinomial_Degree = 3;
    PML_Length = 0.5;
    maxSpeed = 1;
    
    totalBlocks = LN_X*LN_Y*LN_Z;
    PML_Blocks = LN_X*LN_Z;
    
    X_layers = vector<double>(LN_X + 1);
    Y_layers = vector<double>(LN_Y + 1);
    Z_layers = vector<double>(LN_Z + 1);
    speeds   = vector<double>(totalBlocks);
    
    X_layers[0] = 0;
    X_layers[1] = 1;
    Y_layers[0] = 0;
    Y_layers[1] = 1;
    Y_layers[2] = 1 + PML_Length;
    Z_layers[0] = 0;
    Z_layers[1] = 1;
    
    speeds[0] = 1;
    speeds[1] = 1;
    
    // Compute Damping coeff
    Set_DampingCoeff();
}

Vdata FDM::pmlSolver::run(){
    ///// INITIALIZE /////
    
    Vdata u(totalBlocks);
    Vdata u0(totalBlocks);
    Vdata ut0(totalBlocks);
    Vdata u_out(totalBlocks);
    
    Vdata v(PML_Blocks);
    Vdata v_out(PML_Blocks);
    
    Vdata w(PML_Blocks);
    Vdata w_out(PML_Blocks);
    
    Vdata p(PML_Blocks);
    Vdata p_out(PML_Blocks);
    
    Vdata sigma(PML_Blocks);
    
    Vdata a_u(PML_Blocks);
    Vdata b_u(PML_Blocks);
    
    int Nz = 0;
    int Nx = 0;
    int Ny = 0;
    for (int lz = 0; lz < LN_Z; lz++) {
        Nz = static_cast<int>( (Z_layers[lz + 1] - Z_layers[lz])/h + 1 );
        
        for (int ly = 0; ly < LN_Y; ly++) {
            Ny = static_cast<int>( (Y_layers[ly + 1] - Y_layers[ly])/h + 1 );
            
            for (int lx = 0; lx < LN_X; lx++) {
                Nx = static_cast<int>( (X_layers[lx + 1] - X_layers[lx])/h + 1 );
                int index = idx(lx,ly,lz);
                u[index]     = myData(Nx,Ny,Nz);
                u0[index]    = myData(Nx,Ny,Nz);
                ut0[index]   = myData(Nx,Ny,Nz);
                u_out[index] = myData(Nx,Ny,Nz);
            }
        }
    }
    
    
    //initializing PML
    int PML_Ny = static_cast<int>( PML_Length/h + 1);
    for (int lz = 0; lz < LN_Z; lz++) {
        Nz = static_cast<int>( (Z_layers[lz + 1] - Z_layers[lz])/h + 1 );
        
        for (int lx = 0; lx < LN_X; lx++) {
            Nx = static_cast<int>( (X_layers[lx + 1] - X_layers[lx])/h + 1 );
            int index = idx(lx,LN_Y-1,lz);
            int PML_index = PML_idx(lx,lz);
            
            u[index]            = myData(Nx,PML_Ny,Nz);
            u0[index]           = myData(Nx,PML_Ny,Nz);
            ut0[index]          = myData(Nx,PML_Ny,Nz);
            u_out[index]        = myData(Nx,PML_Ny,Nz);
            
            v[PML_index]        = myData(Nx,PML_Ny,Nz);
            v_out[PML_index]    = myData(Nx,PML_Ny,Nz);
            
            w[PML_index]        = myData(Nx,PML_Ny,Nz);
            w_out[PML_index]    = myData(Nx,PML_Ny,Nz);
            
            p[PML_index]        = myData(Nx,PML_Ny,Nz);
            p_out[PML_index]    = myData(Nx,PML_Ny,Nz);
            
            sigma[PML_index]    = sigmaY(u[index]);
            
            a_u[PML_index]      = 1/(1 + 0.5*dt*sigma[PML_index]);
            b_u[PML_index]      = 1 - 0.5*dt*sigma[PML_index];
        }
    }
    
    if (initializer != NULL)
        initializer(u, *this, init_args);
    if(plotter != NULL)
        plotter(u,*this,0,0, plot_args);
    
    ///// TIME MARCH /////
    double t = 0;
    
    int iterNr = static_cast<int>( tEnd/dt );
    for (int k = 0; k < iterNr; k++) {
        t = k*dt;
        
        Tstep(u_out, u, u0, ut0, v_out, v, w_out, w, p_out, p, sigma, a_u, b_u, k, t);
        u0 = u;
        u = u_out;
        v = v_out;
        w = w_out;
        p = p_out;
        
        ///// PLOT /////
        if(k%SPP == 0 && plotter != NULL)
            plotter(u,*this,k+1,t+dt, plot_args);
    }
    
    if(finilizer != NULL)
        finilizer(u,*this,tEnd, finilizer_args);
    
    return u;
    
}

void FDM::pmlSolver::Tstep(Vdata & u_out,const Vdata & u, const Vdata & u0,const Vdata & ut0, Vdata & v_out,const Vdata & v, Vdata & w_out,const Vdata & w, Vdata & p_out,const Vdata & p, const Vdata & sigma, const Vdata & a_u, const Vdata & b_u, int k,double t){
    Vdata rhs4(totalBlocks);
    Vdata F(totalBlocks);
    Vdata Ftt(totalBlocks);
    Vdata F4(totalBlocks);
    Vdata sigma4(PML_Blocks);
    
    for (int ln = 0; ln < totalBlocks; ln++){
        int Nx = u[ln].Nx();
        int Ny = u[ln].Ny();
        int Nz = u[ln].Nz();
        
        F[ln]    = myData(Nx,Ny,Nz);
        Ftt[ln]  = myData(Nx,Ny,Nz);
        F4[ln]   = myData(Nx,Ny,Nz);
        rhs4[ln] = myData(Nx,Ny,Nz);
    }
    
    for (int lz = 0; lz < LN_Z; lz++) {
        for (int lx = 0; lx < LN_X; lx++) {
            int index = idx(lx,LN_Y-1,lz);
            int PML_index = PML_idx(lx,lz);
            int Nx = u[index].Nx();
            int Ny = u[index].Ny();
            int Nz = u[index].Nz();
            
            sigma4[PML_index] = myData(Nx,Ny,Nz);
        }
    }
    
    Fource(F, Ftt, *this, u, k, t, fource_args);
    
    Vdata rhs = RHS(u, sigma, w, p);
    
    if(order > 2){
        rhs4 = RHS(rhs, sigma4, w, p);
        F4   = RHS(F,sigma4, w, p);
        
        for (int ln = 0; ln < totalBlocks; ln++)
            F4[ln] += Ftt[ln];
    }
    
    for (int lz = 0; lz < LN_Z; lz++) {
        for (int ly = 0; ly < LN_Y - 1; ly++) {
            for (int lx = 0; lx < LN_X; lx++) {
                int ln = idx(lx, ly, lz);
                switch (k) {
                    case 0:
                        u_out[ln] = u[ln] + dt*ut0[ln] + (dt*dt/2)*(rhs[ln] + (dt*dt/12)*(rhs4[ln] + F4[ln]) + F[ln]);
                        break;
                    default:
                        //                u_out[ln] = (2*u[ln] - u0[ln] + dt*dt*(rhs[ln] + (dt*dt/12)*(rhs4[ln] + F4[ln]) + F[ln]));
                        u_out[ln] = rhs4[ln];
                        u_out[ln] += F4[ln];
                        u_out[ln] *= (dt*dt/12);
                        u_out[ln] += F[ln];
                        u_out[ln] += rhs[ln];
                        u_out[ln] *= dt*dt;
                        u_out[ln] -= u0[ln];
                        u_out[ln] += 2*u[ln];
                        break;
                }
            }
        }
    }
    
    for (int lz = 0; lz < LN_Z; lz++) {
        for (int lx = 0; lx < LN_X; lx++) {
            
            int ln = idx(lx,LN_Y-1,lz);
            int PML_ln = PML_idx(lx,lz);
            
            switch (k) {
                case 0:
                    u_out[ln] = u[ln] + dt*ut0[ln] + (dt*dt/2)*(rhs[ln] + (dt*dt/12)*(rhs4[ln] + F4[ln]) + F[ln]);
                    break;
                default:
                    //                    u_out[ln] = (2*u[ln] - u0[ln] + dt*dt*(rhs[ln] + (dt*dt/12)*(rhs4[ln] + F4[ln]) + F[ln]));
                    u_out[ln] = rhs4[ln];
                    u_out[ln] += F4[ln];
                    u_out[ln] *= (dt*dt/12);
                    u_out[ln] -= speeds[ln]*(sigma[PML_ln]*v[PML_ln]).D1(h,order, Y_DIR);
                    u_out[ln] += speeds[ln]*(sigma[PML_ln]*w[PML_ln]).D1(h,order, X_DIR);
                    u_out[ln] += speeds[ln]*(sigma[PML_ln]*p[PML_ln]).D1(h,order, Z_DIR);
                    u_out[ln] += F[ln];
                    u_out[ln] += rhs[ln];
                    u_out[ln] *= dt*dt;
                    u_out[ln] -= b_u[PML_ln]*u0[ln];
                    u_out[ln] += 2*u[ln];
                    u_out[ln] *= (a_u[PML_ln]);
                    break;
            }
            
            v_out[PML_ln] = a_u[PML_ln]*(b_u[PML_ln]*v[PML_ln] + (dt/2)*(u_out[ln] + u[ln]).D1(h,order, Y_DIR));
            
            w_out[PML_ln] = (u_out[ln] + u[ln]).D1(h,order, X_DIR);
            w_out[PML_ln] *= (dt/2);
            w_out[PML_ln] += w[PML_ln];
            
            p_out[PML_ln] = (u_out[ln] + u[ln]).D1(h,order, Z_DIR);
            p_out[PML_ln] *= (dt/2);
            p_out[PML_ln] += p[PML_ln];
        }
    }
}

Vdata FDM::pmlSolver::RHS(const Vdata & u, const Vdata & sigma, const Vdata & w, const Vdata & p){
    Vdata rhs(totalBlocks);
    Vdata ux(totalBlocks);
    Vdata uy(totalBlocks);
    Vdata uz(totalBlocks);
    
    for (int ln = 0; ln < totalBlocks;ln++){
        int Nx = u[ln].Nx();
        int Ny = u[ln].Ny();
        int Nz = u[ln].Nz();
        
        rhs[ln] = myData(Nx,Ny,Nz);
        
        ux[ln] = u[ln].D1(h, order, X_DIR);
        uy[ln] = u[ln].D1(h, order, Y_DIR);
        uz[ln] = u[ln].D1(h, order, Z_DIR);
    }
    
    for (int lz = 0; lz < LN_Z;lz++){
        for (int ly = 0; ly < LN_Y - 1; ly++) {
            for (int lx = 0; lx < LN_X; lx++) {
                int index = idx(lx, ly, lz);
                int Nx = u[index].Nx();
                int Ny = u[index].Ny();
                int Nz = u[index].Nz();
                
                myData Bux(Nx,Ny,Nz);
                myData Buy(Nx,Ny,Nz);
                myData Buz(Nx,Ny,Nz);
                
                myData X_Iu1(Nx,Ny,Nz);
                myData X_Iu2(Nx,Ny,Nz);
                myData Y_Iu1(Nx,Ny,Nz);
                myData Y_Iu2(Nx,Ny,Nz);
                myData Z_Iu1(Nx,Ny,Nz);
                myData Z_Iu2(Nx,Ny,Nz);
                
                //Bux setup
                if(lx == 0){
                    for (int k = 0; k < Nz; k++) {
                        for (int j = 0; j < Ny; j++) {
                            Bux(0,j,k) = -1;
                        }
                    }
                }
                
                if(lx == LN_X - 1){
                    for (int k = 0; k < Nz; k++) {
                        for (int j = 0; j < Ny; j++) {
                            Bux(Nx - 1,j,k) = 1;
                        }
                    }
                }
                //
                
                //Buy setup
                if(ly == 0){
                    for (int k = 0; k < Nz; k++) {
                        for (int i = 0; i < Nx; i++) {
                            Buy(i,0,k) = -1;
                        }
                    }
                }
                //
                
                //Buz setup
                if(lz == 0){
                    for (int j = 0; j < Ny; j++) {
                        for (int i = 0; i < Nx; i++) {
                            Buz(i,j,0) = -1;
                        }
                    }
                }
                if(lz == LN_Z - 1){
                    for (int j = 0; j < Ny; j++) {
                        for (int i = 0; i < Nx; i++) {
                            Buz(i,j,Nz - 1) = 1;
                        }
                    }
                }
                //
                
                ///// interface
                // x direction
                if (lx != 0) {
                    int index2 = idx(lx - 1, ly, lz);
                    X_Interfc_L(u[index2], ux[index2], speeds[index2], u[index], ux[index], speeds[index], X_Iu1);
                }
                
                if (lx != LN_X - 1) {
                    int index2 = idx(lx + 1, ly, lz);
                    X_Interfc_R(u[index], ux[index], speeds[index], X_Iu2, u[index2], ux[index2], speeds[index2]);
                }
                
                // y direction
                if (ly != 0) {
                    int index2 = idx(lx, ly - 1, lz);
                    Y_Interfc_L(u[index2],uy[index2],speeds[index2],u[index], uy[index], speeds[index], Y_Iu1);
                }
                
//                if (ly != LN_Y - 1) { //always happens in this loop
                int index2 = idx(lx, ly + 1, lz); 
                Y_Interfc_R(u[index],uy[index],speeds[index],Y_Iu2,u[index2], uy[index2], speeds[index2]);
//                }
                
                // z direction
                if (lz != 0) {
                    int index2 = idx(lx, ly, lz - 1);
                    Z_Interfc_L(u[index2],uz[index2],speeds[index2],u[index], uz[index], speeds[index], Z_Iu1);
                }
                
                if (lz != LN_Z - 1) {
                    int index2 = idx(lx, ly, lz + 1);
                    Z_Interfc_R(u[index],uz[index],speeds[index],Z_Iu2,u[index2], uz[index2], speeds[index2]);
                }
                
                myData Dxxu = u[index].D2(h, order, X_DIR);
                myData Dyyu = u[index].D2(h, order, Y_DIR);
                myData Dzzu = u[index].D2(h, order, Z_DIR);
//                rhs[index] = speeds[index]*(Dxxu + Dyyu + Dzzu) - H1*(speeds[index]*Bux*ux[index] + speeds[index]*Buy*uy[index] + speeds[index]*Buz*uz[index] + X_Iu1 + X_Iu2 + Y_Iu1 + Y_Iu2+ Z_Iu1 + Z_Iu2);
                
                rhs[index] += Dxxu;
                
                rhs[index] += Dyyu;
                
                rhs[index] += Dzzu;
                
                rhs[index] *= speeds[index];

                ux[index] *= speeds[index];
                uy[index] *= speeds[index];
                uz[index] *= speeds[index];

                Bux *= ux[index];
                Buy *= uy[index];
                Buz *= uz[index];
                Bux += Buy;
                Bux += Buz;
                Bux += X_Iu1;
                Bux += X_Iu2;
                Bux += Y_Iu1;
                Bux += Y_Iu2;
                Bux += Z_Iu1;
                Bux += Z_Iu2;
                Bux *= H1;
                
                rhs[index] -= Bux;
            }
        }
    }
    
    for (int lz = 0; lz < LN_Z;lz++){
        for (int lx = 0; lx < LN_X; lx++) {
            int ly = LN_Y - 1;
            int index = idx(lx, ly, lz);
            int Nx = u[index].Nx();
            int Ny = u[index].Ny();
            int Nz = u[index].Nz();
            
            myData Bux(Nx,Ny,Nz);
            myData Buy(Nx,Ny,Nz);
            myData Buz(Nx,Ny,Nz);
            
            myData X_Iu1(Nx,Ny,Nz);
            myData X_Iu2(Nx,Ny,Nz);
            myData Y_Iu1(Nx,Ny,Nz);
            myData Y_Iu2(Nx,Ny,Nz);
            myData Z_Iu1(Nx,Ny,Nz);
            myData Z_Iu2(Nx,Ny,Nz);
            
            //Bux setup
            if(lx == 0){
                for (int k = 0; k < Nz; k++) {
                    for (int j = 0; j < Ny; j++) {
                        Bux(0,j,k) = -1;
                    }
                }
            }
            
            if(lx == LN_X - 1){
                for (int k = 0; k < Nz; k++) {
                    for (int j = 0; j < Ny; j++) {
                        Bux(Nx - 1,j,k) = 1;
                    }
                }
            }
            //
             
            
//            if(ly == LN_Y - 1){ //always happens
            for (int k = 0; k < Nz; k++) {
                for (int i = 0; i < Nx; i++) {
                    Buy(i,Ny - 1,k) = 1;
                }
            }
//            }
            //
                
            //Buz setup
            if(lz == 0){
                for (int j = 0; j < Ny; j++) {
                    for (int i = 0; i < Nx; i++) {
                        Buz(i,j,0) = -1;
                    }
                }
            }
            if(lz == LN_Z - 1){
                for (int j = 0; j < Ny; j++) {
                    for (int i = 0; i < Nx; i++) {
                        Buz(i,j,Nz - 1) = 1;
                    }
                }
            }
            //
                
                ///// interface
                // x direction
            if (lx != 0) {
                int index2 = idx(lx - 1, ly, lz);
                int Pidx = PML_idx(lx, lz);
                int Pidx2 = PML_idx(lx - 1, lz);
                X_Interfc_L(u[index2], ux[index2] + sigma[Pidx2]*w[Pidx2], speeds[index2], u[index], ux[index] + sigma[Pidx]*w[Pidx], speeds[index], X_Iu1);
            }
            
            if (lx != LN_X - 1) {
                int index2 = idx(lx + 1, ly, lz);
                int Pidx = PML_idx(lx, lz);
                int Pidx2 = PML_idx(lx + 1, lz);
                X_Interfc_R(u[index], ux[index] + sigma[Pidx]*w[Pidx], speeds[index], X_Iu2, u[index2], ux[index2] + sigma[Pidx2]*w[Pidx2], speeds[index2]);
            }
            
            // y direction
            if (ly != 0) {
                int index2 = idx(lx, ly - 1, lz);
                Y_Interfc_L(u[index2],uy[index2],speeds[index2],u[index], uy[index], speeds[index], Y_Iu1);
            }
            
//            if (ly != LN_Y - 1) { //never happens
//                int index2 = idx(lx, ly + 1, lz);
//                Y_Interfc_R(u[index],uy[index],speeds[index],Y_Iu2,u[index2], uy[index2], speeds[index2]);
//            }
                
            // z direction
            if (lz != 0) {
                int index2 = idx(lx, ly, lz - 1);
                int Pidx = PML_idx(lx, lz);
                int Pidx2 = PML_idx(lx, lz - 1);
                Z_Interfc_L(u[index2],uz[index2] + sigma[Pidx2]*p[Pidx2],speeds[index2],u[index], uz[index] + sigma[Pidx]*p[Pidx], speeds[index], Z_Iu1);
            }
                
            if (lz != LN_Z - 1) {
                int index2 = idx(lx, ly, lz + 1);
                int Pidx = PML_idx(lx, lz);
                int Pidx2 = PML_idx(lx, lz + 1);
                Z_Interfc_R(u[index],uz[index] + sigma[Pidx]*p[Pidx],speeds[index],Z_Iu2,u[index2], uz[index2] + sigma[Pidx2]*p[Pidx2], speeds[index2]);
            }
                
            myData Dxxu = u[index].D2(h, order, X_DIR);
            myData Dyyu = u[index].D2(h, order, Y_DIR);
            myData Dzzu = u[index].D2(h, order, Z_DIR);
            //                rhs[index] = speeds[index]*(Dxxu + Dyyu + Dzzu) - H1*(speeds[index]*Bux*ux[index] + speeds[index]*Buy*uy[index] + speeds[index]*Buz*uz[index] + X_Iu1 + X_Iu2 + Y_Iu1 + Y_Iu2+ Z_Iu1 + Z_Iu2);
            
            rhs[index] += Dxxu;
            
            rhs[index] += Dyyu;
            
            rhs[index] += Dzzu;
            
            rhs[index] *= speeds[index];
            
            ux[index] *= speeds[index];
            uy[index] *= speeds[index];
            uz[index] *= speeds[index];
            
            int Pidx = PML_idx(lx, lz);
            ux[index] += sigma[Pidx]*w[Pidx];
            uz[index] += sigma[Pidx]*p[Pidx];
            
            Bux *= ux[index];
            Buy *= uy[index];
            Buz *= uz[index];
            Bux += Buy;
            Bux += Buz;
            Bux += X_Iu1;
            Bux += X_Iu2;
            Bux += Y_Iu1;
            Bux += Y_Iu2;
            Bux += Z_Iu1;
            Bux += Z_Iu2;
            Bux *= H1;
            
            rhs[index] -= Bux;
        }
    }

    return rhs;
}

myData FDM::pmlSolver::sigmaY(const myData& u){
    int Nx = u.Nx();
    int Ny = u.Ny();
    int Nz = u.Nz();
    myData sigma(Nx,Ny,Nz);
    
    double sig = 0;
    for (int j = 0; j < Ny; j++) {
        sig = PML_Damping_Coef*pow( (j*h)/PML_Length, PML_Polinomial_Degree);
        for (int k = 0; k < Nz; k++) {
            for (int i = 0; i < Nx; i++) {
                sigma(i,j,k) = sig;
            }
        }
    }
    return sigma;
}

void FDM::pmlSolver::Set_Domain(vector<double> & _X_layers, vector<double> & _Y_layers, vector<double> & _Z_layers,vector<double> & _speeds){
    LN_X = (int)_X_layers.size() - 1;
    LN_Y = (int)_Y_layers.size() - 1;
    LN_Z = (int)_Z_layers.size() - 1;
    
    totalBlocks = LN_X*LN_Y*LN_Z;
    
    if (_speeds.size() != totalBlocks) {
        exit(EXIT_FAILURE);
    }
    
    LN_Y++;
    totalBlocks = LN_X*LN_Y*LN_Z;
    PML_Blocks = LN_X*LN_Z;
    
    speeds = vector<double>(totalBlocks);
    Y_layers = vector<double>(LN_Y + 1);
    X_layers = _X_layers;
    Z_layers = _Z_layers;
    
    for (int ly = 0; ly < LN_Y; ly++) 
        Y_layers[ly] = _Y_layers[ly];

    Y_layers[LN_Y] = Y_layers[LN_Y - 1] + PML_Length;

    for (int lz = 0; lz < LN_Z; lz++){
        for (int lx = 0; lx < LN_X; lx++) {
            for (int ly = 0; ly < LN_Y - 1; ly++){ 
                int ln = idx(lx, ly, lz);
                int tln = lx + ly*LN_X + lz*LN_X*(LN_Y - 1);
                speeds[ln] = _speeds[tln];
            }
            int PMLidx = idx(lx, LN_Y - 1, lz);
            int pidx = idx(lx, LN_Y - 2, lz);
            speeds[PMLidx] = speeds[pidx];
        }
    }

    
    for (int ln = 0; ln < totalBlocks; ln++) {
        if (maxSpeed < speeds[ln])
            maxSpeed = speeds[ln];
    }
}

void FDM::pmlSolver::Set_PML( double _PML_Length, int _dampingPDegree){
    PML_Damping_Coef = _dampingPDegree;
    PML_Length = _PML_Length;
    
    Y_layers[LN_Y] = Y_layers[LN_Y - 1] + PML_Length;
}

void FDM::pmlSolver::Set_DampingCoeff(bool relativeToH ,double _DampingCoeff){
    if (!relativeToH) {
        PML_Damping_Coef = _DampingCoeff;
        return;
    }
    double C_0 = 0.0001;
    double tol = pow((C_0*h),order);
    PML_Damping_Coef = maxSpeed*(PML_Polinomial_Degree+1)/(2*PML_Length)*(log(1/tol)/log(2));
}