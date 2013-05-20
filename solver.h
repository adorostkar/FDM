////  Created by Ali Dorostkar on 3/24/12.
////  Copyright (c) 2012 Uppsala. All rights reserved.
////

// Solver class
// -- Solver options:
// ---- Order
// ---- Number of layers in each direction
// ---- Spatial step size
// ---- time step size
// ---- final time
// ---- Penalty terms(SAT arguments)
// ---- Speeds of each layer
// ---- PML length
// ---- PML damping coefficient
// ---- PML polynomial degree
// ---- Pointer to plot function
// ---- pointer to initializer
// ---- pointer to finalizer

// TODO - Solve speeds issue, a better way to get them.
// TODO - Add options class which configures the solver.
// TODO - Change to generic coding for the dimensions.
// TODO - Change the structure of run method from for loop to while
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "myData.h"

#define Vdata vector<myData>

using namespace std;
#ifndef PDE_solver_h
#define PDE_solver_h

namespace FDM {
    class baseSolver{
    protected:
        int order;
        int LN_X;
        int LN_Y;
        int LN_Z;
        int SPP;
        double h;
        double H1;
        double dt;
        double tEnd;
        double Ta_0;
        int totalBlocks;
        std::vector<double> X_layers;
        std::vector<double> Y_layers;
        std::vector<double> Z_layers;
        std::vector<double> speeds;
        
        ///// SOURCE AND PLOT HANDELING /////
        void * plot_args;
        void (*plotter)(const Vdata &, baseSolver& ,int, double , void *);
        void * fource_args;
        void (*Fource)(Vdata &, Vdata &, baseSolver &, const Vdata&, int, double, void *);
        void * init_args;
        void (*initializer)(Vdata &, baseSolver &, void *);
        void * finilizer_args;
        void (*finilizer)(Vdata &, baseSolver &, double, void *);
        
        void X_Interfc_L(const myData & u1, const myData & u1x,const double c1 ,const myData & u2, const myData & u2x,const double c2,myData & I2);
        void X_Interfc_R(const myData & u1, const myData & u1x,const double c1,myData & I1 ,const myData & u2, const myData & u2x,const double c2);
        
        void Y_Interfc_L(const myData & u1, const myData & u1y,const double c1 ,const myData & u2, const myData & u2y,const double c2,myData & I2);
        void Y_Interfc_R(const myData & u1, const myData & u1y,const double c1,myData & I1 ,const myData & u2, const myData & u2y,const double c2);
        
        void Z_Interfc_L(const myData & u1, const myData & u1z,const double c1 ,const myData & u2, const myData & u2z,const double c2,myData & I2);
        void Z_Interfc_R(const myData & u1, const myData & u1z,const double c1,myData & I1 ,const myData & u2, const myData & u2z,const double c2);
        
        void Set_H_INV();
        
    public:
        baseSolver(double _h = 0.01, double _dt = 0.001, int _order = 2);
        
        ///// METHODS /////
        const double X(const int i,int ln){ return X_layers[ln] + i*h;};
        const double Y(const int i,int ln){ return Y_layers[ln] + i*h;};
        const double Z(const int i,int ln){ return Z_layers[ln] + i*h;};
        const int Order(){ return order;};
        const int X_LayerNumber(){ return LN_X;};
        const int Y_LayerNumber(){ return LN_Y;};
        const int Z_LayerNumber(){ return LN_Z;};
        const double dx(){ return h;};
        const double Dt(){ return dt;};
        const double EndTime(){ return tEnd;};
        const int TotalTimeIteration(){ return static_cast<int>(tEnd/dt) - 1;};
        const double H_INV(){ return H1;};
        const double TA0(){return Ta_0;};
        
        void Set_Plotter(void (*SP)(const Vdata &, baseSolver&,int, double , void *),void *);
        void Set_FourceFunction(void (*SF)(Vdata &, Vdata &, baseSolver & ,const Vdata&, int, double, void *), void *);
        void Set_Initializer(void (*SI)(Vdata &, baseSolver &, void *), void *);
        void Set_Finilizer(void (*FI)(Vdata &, baseSolver &, double, void *), void *);
        void Set_StepPerPlot(int _SPP){SPP = _SPP;};
        void Set_dx(double _h){h = _h;Set_H_INV();};
        void Set_dt(double _dt){dt = _dt;};
        void Set_Ta0(double _Ta_0){Ta_0 = _Ta_0;};
        void Set_Order(int _order){order = _order;Set_H_INV();};
        void Set_EndTime(double _t){tEnd = _t;};
    };
    
    class solver : public baseSolver{
    private:
        
        ///// METHODS /////
        void Tstep(Vdata & u_out,const Vdata & u, const Vdata & u0,const Vdata & ut0,int k, double t);
        
        Vdata RHS(const Vdata & u);
    public:
        ///// CONSTRUCTORS /////
        solver(double _h = 0.01, double _dt = 0.001, int _order = 2);
        void Set_Domain(vector<double> & _X_layers, vector<double> & _Y_layers, vector<double> & _Z_layers,vector<double> & _speeds);
        
        Vdata run();
    };
    
    class pmlSolver : public baseSolver{
    private:
        int PML_Blocks;                 // Number of PML blocks
        
        /// PML VAR ///
        double  PML_Length;             // length of the PML
        double  maxSpeed;               // Maximum speed between layers
        double  PML_Damping_Coef;       // PML damping coefficient
        int     PML_Polinomial_Degree;  // PML polynomial degree
        
        ///// METHODS /////
        // Time stepping method
        void Tstep(Vdata & u_out,const Vdata & u, const Vdata & u0,const Vdata & ut0, Vdata & v_out,const Vdata & v, Vdata & w_out,const Vdata & w, Vdata & p_out,const Vdata & p, const Vdata & sigma, const Vdata & a_u, const Vdata & b_u, int k,double t);
        
        // Computes right hand side of the discreet problem consisting SAT terms
        Vdata RHS(const Vdata & u, const Vdata & sigma, const Vdata & w, const Vdata & p);
    public:
        ///// CONSTRUCTORS /////
        pmlSolver(double _h = 0.01, double _dt = 0.001, int _order = 2);
        
        ///// METHODS /////
        void Set_Domain(vector<double> & _X_layers, vector<double> & _Y_layers, vector<double> & _Z_layers,vector<double> & _speeds);
        
        // Compute sigma in Y direction
        myData sigmaY(const myData& u);
        
        // Set PML length and degree of freedom
        void Set_PML( double _PML_Length, int _dampingPDegree);
        
        // If reletive to h, compute. If not set to value with diffult equal zero
        void Set_DampingCoeff(bool relativeToH = true ,double _DampingCoeff = 0);
        
        // Run experiment and return final result
        Vdata run();
    };
}

#endif
