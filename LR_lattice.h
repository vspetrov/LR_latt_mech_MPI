#ifndef LR_LATTICE_H
#define LR_LATTICE_H
#include "parallel.h"
extern int Size; //Lattice dimensions


struct CellVariables //Phase space variables of one cell
{
	double m;
	double h;
	double j;
	double d;
	double f;
	double X;
	double Cai;
    double tension;
    double s_tension;
};

extern double *V;
extern double *gk1;
extern double *Gk;
extern double *Gsi;
extern double *type;

extern CellVariables *Cell; 
extern double *I_ext;//External current
extern const double dt;//Integrating step (variable)
extern double D;//coupling parameter 
extern int DrawNum;
double Coupling(int &I, int &J); // This function returns the value of coupling D*(...) of the cell [i,j]
extern void LattInit(const int &size);//initializing function
void OdeSolve(int &ii, int &jj);//this function solves the system with NO coupling on the interval 'dt'
inline int Substeps(double &vd);//devides step length due to value of Voltage (V)
extern void SolveEquations(double MaxTime /*time of calculations*/, MPI_Comm &GridComm, int &GridSize, MPI_Datatype &column, MPI_Datatype &square);//Solves the task
extern void LattInit(int &size, char a[]);
void Latt_Init1(int &size, char a[]);
#endif 
