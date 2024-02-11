#pragma once
#include "fluid_mesh.h"
#include <cmath>
#include <assert.h>

extern int K;
extern double Gamma;
extern double Mu;
extern double Nu;
extern double T_inf; //the real temperature
extern double c1_euler;
extern double c2_euler;
extern bool is_Prandtl_fix;
extern double Pr;
extern double R_gas; //the gas constant for real ideal air
enum TAU_TYPE { Euler, NS};
extern TAU_TYPE tau_type;

//basic gks function
// to store the moment
class MMDF1d
{
private:
	double u;
	double lambda;

public:
	double uwhole[10];
	double uplus[10];
	double uminus[10];
	double upxi[10][4];
	double unxi[10][4];
	double uxi[10][4];
	double xi2;
	double xi4;
	double xi6;
	MMDF1d();
	MMDF1d(double u_in, double lambda_in);
	void calcualte_MMDF1d();
};

// to calculate the microsolpe moment
void G(int no_u, int no_xi, double* psi, double a[3], MMDF1d m);
void GL(int no_u, int no_xi, double* psi, double a[3], MMDF1d m);
void GR(int no_u, int no_xi, double* psi, double a[3], MMDF1d m);

void Microslope(double* a, double der[3], double prim[3]);

double Get_Tau_NS(double density0, double lambda0);

double Get_Tau(double density_left, double density_right, double density0, double lambda_left, double lambda_right, double lambda0, double dt);

double Alpha(double lambda, double u);

double Beta(double lambda, double u);

void Convar_to_primvar_1D(double* primvar, double convar[3]);

void Convar_to_ULambda_1d(double* primvar, double convar[3]);

void Primvar_to_convar_1D(double* convar, double primvar[3]);

double DensityU(double density, double u);

double DensityE(double density, double u, double pressure);

double Pressure(double density, double densityu, double densityE);

double U(double density, double q_densityu);

double Lambda(double density, double u, double densityE);

Flux1d** Setflux_array(Block1d block);

void SetUniformMesh(Block1d block, Fluid1d* fluids, Interface1d* interfaces, Flux1d** fluxes);

void CopyFluid_new_to_old(Fluid1d* fluids, Block1d block);