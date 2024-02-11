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