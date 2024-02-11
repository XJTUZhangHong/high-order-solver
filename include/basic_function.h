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

double Pressure(double density, double densityu, double densityE);

double U(double density, double q_densityu);