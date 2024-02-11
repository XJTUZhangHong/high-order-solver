#pragma once
#include "basic_function.h"

// one-dimensional problem
enum GKS1d_type{nothing, kfvs1st, kfvs2nd, gks1st, gks2nd};
extern GKS1d_type gks1dsolver;

void Calculate_flux(Flux1d** fluxes, Interface1d* interfaces, Block1d &block, int stage);
typedef void(*Flux_function)(Flux1d &flux, Interface1d& interface, double dt);
extern Flux_function flux_function;
void LF(Flux1d& flux, Interface1d& interface, double dt);
void get_Euler_flux(double p[3], double* flux);
