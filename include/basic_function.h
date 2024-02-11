#pragma once
#include "fluid_mesh.h"
#include <cmath>
#include <assert.h>

extern int K;
extern double Gamma;
extern double Mu;
extern double Nu;
extern double T_inf; //the real temperature
extern double R_gas; //the gas constant for real ideal air
enum TAU_TYPE { Euler, NS};
extern TAU_TYPE tau_type;