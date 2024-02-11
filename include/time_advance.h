#pragma once
#include "basic_function.h"

// one-dimensional problem
typedef void(*TimeMarchingCoefficient)(Block1d &block);
extern TimeMarchingCoefficient timecoe_list;
void S1O1(Block1d& block);
void S1O2(Block1d& block);
void S2O4(Block1d& block);
void RK2(Block1d& block);
void RK3(Block1d& block);
void RK4(Block1d &block);

void Initial_stages(Block1d &block);

double Get_CFL(Block1d& block, Fluid1d* fluids, double tstop);

double Dtx(double dtx, double dx, double CFL, double convar[3]);

void Update(Fluid1d* fluids, Flux1d** fluxes, Block1d block, int stage);