#pragma once
#include "1Dproblem.h"

void PlanarShock();

void PlanarSheer();

void IC_for_riemann_2d(Fluid2d* fluid, double* zone1, double* zone2, double* zone3, double* zone4, Block2d block, double* discon);

void High_mach_astrophusical_jet();

void IC_for_astrophusical_jet(Fluid2d* fluid, Block2d block);

void High_mach_astrophusical_jet();

void IC_for_astrophusical_jet(Fluid2d* fluid, Block2d block);

void inflow_boundary_left(Fluid2d* fluids, Block2d block);

void RT_instability();

void IC_for_RT_instability(Fluid2d* fluid, Block2d block);

void doubleMach();

void ICforDoubleMach(Fluid2d* fluid, Block2d block);

void boundaryforDoubleMach(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);

void viscous_sod_shock_problem();

void IC_for_vishocktube(int order, Fluid2d* fluid, Block2d block);

void accuracy_sinwave_2d();

void sinwave_2d(double& CFL, double& dt_ratio, int& mesh_number, double* error);

void ICfor_sinwave_2d(Fluid2d* fluids, Block2d block);

void error_for_sinwave_2d(Fluid2d* fluids, Block2d block, double tstop, double* error);

