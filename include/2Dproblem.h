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



