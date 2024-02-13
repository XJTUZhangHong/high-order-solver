#pragma once
#include "1Dproblem.h"

void PlanarShock();

void IC_for_riemann_2d(Fluid2d* fluid, double* zone1, double* zone2, double* zone3, double* zone4, Block2d block, double* discon);