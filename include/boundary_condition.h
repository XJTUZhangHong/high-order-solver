#pragma once
#include "basic_function.h"

// one-dimensional problem
typedef void(*BoundaryCondition) (Fluid1d *fluids, Block1d block, Fluid1d bcvalue);
void free_boundary_left(Fluid1d *fluids, Block1d block, Fluid1d bcvalue);
void free_boundary_right(Fluid1d *fluids, Block1d block, Fluid1d bcvalue);