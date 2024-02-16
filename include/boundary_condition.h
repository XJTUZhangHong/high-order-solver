#pragma once
#include "basic_function.h"

// one-dimensional problem
typedef void(*BoundaryCondition) (Fluid1d *fluids, Block1d block, Fluid1d bcvalue);
void free_boundary_left(Fluid1d *fluids, Block1d block, Fluid1d bcvalue);
void free_boundary_right(Fluid1d *fluids, Block1d block, Fluid1d bcvalue);
void reflection_boundary_left(Fluid1d* fluids, Block1d block, Fluid1d bcvalue);
void reflection_boundary_right(Fluid1d* fluids, Block1d block, Fluid1d bcvalue);
void periodic_boundary_left(Fluid1d* fluids, Block1d block, Fluid1d bcvalue);
void periodic_boundary_right(Fluid1d* fluids, Block1d block, Fluid1d bcvalue);
// two-dimensional problem
typedef void(*BoundaryCondition2d) (Fluid2d *fluids, Block2d block, Fluid2d bcvalue);
void free_boundary_left(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void free_boundary_right(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void free_boundary_down(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void free_boundary_up(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void periodic_boundary_left(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void periodic_boundary_right(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void periodic_boundary_down(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void periodic_boundary_up(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void RT_boundary(Fluid2d* fluids, Block2d block);

