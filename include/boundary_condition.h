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
void inflow_boundary_left(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void inflow_boundary_right(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void inflow_boundary_up(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void noslip_adiabatic_boundary_left(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void noslip_adiabatic_boundary_right(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void noslip_adiabatic_boundary_down(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void noslip_adiabatic_boundary_up(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void reflection_boundary_up(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void reflection_boundary_down(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void reflection_boundary_right(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void reflection_boundary_left(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void RT_boundary(Fluid2d* fluids, Block2d block);
// three-dimensional problem
typedef void(*BoundaryCondition3d) (Fluid3d *fluids, Block3d block, Fluid3d bcvalue);
void free_boundary_xleft(Fluid3d *fluids, Block3d block, Fluid3d bcvalue);
void free_boundary_xright(Fluid3d *fluids, Block3d block, Fluid3d bcvalue);
void free_boundary_yleft(Fluid3d *fluids, Block3d block, Fluid3d bcvalue);
void free_boundary_yright(Fluid3d *fluids, Block3d block, Fluid3d bcvalue);
void free_boundary_zleft(Fluid3d *fluids, Block3d block, Fluid3d bcvalue);
void free_boundary_zright(Fluid3d *fluids, Block3d block, Fluid3d bcvalue);

