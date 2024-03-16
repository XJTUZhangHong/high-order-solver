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

// two-dimensioal problem
typedef void(*TimeMarchingCoefficient_2d)(Block2d &block);
extern TimeMarchingCoefficient_2d timecoe_list_2d;

void S1O1_2D(Block2d &block);
void S1O2_2D(Block2d& block);
void RK2_2D (Block2d& block);
void S2O4_2D(Block2d &block);
void RK4_2D(Block2d& block);
void RK3_2D(Block2d& block);
void Initial_stages(Block2d &block);

double Get_CFL(Block2d& block, Fluid2d* fluids, double tstop);

double Dtx(double dtx, double dx, double CFL, double density, double u,double v, double pressure);

void Update(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, int stage);

void Update_with_gauss(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, int stage);

void Update_with_gauss_RK3(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, int stage);

void Update_with_source(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, int stage);

// three-dimensional problem
typedef void(*TimeMarchingCoefficient_3d)(Block3d &block);
extern TimeMarchingCoefficient_3d timecoe_list_3d;

void S1O1_3D(Block3d &block);
void S1O2_3D(Block3d &block);
void S1O3_3D(Block3d &block);
void S2O4_3D(Block3d &block);
void RK2_3D(Block3d &block);
void RK3_3D(Block3d& block);
void RK4_3D(Block3d &block);
void Initial_stages(Block3d &block);

double Get_CFL(Block3d &block, Fluid3d *fluids, double tstop);

double Dtx_3d(double dtx, double dx, double CFL, double density, double u, double v,double w, double pressure);

void Update_gauss(Fluid3d *fluids, Flux3d_gauss** xfluxes, Flux3d_gauss** yfluxes, Flux3d_gauss** zfluxes, Block3d block, int stage);