#pragma once
#include "time_advance.h"
#include "reconstruction.h"

// one-dimensional problem
enum GKS1d_type{nothing, kfvs1st, kfvs2nd, gks1st, gks2nd, pp_gks}; // positive-preserving gks
extern GKS1d_type gks1dsolver;

void Calculate_flux(Flux1d** fluxes, Interface1d* interfaces, Block1d &block, int stage);
typedef void(*Flux_function)(Flux1d &flux, Interface1d& interface, double dt);
extern Flux_function flux_function;
void GKS(Flux1d& flux, Interface1d& interface, double dt);
void LF(Flux1d& flux, Interface1d& interface, double dt);
void get_Euler_flux(double p[3], double* flux);

// two-dimensional problem
enum GKS2d_type { nothing_2d, kfvs1st_2d, kfvs2nd_2d, gks1st_2d, gks2nd_2d, pp_gks_2d};
extern GKS2d_type gks2dsolver;

void Calculate_flux(Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Interface2d* xinterfaces, Interface2d* yinterfaces, Block2d block, int stage);
typedef void(*Flux_function_2d)(Flux2d &flux, Recon2d & re, double dt);
extern Flux_function_2d flux_function_2d;
void GKS2D_smooth(Flux2d &flux, Recon2d& interface, double dt);
void GKS2D(Flux2d &flux, Recon2d& interface, double dt);
void LF2D(Flux2d& flux, Recon2d& interface, double dt);
void get_flux(double p[4], double* flux);
void NS_by_central_difference_prim_2D(Flux2d& flux, Recon2d& interface, double dt);

// three-dimensional problem
enum GKS3d_type {
	nothing_3d, kfvs1st_3d, kfvs2nd_3d, gks1st_3d, gks2nd_3d
};
extern GKS3d_type gks3dsolver;

typedef void(*Flux_function_3d)(Flux3d &flux, Recon3d & re, double area, double dt);
extern Flux_function_3d flux_function_3d;
void Calculate_flux(Flux3d_gauss** xfluxes, Flux3d_gauss** yfluxes, Flux3d_gauss** zfluxes,
	Interface3d* xinterfaces, Interface3d* yinterfaces, Interface3d* zinterfaces,
	Block3d block, int stage);
void GKS3D(Flux3d &flux, Recon3d& interface, double area, double dt);
void GKS3D_speed(Flux3d &flux, Recon3d& interface, double area, double dt);