#pragma once
#include "flux_function.h"

enum Reconstruction_variable{conservative,characteristic};
extern Reconstruction_variable reconstruction_variable;
enum WENOtype { linear, wenojs, wenoz };
extern WENOtype wenotype;
extern bool is_reduce_order_warning;
extern bool is_using_df_factor;

// one-dimensional problem
void Reconstruction_within_cell(Interface1d *interfaces, Fluid1d *fluids, Block1d block);
typedef void(*Reconstruction_within_Cell)(Point1d &left, Point1d &right, Fluid1d *fluids, Block1d block);
extern Reconstruction_within_Cell cellreconstruction;

void Check_Order_Reduce(Point1d& left, Point1d& right, Fluid1d& fluid);

void Check_Order_Reduce_by_Lambda_1D(bool& order_reduce, double* convar);

// interface left & right reconstruction
void Reconstruction_within_cell(Interface1d* interfaces, Fluid1d* fluids, Block1d block);

void Vanleer(Point1d& left, Point1d& right, Fluid1d* fluids, Block1d block);

void WENO5_AO(Point1d& left, Point1d& right, Fluid1d* fluids, Block1d block);

void weno_5th_ao_left(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double h);

void weno_5th_ao_right(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double h);

double Calculate_alpha_k_1d(double* prim_left, double* prim_right);

void Update_alpha(Interface1d* interfaces, Fluid1d* fluids, Block1d block);

void WENO5_AO_with_DF(Point1d& left, Point1d& right, Fluid1d* fluids, Block1d block);

void weno_5th_ao_with_df_left(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double* df, double h);

void weno_5th_ao_with_df_right(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double* df, double h);

// interface center reconstruction
void Reconstruction_forg0(Interface1d *interfaces, Fluid1d *fluids, Block1d block);
typedef void (*Reconstruction_forG0)(Interface1d &interfaces, Fluid1d *fluids);
extern Reconstruction_forG0 g0reconstruction;
void Center_collision(Interface1d &interfaces, Fluid1d *fluids);


// two-dimensional problem
void Reconstruction_within_cell(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
typedef void(*Reconstruction_within_Cell_2D_normal)(Interface2d &left, Interface2d &right, Interface2d &down, Interface2d &up, Fluid2d *fluids, Block2d block);
extern Reconstruction_within_Cell_2D_normal cellreconstruction_2D_normal;
void Check_Order_Reduce_by_Lambda_2D(bool& order_reduce, double* convar);
void WENO5_AO_normal(Interface2d &left, Interface2d &right, Interface2d &down, Interface2d &up, Fluid2d *fluids, Block2d block);
void WENO5_AO(Point2d &left, Point2d &right, double * wn2, double * wn1, double * w, double * wp1, double * wp2, double h);

typedef void(*Reconstruction_within_Cell_2D_tangent)(Interface2d *left, Interface2d *right, Interface2d *down, Interface2d *up, Fluid2d *fluids, Block2d block);
extern Reconstruction_within_Cell_2D_tangent cellreconstruction_2D_tangent;
void WENO5_AO_tangent(Interface2d *left, Interface2d *right, Interface2d *down, Interface2d *up, Fluid2d *fluids, Block2d block);
void WENO5_AO_tangential(Recon2d *re, Recon2d &wn2, Recon2d &wn1, Recon2d &w0, Recon2d &wp1, Recon2d &wp2, double h);
void weno_5th_ao_2gauss(double &g1, double &g1x, double &g1xx, double &g2, double &g2x, double &g2xx, double wn2, double wn1, double w0, double wp1, double wp2, double h, int order);

void Reconstruction_forg0(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
typedef void(*Reconstruction_forG0_2D_normal)(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
extern Reconstruction_forG0_2D_normal g0reconstruction_2D_normal;
void Center_do_nothing_normal(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block);

typedef void(*Reconstruction_forG0_2D_tangent)(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
extern Reconstruction_forG0_2D_tangent g0reconstruction_2D_tangent;
void Center_all_collision_multi(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
void Center_all_collision_2d_multi(Recon2d &gauss);