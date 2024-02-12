#pragma once
#include "basic_function.h"

enum Reconstruction_variable{conservative,characteristic};
extern Reconstruction_variable reconstruction_variable;
enum WENOtype { linear, wenojs, wenoz };
extern WENOtype wenotype;
extern bool is_reduce_order_warning;

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

void weno_5th_ao_left(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double* beta, double h);

void weno_5th_ao_right(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double* beta, double h);

// interface center reconstruction
void Reconstruction_forg0(Interface1d *interfaces, Fluid1d *fluids, Block1d block);
typedef void (*Reconstruction_forG0)(Interface1d &interfaces, Fluid1d *fluids);
extern Reconstruction_forG0 g0reconstruction;
void Center_collision(Interface1d &interfaces, Fluid1d *fluids);