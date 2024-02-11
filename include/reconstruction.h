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

// interface center reconstruction
void Reconstruction_forg0(Interface1d *interfaces, Fluid1d *fluids, Block1d block);
typedef void (*Reconstruction_forG0)(Interface1d &interfaces, Fluid1d *fluids);
extern Reconstruction_forG0 g0reconstruction;
void Center_collision(Interface1d &interfaces, Fluid1d *fluids);