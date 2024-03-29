#pragma once
#include "flux_function.h"

enum Reconstruction_variable{conservative,characteristic};
extern Reconstruction_variable reconstruction_variable;
enum WENOtype { linear, wenojs, wenoz };
extern WENOtype wenotype;
extern double df_thres;
extern bool is_reduce_order_warning;
extern bool is_using_df_factor;
extern bool smooth;

// one-dimensional problem
void Reconstruction_within_cell(Interface1d *interfaces, Fluid1d *fluids, Block1d block);
typedef void(*Reconstruction_within_Cell)(Point1d &left, Point1d &right, Fluid1d *fluids, Block1d block);
extern Reconstruction_within_Cell cellreconstruction;

void Check_Order_Reduce(Point1d& left, Point1d& right, Fluid1d& fluid);

void Check_Order_Reduce_by_Lambda_1D(bool& order_reduce, double* convar);

// interface left & right reconstruction
void Reconstruction_within_cell(Interface1d* interfaces, Fluid1d* fluids, Block1d block);
double Calculate_alpha_k_1d(double* prim_left, double* prim_right);
void Update_alpha(Interface1d* interfaces, Fluid1d* fluids, Block1d block);
void Vanleer(Point1d& left, Point1d& right, Fluid1d* fluids, Block1d block);
void WENO5_AO(Point1d& left, Point1d& right, Fluid1d* fluids, Block1d block);
void weno_5th_ao_left(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double h);
void weno_5th_ao_right(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double h);
void WENO5_AO_with_DF(Point1d& left, Point1d& right, Fluid1d* fluids, Block1d block);
void weno_5th_ao_with_df_left(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double* df, double h);
void weno_5th_ao_with_df_right(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double* df, double h);
void WENO7_AO_with_DF(Point1d& left, Point1d& right, Fluid1d* fluids, Block1d block);
void weno_7th_ao_with_df_left(double& var, double& der1, double& der2, double wn3, double wn2, double wn1, double w0, double wp1, double wp2, double wp3, double* df, double h);
void weno_7th_ao_with_df_right(double& var, double& der1, double& der2, double wn3, double wn2, double wn1, double w0, double wp1, double wp2, double wp3, double* df, double h);
void WENO9_AO_with_DF(Point1d& left, Point1d& right, Fluid1d* fluids, Block1d block);
void weno_9th_ao_with_df_left(double& var, double& der1, double& der2, double wn4, double wn3, double wn2, double wn1, double w0, double wp1, double wp2, double wp3, double wp4, double* df, double h);
void weno_9th_ao_with_df_right(double& var, double& der1, double& der2, double wn4, double wn3, double wn2, double wn1, double w0, double wp1, double wp2, double wp3, double wp4, double* df, double h);

// interface center reconstruction
void Reconstruction_forg0(Interface1d *interfaces, Fluid1d *fluids, Block1d block);
typedef void (*Reconstruction_forG0)(Interface1d &interfaces, Fluid1d *fluids);
extern Reconstruction_forG0 g0reconstruction;
void Center_collision(Interface1d &interfaces, Fluid1d *fluids);

// two-dimensional problem
double Calculate_alpha_k(double* prim_left, double* prim_right);
void Calculate_alpha(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d& fluids);
void Update_alpha(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block);

void Reconstruction_within_cell(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
typedef void(*Reconstruction_within_Cell_2D_normal)(Interface2d &left, Interface2d &right, Interface2d &down, Interface2d &up, Fluid2d *fluids, Block2d block);
extern Reconstruction_within_Cell_2D_normal cellreconstruction_2D_normal;
void Check_Order_Reduce_by_Lambda_2D(bool& order_reduce, double* convar);
void First_order_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block);
void first_order(Point2d& left, Point2d& right, double* normal_l, double* normal_r, double* w);

void WENO5_AO_normal(Interface2d &left, Interface2d &right, Interface2d &down, Interface2d &up, Fluid2d *fluids, Block2d block);
void WENO5_AO(Point2d &left, Point2d &right, double * wn2, double * wn1, double * w, double * wp1, double * wp2, double h);
void WENO5_AO_with_df_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block);
void WENO5_AO_with_df(Point2d& left, Point2d& right, double* alpha, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double h);
void WENO7_AO_with_df_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block);
void WENO7_AO_with_df(Point2d& left, Point2d& right, double* alpha, double* wn3, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double* wp3, double h);
void WENO9_AO_with_df_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block);
void WENO9_AO_with_df(Point2d& left, Point2d& right, double* alpha, double* wn4, double* wn3, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double* wp3, double* wp4, double h);

typedef void(*Reconstruction_within_Cell_2D_tangent)(Interface2d *left, Interface2d *right, Interface2d *down, Interface2d *up, Fluid2d *fluids, Block2d block);
extern Reconstruction_within_Cell_2D_tangent cellreconstruction_2D_tangent;
void First_order_tangent(Interface2d* left, Interface2d* right, Interface2d* down, Interface2d* up, Fluid2d* fluids, Block2d block);
void first_order_tangent(Point2d& gauss, Point2d& w0);
void WENO5_AO_tangent(Interface2d *left, Interface2d *right, Interface2d *down, Interface2d *up, Fluid2d *fluids, Block2d block);
void WENO5_AO_tangential(Recon2d *re, Recon2d &wn2, Recon2d &wn1, Recon2d &w0, Recon2d &wp1, Recon2d &wp2, double h);
void weno_5th_ao_2gauss(double &g1, double &g1x, double &g1xx, double &g2, double &g2x, double &g2xx, double wn2, double wn1, double w0, double wp1, double wp2, double h, int order);
void WENO5_AO_with_df_tangent(Interface2d* left, Interface2d* right, Interface2d* down, Interface2d* up, Fluid2d* fluids, Block2d block);
void weno_5th_ao_with_df_tangential(Recon2d* re, Recon2d& wn2, Recon2d& wn1, Recon2d& w0, Recon2d& wp1, Recon2d& wp2, double* alpha1, double* alpha2, double h);
void weno_5th_ao_with_df_2gauss(double& g1, double& g1x, double& g1xx, double& g2, double& g2x, double& g2xx, double* df, double wn2, double wn1, double w0, double wp1, double wp2, double h, int order);
void WENO7_AO_with_df_tangent(Interface2d* left, Interface2d* right, Interface2d* down, Interface2d* up, Fluid2d* fluids, Block2d block);
void weno_7th_ao_with_df_tangential(Recon2d* re, Recon2d& wn3, Recon2d& wn2, Recon2d& wn1, Recon2d& w0, Recon2d& wp1, Recon2d& wp2, Recon2d& wp3, double* alpha1, double* alpha2, double h);
void weno_7th_ao_with_df_2gauss(double& g1, double& g1x, double& g2, double& g2x, double& g3, double& g3x, double* df, double wn3, double wn2, double wn1, double w0, double wp1, double wp2, double wp3, double h, int order);
void WENO9_AO_with_df_tangent(Interface2d* left, Interface2d* right, Interface2d* down, Interface2d* up, Fluid2d* fluids, Block2d block);
void weno_9th_ao_with_df_tangential(Recon2d* re, Recon2d& wn4, Recon2d& wn3, Recon2d& wn2, Recon2d& wn1, Recon2d& w0, Recon2d& wp1, Recon2d& wp2, Recon2d& wp3, Recon2d& wp4, double* alpha1, double* alpha2, double h);
void weno_9th_ao_with_df_2gauss(double& g1, double& g1x, double& g2, double& g2x, double& g3, double& g3x, double& g4, double& g4x, double* df, double wn4, double wn3, double wn2, double wn1, double w0, double wp1, double wp2, double wp3, double wp4, double h, int order);
void Polynomial_9th(double* p, double* px, double* df, double wn4, double wn3, double wn2, double wn1, double w0, double wp1, double wp2, double wp3, double wp4, double x, double h);

void Reconstruction_forg0(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
typedef void(*Reconstruction_forG0_2D_normal)(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
extern Reconstruction_forG0_2D_normal g0reconstruction_2D_normal;
void Center_do_nothing_normal(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block);

typedef void(*Reconstruction_forG0_2D_tangent)(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
extern Reconstruction_forG0_2D_tangent g0reconstruction_2D_tangent;
void Center_all_collision_multi(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
void Center_all_collision_2d_multi(Recon2d &gauss);

// three-dimensional problem
enum G0_construct_type { collisionn, collisionnless,all_collisionn };
extern G0_construct_type g0type;

void Reconstruction_within_cell(Interface3d *xinterfaces, Interface3d *yinterfaces, Interface3d *zinterfaces, Fluid3d *fluids, Block3d block);

typedef void(*Reconstruction_within_Cell_3D_of_face_value)
(Interface3d &left, Interface3d &right, Interface3d &down, Interface3d &up, Interface3d &back, Interface3d &front, Fluid3d *fluids, Block3d block);
extern Reconstruction_within_Cell_3D_of_face_value cellreconstruction_3D_of_face_value;
void WENO5_AO_splitting_3d(Interface3d &left, Interface3d &right, Interface3d &down, Interface3d &up, Interface3d &back, Interface3d &front, Fluid3d *fluids, Block3d block);
void WENO5_AO(Point3d &left, Point3d &right, double * wn2, double * wn1, double * w, double * wp1, double * wp2, double h);

typedef void(*Reconstruction_within_Cell_3D_of_line_value)
(Interface3d *xinterface, Interface3d *yinterface, Interface3d *zinterface, Block3d block);
extern Reconstruction_within_Cell_3D_of_line_value cellreconstruction_3D_of_line_value;
void WENO5_AO_of_line_value(Interface3d *right, Interface3d *up, Interface3d *front, Block3d block);
void WENO5_AO_for_line_value(Recon3d *line, Recon3d& wn2, Recon3d& wn1, Recon3d& w0, Recon3d& wp1, Recon3d& wp2, double h);

typedef void(*Reconstruction_within_Cell_3D_of_point_value)
(Interface3d *xinterface, Interface3d *yinterface, Interface3d *zinterface, Block3d block);
extern Reconstruction_within_Cell_3D_of_point_value cellreconstruction_3D_of_point_value;
void WENO5_AO_of_point_value(Interface3d *right, Interface3d *up, Interface3d *front, Block3d block);
void WENO5_AO_for_point_value(Recon3d *point, Recon3d& wn2, Recon3d& wn1, Recon3d& w0, Recon3d& wp1, Recon3d& wp2, double h);

void Reconstruction_forg0(Interface3d *xinterfaces, Interface3d *yinterfaces, Interface3d *zinterfaces, Fluid3d *fluids, Block3d block);

typedef void(*Reconstruction_forG0_3D_of_face_value)(Interface3d *xinterfaces, Interface3d *yinterfaces, Interface3d *zinterfaces, Fluid3d *fluids, Block3d block);
extern Reconstruction_forG0_3D_of_face_value g0reconstruction_3D_of_face_value;
void Do_nothing_splitting_3d(Interface3d *xinterfaces, Interface3d *yinterfaces, Interface3d *zinterfaces, 
	Fluid3d *fluids, Block3d block);

typedef void(*Reconstruction_forG0_3D_of_line_value)
(Interface3d *xinterfaces, Interface3d *yinterfaces, Interface3d *zinterfaces, Block3d block);
extern Reconstruction_forG0_3D_of_line_value g0reconstruction_3D_of_line_value;
void Do_nothing_reconstruction_of_line_value
	(Interface3d *xinterfaces, Interface3d *yinterfaces, Interface3d *zinterfaces, Block3d block);

typedef void(*Reconstruction_forG0_3D_of_point_value)
(Interface3d *xinterfaces, Interface3d *yinterfaces, Interface3d *zinterfaces, Block3d block);
extern Reconstruction_forG0_3D_of_point_value g0reconstruction_3D_of_point_value;
void Do_nothing_reconstruction_of_point_value
	(Interface3d *xinterfaces, Interface3d *yinterfaces, Interface3d *zinterfaces, Block3d block);