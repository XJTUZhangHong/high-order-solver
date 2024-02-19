#pragma once
#include <omp.h>
#include <cmath>
#include <iostream>
using namespace std;

// one dimensional problem
class Block1d
{
public:
	//bool uniform;
	int ghost;
	int nodex; // mesh number
	int nx;
	int nodex_begin;
	int nodex_end;
	double dx;
	double left;
	double right;
	int stages;
	double timecoefficient[5][5][2];
	double t; //current simulation time
	double CFL; //cfl number, actually just a amplitude factor
	double dt;//current_global time size
	int step; //start from which step
};

// to remeber the cell avg values
class Fluid1d
{
public:
	double primvar[3];
	double convar[3];
	double convar_old[3];
	double convar_stage1[3];
	double convar_stage2[3];
	double cx; //center coordinate in x direction
	double dx; //the mesh size dx
	double alpha = 1.0;
};

// to remember the fluid information in a fixed point,
// such as reconstructioned value, or middle point value
class Point1d
{
public:
	double convar[3];
	double convar_old[3]; //this one is for compact point-wise reconstruction
	double der1[3];
	double x; // coordinate of a point
};

// remember the flux in a fixed interface,
// for RK method, we only need f
// for 2nd der method, we need derf
class Flux1d
{
public:
	double F[3]; //total flux in dt time
	double f[3]; // the f0 in t=0
	double derf[4]; // the f_t in t=0
};

// every interfaces have left center and right value.
// and center value for GKS
// and a group of flux for RK method or multi-stage GKS method
class Interface1d
{
public:
	Point1d left;
	Point1d center;
	Point1d right;
	Flux1d* flux;
	double x; // coordinate of the interface, equal to point1d.x
};

// two-dimensional problem
extern int gausspoint;
extern double *gauss_loc;
extern double *gauss_weight;

class Block2d
{
	//geometry information of 2D block
public:
	int index;
	bool uniform;
	int ghost;
	int nodex;
	int nodey;
	int nx; //total = real + ghost
	int ny; //total = real + ghost
	int xcell_begin;
	int xcell_end;
	int ycell_begin;
	int ycell_end;
	int xinterface_begin_n;
	int xinterface_begin_t;
	int xinterface_end_n;
	int xinterface_end_t;
	int yinterface_begin_n;
	int yinterface_begin_t;
	int yinterface_end_n;
	int yinterface_end_t;
	double dx;
	double dy;
	double overdx;
	double overdy;
	double left;
	double right;
	double down;
	double up;
	int stages;
	double timecoefficient[5][5][2];
	double t; //current simulation time
	double CFL; //cfl number, actually just a amplitude factor
	double dt;//current_global time size
	int step; //start from which step
	double T; //stop time
	int outputperstep;
};

// to remember the fluid information in a fixed point,
// such as reconstructioned value, or middle point value
class Point2d
{
	// Note : store values relating to one point
public:
	double convar[4];
	double primvar[4];
	double der1x[4];
	double der1y[4];
	double x;
	double y;
	double normal[2];
	bool is_reduce_order;
};

// to remeber the cell avg values
//typedef class Fluid2d
class Fluid2d
{
	// Note : store the cell-averaged values
	// do not have the global geometry/block information
public:
	double primvar[4];
	double exact;
	double convar[4];
	double convar_old[4];
	double convar_stage1[4];
	double convar_stage2[2];
	double res[4];
	double xindex;
	double yindex;
	double coordx;
	double coordy;
	double dx;
	double dy;
	double area;
	bool is_hweno;
	double node[8];
	bool boundary;
	double alpha_x = 1.0;
	double alpha_y = 1.0;
	string direction;   // x direction or y direction
	bool is_left = false;
	double Sw1[2][4];
	double Lw1[2][4];
	double derSw1[2][4];
	double derLw1[2][4];
};

// remember the flux in a fixed interface,
// for RK method, we only need f
// for 2nd der method, we need derf
class Flux2d
{
public:
	double x[4];
	double f[4];
	double derf[4];
	double length;
	double normal[2];
};

// every interfaces have left center and right value.
// and center value for GKS
// and a group of flux for RK method or multi-stage GKS method
class Recon2d
{
	// Note : reconstruction can obtain left, cenrer, and right value
public:
	Point2d left;
	Point2d center;
	Point2d right;
	double x;
	double y;
	double normal[2];
};

class Interface2d
{
	// Note : for one interface, reconstruction can be for line-averagend or gauss point
	// These reconstructed values can be used for the calculation of the flux at gauss points (stored in Flux2d_gauss) of the interface
public:

	Recon2d line;
	Recon2d* gauss;
	double x;
	double y;
	double normal[2];
	double length;
};

class Flux2d_gauss
{
	// Note : Flux2d store the variables relating to flux
	// Flux2d_gauss, includes several gauss points
	// Final, one interface should have one Flux2d_gauss
public:
	Flux2d* gauss;
};

void SetGuassPoint();

// three-dimensional problem
extern double **gauss_loc_3d;
extern bool is_multi_dimensional;

class Block3d
{
public:
	bool uniform;
	int ghost;
	int nodex; int nodey; int nodez;
	int nx; int ny; int nz;
	double dx; double dy; double dz;
	double overdx; double overdy; double overdz;
	double xarea; double yarea; double zarea;
	double left; double right;
	double down; double up;
	double front; double back;
	int stages;
	double timecoefficient[5][5][3];
	double timecoefficient_hweno[5][5][3];
	double t; //current simulation time
	double CFL; //cfl number, actually just a amplitude factor
	double dt;//current_global time size
	double tstop;//stop time
	int step; //start from which step
	double xcell_begin; double xcell_end; double ycell_begin; double ycell_end; double zcell_begin; double zcell_end;
	double xinterface_begin_n; double xinterface_end_n; double xinterface_begin_t; double xinterface_end_t;
	double yinterface_begin_n; double yinterface_end_n; double yinterface_begin_t; double yinterface_end_t;
	double zinterface_begin_n; double zinterface_end_n; double zinterface_begin_t; double zinterface_end_t;
	int index;
	int outputstep_fluid;
	int outputstep_res;
	int outputstep_wall_info;
	int outputstep_continue_file;
	int output_time_step;
	int outputstep_extract_point_data;
};

class Fluid3d
{
public:
	double primvar[6];
	double convar[6];
	double convar_old[6];
	int index[3];
	double loc[3];
	double base[10]; //high-order reconstruction base1 2 3 is location for uns mesh
	double *geo;
	double node[24];
	double vol;
	double avg_derx[6];
	double avg_dery[6];
	double avg_derz[6];
	double* beta;
	double R_c;
	double* dt;
	double* Res;
	double* C;
	double* DW;
	double slope_compress;
	double vec_lambda[2];
	double *weno_weight;
	double wall_dis;
	double* source;
};

class Point3d
{
public:
	double convar[6];
	double der1x[6];
	double der1y[6];
	double der1z[6];
};

class Recon3d
{
public:
	Point3d left;
	Point3d center;
	Point3d right;
	bool order_reduce;
};

class Recon3d_lr
{
public:
	Point3d left;
	Point3d right;
	bool order_reduce;
};

class Interface3d
{
public:

	Recon3d face;
	Recon3d *line;
	Recon3d *gauss;
	double loc[12];
	double weight[4];
	double normal[12];
	int index[3];
	double length[2];
	double area;
	int bc_type;
};

class Flux3d
{
public:
	double x[6];
	double f[6];
	double derf[6];
	double d[30];
	double w[30];
	double derw[30];
};

void SetGuassPoint_3D();