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
	double alpha = 1.0; // compression factor
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