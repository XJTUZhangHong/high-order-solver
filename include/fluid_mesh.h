#pragma once
#include<omp.h>
#include<iostream>
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