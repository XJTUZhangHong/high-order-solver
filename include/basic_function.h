#pragma once
#include "fluid_mesh.h"
#include <assert.h>

extern int K;
extern double Gamma;
extern double Mu;
extern double Nu;
extern double T_inf; //the real temperature
extern double c1_euler;
extern double c2_euler;
extern bool is_Prandtl_fix;
extern double Pr;
extern double R_gas; //the gas constant for real ideal air
enum TAU_TYPE { Euler, NS};
extern TAU_TYPE tau_type;

//basic gks function
// to store the moment
class MMDF1d
{
private:
	double u;
	double lambda;

public:
	double uwhole[10];
	double uplus[10];
	double uminus[10];
	double upxi[10][4];
	double unxi[10][4];
	double uxi[10][4];
	double xi2;
	double xi4;
	double xi6;
	MMDF1d();
	MMDF1d(double u_in, double lambda_in);
	void calcualte_MMDF1d();
};

class MMDF
{
private:
	double u;
	double v;
	double lambda;

public:
	double uwhole[7];
	double uplus[7];
	double uminus[7];
	double vwhole[7];
	double upvxi[7][7][3];
	double unvxi[7][7][3];
	double uvxi[7][7][3];
	double xi2;
	double xi4;
	MMDF();
	MMDF(double u_in, double v_in, double lambda_in);
	void calcualte_MMDF();
};

// to calculate the microsolpe moment
void G(int no_u, int no_xi, double* psi, double a[3], MMDF1d m);
void GL(int no_u, int no_xi, double* psi, double a[3], MMDF1d m);
void GR(int no_u, int no_xi, double* psi, double a[3], MMDF1d m);

void Microslope(double* a, double der[3], double prim[3]);

void Collision(double* w0, double left, double right, MMDF& m2, MMDF& m3);

void A(double* a, double der[4], double prim[4]);

void GL_address(int no_u, int no_v, int no_xi, double* psi, double a[4], MMDF& m);

void GR_address(int no_u, int no_v, int no_xi, double* psi, double a[4], MMDF& m);

void G_address(int no_u, int no_v, int no_xi, double* psi, double a[4], MMDF& m);

double Get_Tau_NS(double density0, double lambda0);

double Get_Tau(double density_left, double density_right, double density0, double lambda_left, double lambda_right, double lambda0, double dt);

double Alpha(double lambda, double u);

double Beta(double lambda, double u);

void Global_to_Local(double* change, double* origin, double* normal);

void Array_zero(double* target, int dim);

bool negative_density_or_pressure(double* primvar);

void Convar_to_Primvar(Fluid2d* fluids, Block2d block);

void Convar_to_primvar_1D(double* primvar, double convar[3]);

void Convar_to_primvar_2D(double *primvar, double convar[4]);

void Convar_to_ULambda_1d(double* primvar, double convar[3]);

void Convar_to_ULambda_2d(double* primvar, double convar[4]);

void Primvar_to_convar_1D(double* convar, double primvar[3]);

void Primvar_to_convar_2D(double *convar, double primvar[4]);

void Convar_to_char1D(double* character, double primvar[3], double convar[3]);

void Char_to_convar1D(double* convar, double primvar[3], double charvar[3]);

void Convar_to_char(double *character, double *primvar, double convar[4]);

void Char_to_convar(double *convar, double *primvar, double character[4]);

double DensityU(double density, double u);

double DensityE(double density, double u, double pressure);

double Pressure(double density, double densityu, double densityE);

double Q_densityu(double density, double u);

double Q_densityv(double density, double v);

double Q_densityE(double density, double u, double v, double pressure);

double U(double density, double q_densityu);

double V(double density, double q_densityv);

double Pressure(double density, double q_densityu, double q_densityv, double q_densityE);

double Lambda(double density, double u, double densityE);

double Lambda(double density, double u, double v, double densityE);

void Copy_Array(double* target, double* origin, int dim);

void YchangetoX(double* fluidtmp, double* fluid);

void XchangetoY(double* fluidtmp, double* fluid);

void Local_to_Global(double *change,double *normal);

Flux1d** Setflux_array(Block1d block);

void SetUniformMesh(Block1d block, Fluid1d* fluids, Interface1d* interfaces, Flux1d** fluxes);

Fluid2d* Setfluid(Block2d& block);

Interface2d* Setinterface_array(Block2d block);

Flux2d_gauss** Setflux_gauss_array(Block2d block);

void SetUniformMesh(Block2d& block, Fluid2d* fluids, Interface2d* xinterfaces, Interface2d* yinterfaces, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes);

void Copy_geo_from_interface_to_line(Interface2d& interface);

void Copy_geo_from_interface_to_flux(Interface2d& interface, Flux2d_gauss* flux, int stages);

void Set_Gauss_Coordinate(Interface2d& xinterface, Interface2d& yinterface);

void Set_Gauss_Intergation_Location_x(Point2d& xgauss, int index, double h);

void Set_Gauss_Intergation_Location_y(Point2d& ygauss, int index, double h);

void CopyFluid_new_to_old(Fluid1d* fluids, Block1d block);

void CopyFluid_new_to_old(Fluid2d* fluids, Block2d block);