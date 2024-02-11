#include "basic_function.h"

int K = 0; //initialization  // r=1.4，K=（4-2r）/（r-1）=3
double Gamma = 2.0 / (K + 2) + 1.0; // initialization
double T_inf = 285; // the real temperature, 20 Celcius degree
double R_gas = 8.31441 / (28.959e-3); // initialization // the gas constant for real ideal air
double Mu = -1.0; // mu = rho*u*L/Re
double Nu = -1.0; // Nu = u*L/Re
double c1_euler = -1.0; //initialization
double c2_euler = -1.0; //initialization
bool is_Prandtl_fix = false; //initialization
double Pr = 1.0; //initialization
TAU_TYPE tau_type=Euler; // viscous problem or non-viscous problem

void Convar_to_primvar_1D(double* primvar, double convar[3])
{
	primvar[0] = convar[0];
	primvar[1] = U(convar[0], convar[1]);
	primvar[2] = Pressure(convar[0], convar[1], convar[2]);

}

double Pressure(double density, double densityu, double densityE)
{
	return (Gamma - 1) * (densityE - 0.5 * densityu * densityu / density);
}

double U(double density, double q_densityu)
{
	return q_densityu / density;
}