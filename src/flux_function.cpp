#include "flux_function.h"

GKS1d_type gks1dsolver = nothing; //initialization
Flux_function flux_function = LF; //initialization

void LF(Flux1d& flux, Interface1d& interface, double dt)
{
	double pl[3], pr[3];
	Convar_to_primvar_1D(pl, interface.left.convar);
	Convar_to_primvar_1D(pr, interface.right.convar);

	//Reimann invariants, u +/- 2c/(r-1)
	//Sound speed, c = sqrt(r*p/density)
	double k[2];
	k[0] = abs(pl[1]) + sqrt(Gamma * pl[2] / pl[0]); //abs(u)+c, left
	k[1] = abs(pr[1]) + sqrt(Gamma * pr[2] / pr[0]); //abs(u)+c, right

	double beta = k[0]; // beta, means, abs(partialF/partialX)
	if (k[1] > k[0]) { beta = k[1]; }
	double flux_l[3], flux_r[3];
	get_Euler_flux(pl, flux_l);
	get_Euler_flux(pr, flux_r);

	for (int m = 0; m < 3; m++)
	{
		flux.f[m] = 0.5 * ((flux_l[m] + flux_r[m]) - beta * (interface.right.convar[m] - interface.left.convar[m]));
		flux.f[m] *= dt;
	}
}

void get_Euler_flux(double p[3], double* flux)
{
	//the flux of Euler equation, density*u, density*u*u+p, (p+0.5*density*u*u+p/(r-1))*u
	//p[0 1 2], means, density, u, p
	flux[0] = p[0] * p[1];
	flux[1] = p[0] * p[1] * p[1] + p[2];
	double ENERGS = 0.5 * (p[1] * p[1]) * p[0] + p[2] / (Gamma - 1.0);
	flux[2] = p[1] * (ENERGS + p[2]);
}