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

// to prepare the basic element for moment calculation
MMDF1d::MMDF1d() { u = 1.0; lambda = 1.0; };
MMDF1d::MMDF1d(double u_in, double lambda_in)
{
	u = u_in;
	lambda = lambda_in;
	calcualte_MMDF1d();
}
void MMDF1d::calcualte_MMDF1d()
{
	uwhole[0] = 1;
	uwhole[1] = u;
	uplus[0] = 0.5 * Alpha(lambda, -u);
	uminus[0] = 0.5 * Alpha(lambda, u);
	uplus[1] = u * uplus[0] + 0.5 * Beta(lambda, u);
	uminus[1] = u * uminus[0] - 0.5 * Beta(lambda, u);
	for (int i = 2; i <= 9; i++)
	{
		uwhole[i] = u * uwhole[i - 1] + 0.5 * (i - 1) / lambda * uwhole[i - 2];
	}
	for (int i = 2; i <= 9; i++)
	{
		uplus[i] = u * uplus[i - 1] + 0.5 * (i - 1) / lambda * uplus[i - 2];
		uminus[i] = u * uminus[i - 1] + 0.5 * (i - 1) / lambda * uminus[i - 2];
	}
	xi2 = 0.5 * K / lambda;
	xi4 = 0.25 * (K * K + 2 * K) / (lambda * lambda);
	//xi6?? how to calculate
	xi6 = 0.5 * (K + 4) / lambda * xi4;

	for (int i = 0; i < 10; i++)
	{
		for (int k = 0; k < 4; k++)
		{
			if ((i + 2 * k) <= 9)
			{
				if (k == 0)
				{
					upxi[i][k] = uplus[i];
					unxi[i][k] = uminus[i];
				}
				if (k == 1)
				{
					upxi[i][k] = uplus[i] * xi2;
					unxi[i][k] = uminus[i] * xi2;
				}
				if (k == 2)
				{
					upxi[i][k] = uplus[i] * xi4;
					unxi[i][k] = uminus[i] * xi4;
				}
				if (k == 3)
				{
					upxi[i][k] = uplus[i] * xi6;
					unxi[i][k] = uminus[i] * xi6;
				}
			}
		}
	}
	for (int i = 0; i < 10; i++)
	{
		for (int k = 0; k < 4; k++)
		{
			if ((i + 2 * k) <= 9)
			{
				if (k == 0)
				{
					uxi[i][k] = uwhole[i];
				}
				if (k == 1)
				{
					uxi[i][k] = uwhole[i] * xi2;
				}
				if (k == 2)
				{
					uxi[i][k] = uwhole[i] * xi4;
				}
				if (k == 3)
				{
					uxi[i][k] = uwhole[i] * xi6;
				}
			}
		}
	}
}

//a general G function
void G(int no_u, int no_xi, double* psi, double a[3], MMDF1d m)
{
	psi[0] = a[0] * m.uxi[no_u][no_xi] + a[1] * m.uxi[no_u + 1][no_xi] + a[2] * 0.5 * (m.uxi[no_u + 2][no_xi] + m.uxi[no_u][no_xi + 1]);
	psi[1] = a[0] * m.uxi[no_u + 1][no_xi] + a[1] * m.uxi[no_u + 2][no_xi] + a[2] * 0.5 * (m.uxi[no_u + 3][no_xi] + m.uxi[no_u + 1][no_xi + 1]);
	psi[2] = 0.5 * (a[0] * (m.uxi[no_u + 2][no_xi] + m.uxi[no_u][no_xi + 1]) +
		a[1] * (m.uxi[no_u + 3][no_xi] + m.uxi[no_u + 1][no_xi + 1]) +
		a[2] * 0.5 * (m.uxi[no_u + 4][no_xi] + m.uxi[no_u][no_xi + 2] + 2 * m.uxi[no_u + 2][no_xi + 1]));

}

void GL(int no_u, int no_xi, double* psi, double a[3], MMDF1d m)
{
	psi[0] = a[0] * m.upxi[no_u][no_xi] + a[1] * m.upxi[no_u + 1][no_xi] + a[2] * 0.5 * (m.upxi[no_u + 2][no_xi] + m.upxi[no_u][no_xi + 1]);
	psi[1] = a[0] * m.upxi[no_u + 1][no_xi] + a[1] * m.upxi[no_u + 2][no_xi] + a[2] * 0.5 * (m.upxi[no_u + 3][no_xi] + m.upxi[no_u + 1][no_xi + 1]);
	psi[2] = 0.5 * (a[0] * (m.upxi[no_u + 2][no_xi] + m.upxi[no_u][no_xi + 1]) +
		a[1] * (m.upxi[no_u + 3][no_xi] + m.upxi[no_u + 1][no_xi + 1]) +
		a[2] * 0.5 * (m.upxi[no_u + 4][no_xi] + m.upxi[no_u][no_xi + 2] + 2 * m.upxi[no_u + 2][no_xi + 1]));
}

void GR(int no_u, int no_xi, double* psi, double a[3], MMDF1d m)
{
	psi[0] = a[0] * m.unxi[no_u][no_xi] + a[1] * m.unxi[no_u + 1][no_xi] + a[2] * 0.5 * (m.unxi[no_u + 2][no_xi] + m.unxi[no_u][no_xi + 1]);
	psi[1] = a[0] * m.unxi[no_u + 1][no_xi] + a[1] * m.unxi[no_u + 2][no_xi] + a[2] * 0.5 * (m.unxi[no_u + 3][no_xi] + m.unxi[no_u + 1][no_xi + 1]);
	psi[2] = 0.5 * (a[0] * (m.unxi[no_u + 2][no_xi] + m.unxi[no_u][no_xi + 1]) +
		a[1] * (m.unxi[no_u + 3][no_xi] + m.unxi[no_u + 1][no_xi + 1]) +
		a[2] * 0.5 * (m.unxi[no_u + 4][no_xi] + m.unxi[no_u][no_xi + 2] + 2 * m.unxi[no_u + 2][no_xi + 1]));
}

//solution of matrix equation b=Ma 1D
void Microslope(double* a, double der[3], double prim[3])
{
	double R4, R2;
	R4 = der[2] / prim[0] - 0.5 * (prim[1] * prim[1] + 0.5 * (K + 1) / prim[2]) * der[0] / prim[0];
	R2 = (der[1] - prim[1] * der[0]) / prim[0];
	a[2] = 4 * prim[2] * prim[2] / (K + 1) * (2 * R4 - 2 * prim[1] * R2);
	a[1] = 2 * prim[2] * R2 - prim[1] * a[2];
	a[0] = der[0] / prim[0] - prim[1] * a[1] - 0.5 * a[2] * (prim[1] * prim[1] + 0.5 * (K + 1) / prim[2]);
}

double Alpha(double lambda, double u)
{
	return erfc(sqrt(lambda) * u);
}

double Beta(double lambda, double u)
{
	double pi = 3.14159265358979323846;
	return exp(-lambda * u * u) / sqrt(pi * lambda);
}

void Convar_to_primvar_1D(double* primvar, double convar[3])
{
	primvar[0] = convar[0];
	primvar[1] = U(convar[0], convar[1]);
	primvar[2] = Pressure(convar[0], convar[1], convar[2]);

}

void Convar_to_ULambda_1d(double* primvar, double convar[3])
{
	primvar[0] = convar[0];
	primvar[1] = U(convar[0], convar[1]);
	primvar[2] = Lambda(convar[0], primvar[1], convar[2]);
}

void Primvar_to_convar_1D(double* convar, double primvar[3])
{
	convar[0] = primvar[0];
	convar[1] = DensityU(primvar[0], primvar[1]);
	convar[2] = DensityE(primvar[0], primvar[1], primvar[2]);
}
double DensityU(double density, double u)
{
	return density * u;
}
double DensityE(double density, double u, double pressure)
{
	return density * (pressure / density / (Gamma - 1) + 0.5 * (u * u));
}

double Pressure(double density, double densityu, double densityE)
{
	return (Gamma - 1) * (densityE - 0.5 * densityu * densityu / density);
}

double U(double density, double q_densityu)
{
	return q_densityu / density;
}

double Lambda(double density, double u, double densityE)
{
	return (K + 1.0) * 0.25 * (density / (densityE - 0.5 * density * (u * u)));
}

Flux1d** Setflux_array(Block1d block)
{
	Flux1d** fluxes = new Flux1d * [block.nx + 1];  // dynamic variable (since block.nx is not determined)

	for (int i = 0; i <= block.nx; i++)
	{
		// for m th step time marching schemes, m subflux needed
		fluxes[i] = new Flux1d[block.stages];
		for (int j = 0; j < block.stages; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				fluxes[i][j].f[k] = 0.0;
				fluxes[i][j].derf[k] = 0.0;
			}
		}
	}
	if (fluxes == 0)
	{
		cout << "memory allocation failed for muli-stage flux";
		return NULL;
	}
	cout << "the memory for muli-stage flux has been allocated..." << endl;
	return fluxes;
}

void SetUniformMesh(Block1d block, Fluid1d* fluids, Interface1d* interfaces, Flux1d** fluxes)
{
	//cell avg information
	for (int i = 0; i < block.nx; i++)
	{
		fluids[i].dx = block.dx; //cell size
		fluids[i].cx = block.left + (i + 0.5 - block.ghost) * block.dx; //cell center location
	}
	// interface information
	for (int i = 0; i <= block.nx; i++)
	{
		interfaces[i].x = block.left + (i - block.ghost) * block.dx;
		interfaces[i].left.x = interfaces[i].x;
		interfaces[i].right.x = interfaces[i].x;
		interfaces[i].center.x = interfaces[i].x;
		interfaces[i].flux = fluxes[i];
	}
}

void CopyFluid_new_to_old(Fluid1d* fluids, Block1d block)
{
#pragma omp parallel  for
	for (int i = block.ghost; i < block.ghost + block.nodex; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			fluids[i].convar_old[j] = fluids[i].convar[j];
		}
	}
}
