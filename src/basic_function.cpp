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

//moments of the maxwellian distribution function
MMDF::MMDF() { u = 1.0; v = 1.0; lambda = 1.0; };
MMDF::MMDF(double u_in, double v_in, double lambda_in)
{
	u = u_in;
	v = v_in;
	lambda = lambda_in;
	calcualte_MMDF();
}
void MMDF::calcualte_MMDF()
{
	//the moment of u in whole domain
	uwhole[0] = 1;
	uwhole[1] = u;
	//the moment of u in half domain (>0 && <0)
	uplus[0] = 0.5 * Alpha(lambda, -u);
	uminus[0] = 1.0 - uplus[0];
	uplus[1] = u * uplus[0] + 0.5 * Beta(lambda, u);
	uminus[1] = u - uplus[1];
	double overlambda = 1.0 / lambda;
	for (int i = 2; i <= 6; i++)
	{
		uwhole[i] = u * uwhole[i - 1] + 0.5 * (i - 1) * overlambda * uwhole[i - 2];
	}
	for (int i = 2; i <= 6; i++)
	{
		uplus[i] = u * uplus[i - 1] + 0.5 * (i - 1) * overlambda * uplus[i - 2];
		uminus[i] = uwhole[i] - uplus[i];
	}
	//the moment of v in whole domain
	//the moment of v in half domain is useless
	vwhole[0] = 1;
	vwhole[1] = v;
	for (int i = 2; i <= 6; i++)
	{
		vwhole[i] = v * vwhole[i - 1] + 0.5 * (i - 1) * overlambda * vwhole[i - 2];
	}
	//the moment of xi (kesi)
	xi2 = 0.5 * K * overlambda;
	xi4 = 0.25 * (K * K + 2 * K) * overlambda * overlambda;

	//xi6 = 0.5*(K + 4) / lambda*xi4; // no use

	//get other variables by the above mement variables
	//i, j, k, means i power of u, j power of v, and 2k power of xi
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 7; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				if ((i + j + 2 * k) <= 6)
				{
					if (k == 0)
					{
						// u, half domain > 0; v, whole domain; xi, equals 0
						upvxi[i][j][k] = uplus[i] * vwhole[j];
						// u, half domain < 0; v, whole domain; xi, equals 0
						unvxi[i][j][k] = uminus[i] * vwhole[j];
						// u, whole domain ; v, whole domain; xi, equals 0
						uvxi[i][j][k] = upvxi[i][j][k] + unvxi[i][j][k];
						// = (uplus[i] + uminus[i]) * vwhole[j] = uwhole[i] * vwhole[j]
					}
					if (k == 1)
					{
						// u, half domain > 0; v, whole domain; xi, equals 2
						upvxi[i][j][k] = uplus[i] * vwhole[j] * xi2;
						// u, half domain < 0; v, whole domain; xi, equals 2
						unvxi[i][j][k] = uminus[i] * vwhole[j] * xi2;
						// u, whole domain ; v, whole domain; xi, equals 2
						uvxi[i][j][k] = upvxi[i][j][k] + unvxi[i][j][k];
					}
					if (k == 2)
					{
						// u, half domain > 0; v, whole domain; xi, equals 4
						upvxi[i][j][k] = uplus[i] * vwhole[j] * xi4;
						// u, half domain < 0; v, whole domain; xi, equals 4
						unvxi[i][j][k] = uminus[i] * vwhole[j] * xi4;
						// u, whole domain ; v, whole domain; xi, equals 4
						uvxi[i][j][k] = upvxi[i][j][k] + unvxi[i][j][k];
					}

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

void Collision(double* w0, double left, double right, MMDF& m2, MMDF& m3)
{
	// get the equilibrium variables by collision
	w0[0] = left * m2.uplus[0] + right * m3.uminus[0];
	w0[1] = left * m2.uplus[1] + right * m3.uminus[1];
	w0[2] = left * m2.vwhole[1] * m2.uplus[0] + right * m3.vwhole[1] * m3.uminus[0];
	w0[3] = 0.5 * left * (m2.uplus[2] + m2.uplus[0] * m2.vwhole[2] + m2.uplus[0] * m2.xi2) +
		0.5 * right * (m3.uminus[2] + m3.uminus[0] * m3.vwhole[2] + m3.uminus[0] * m3.xi2);
}

void A(double* a, double der[4], double prim[4])
{
	double R4, R3, R2;
	double overden = 1.0 / prim[0];
	//double overlambda = 1.0 / prim[3];
	R4 = der[3] * overden - 0.5 * (prim[1] * prim[1] + prim[2] * prim[2] + 0.5 * (K + 2) / prim[3]) * der[0] * overden;
	R3 = (der[2] - prim[2] * der[0]) * overden;
	R2 = (der[1] - prim[1] * der[0]) * overden;
	a[3] = (4.0 / (K + 2)) * prim[3] * prim[3] * (2 * R4 - 2 * prim[1] * R2 - 2 * prim[2] * R3);
	a[2] = 2 * prim[3] * R3 - prim[2] * a[3];
	a[1] = 2 * prim[3] * R2 - prim[1] * a[3];
	a[0] = der[0] * overden - prim[1] * a[1] - prim[2] * a[2] - 0.5 * a[3] * (prim[1] * prim[1] + prim[2] * prim[2] + 0.5 * (K + 2) / prim[3]);
}

void GL_address(int no_u, int no_v, int no_xi, double* psi, double a[4], MMDF& m)
{
	psi[0] = a[0] * m.upvxi[no_u][no_v][no_xi] + a[1] * m.upvxi[no_u + 1][no_v][no_xi] + a[2] * m.upvxi[no_u][no_v + 1][no_xi] + a[3] * 0.5 * (m.upvxi[no_u + 2][no_v][no_xi] + m.upvxi[no_u][no_v + 2][no_xi] + m.upvxi[no_u][no_v][no_xi + 1]);
	psi[1] = a[0] * m.upvxi[no_u + 1][no_v][no_xi] + a[1] * m.upvxi[no_u + 2][no_v][no_xi] + a[2] * m.upvxi[no_u + 1][no_v + 1][no_xi] + a[3] * 0.5 * (m.upvxi[no_u + 3][no_v][no_xi] + m.upvxi[no_u + 1][no_v + 2][no_xi] + m.upvxi[no_u + 1][no_v][no_xi + 1]);
	psi[2] = a[0] * m.upvxi[no_u][no_v + 1][no_xi] + a[1] * m.upvxi[no_u + 1][no_v + 1][no_xi] + a[2] * m.upvxi[no_u][no_v + 2][no_xi] + a[3] * 0.5 * (m.upvxi[no_u + 2][no_v + 1][no_xi] + m.upvxi[no_u][no_v + 3][no_xi] + m.upvxi[no_u][no_v + 1][no_xi + 1]);
	psi[3] = 0.5 * (a[0] * (m.upvxi[no_u + 2][no_v][no_xi] + m.upvxi[no_u][no_v + 2][no_xi] + m.upvxi[no_u][no_v][no_xi + 1]) +
		a[1] * (m.upvxi[no_u + 3][no_v][no_xi] + m.upvxi[no_u + 1][no_v + 2][no_xi] + m.upvxi[no_u + 1][no_v][no_xi + 1]) +
		a[2] * (m.upvxi[no_u + 2][no_v + 1][no_xi] + m.upvxi[no_u][no_v + 3][no_xi] + m.upvxi[no_u][no_v + 1][no_xi + 1]) +
		a[3] * 0.5 * (m.upvxi[no_u + 4][no_v][no_xi] + m.upvxi[no_u][no_v + 4][no_xi] + m.upvxi[no_u][no_v][no_xi + 2] + 2 * m.upvxi[no_u + 2][no_v + 2][no_xi] + 2 * m.upvxi[no_u + 2][no_v][no_xi + 1] + 2 * m.upvxi[no_u][no_v + 2][no_xi + 1]));
}

void GR_address(int no_u, int no_v, int no_xi, double* psi, double a[4], MMDF& m)
{
	// Similar to the GL_address
	psi[0] = a[0] * m.unvxi[no_u][no_v][no_xi] + a[1] * m.unvxi[no_u + 1][no_v][no_xi] + a[2] * m.unvxi[no_u][no_v + 1][no_xi] + a[3] * 0.5 * (m.unvxi[no_u + 2][no_v][no_xi] + m.unvxi[no_u][no_v + 2][no_xi] + m.unvxi[no_u][no_v][no_xi + 1]);
	psi[1] = a[0] * m.unvxi[no_u + 1][no_v][no_xi] + a[1] * m.unvxi[no_u + 2][no_v][no_xi] + a[2] * m.unvxi[no_u + 1][no_v + 1][no_xi] + a[3] * 0.5 * (m.unvxi[no_u + 3][no_v][no_xi] + m.unvxi[no_u + 1][no_v + 2][no_xi] + m.unvxi[no_u + 1][no_v][no_xi + 1]);
	psi[2] = a[0] * m.unvxi[no_u][no_v + 1][no_xi] + a[1] * m.unvxi[no_u + 1][no_v + 1][no_xi] + a[2] * m.unvxi[no_u][no_v + 2][no_xi] + a[3] * 0.5 * (m.unvxi[no_u + 2][no_v + 1][no_xi] + m.unvxi[no_u][no_v + 3][no_xi] + m.unvxi[no_u][no_v + 1][no_xi + 1]);
	psi[3] = 0.5 * (a[0] * (m.unvxi[no_u + 2][no_v][no_xi] + m.unvxi[no_u][no_v + 2][no_xi] + m.unvxi[no_u][no_v][no_xi + 1]) +
		a[1] * (m.unvxi[no_u + 3][no_v][no_xi] + m.unvxi[no_u + 1][no_v + 2][no_xi] + m.unvxi[no_u + 1][no_v][no_xi + 1]) +
		a[2] * (m.unvxi[no_u + 2][no_v + 1][no_xi] + m.unvxi[no_u][no_v + 3][no_xi] + m.unvxi[no_u][no_v + 1][no_xi + 1]) +
		a[3] * 0.5 * (m.unvxi[no_u + 4][no_v][no_xi] + m.unvxi[no_u][no_v + 4][no_xi] + m.unvxi[no_u][no_v][no_xi + 2] + 2 * m.unvxi[no_u + 2][no_v + 2][no_xi] + 2 * m.unvxi[no_u + 2][no_v][no_xi + 1] + 2 * m.unvxi[no_u][no_v + 2][no_xi + 1]));
}

void G_address(int no_u, int no_v, int no_xi, double* psi, double a[4], MMDF& m)
{

	psi[0] = a[0] * m.uvxi[no_u][no_v][no_xi] + a[1] * m.uvxi[no_u + 1][no_v][no_xi] + a[2] * m.uvxi[no_u][no_v + 1][no_xi] + a[3] * 0.5 * (m.uvxi[no_u + 2][no_v][no_xi] + m.uvxi[no_u][no_v + 2][no_xi] + m.uvxi[no_u][no_v][no_xi + 1]);
	psi[1] = a[0] * m.uvxi[no_u + 1][no_v][no_xi] + a[1] * m.uvxi[no_u + 2][no_v][no_xi] + a[2] * m.uvxi[no_u + 1][no_v + 1][no_xi] + a[3] * 0.5 * (m.uvxi[no_u + 3][no_v][no_xi] + m.uvxi[no_u + 1][no_v + 2][no_xi] + m.uvxi[no_u + 1][no_v][no_xi + 1]);
	psi[2] = a[0] * m.uvxi[no_u][no_v + 1][no_xi] + a[1] * m.uvxi[no_u + 1][no_v + 1][no_xi] + a[2] * m.uvxi[no_u][no_v + 2][no_xi] + a[3] * 0.5 * (m.uvxi[no_u + 2][no_v + 1][no_xi] + m.uvxi[no_u][no_v + 3][no_xi] + m.uvxi[no_u][no_v + 1][no_xi + 1]);
	psi[3] = 0.5 * (a[0] * (m.uvxi[no_u + 2][no_v][no_xi] + m.uvxi[no_u][no_v + 2][no_xi] + m.uvxi[no_u][no_v][no_xi + 1]) +
		a[1] * (m.uvxi[no_u + 3][no_v][no_xi] + m.uvxi[no_u + 1][no_v + 2][no_xi] + m.uvxi[no_u + 1][no_v][no_xi + 1]) +
		a[2] * (m.uvxi[no_u + 2][no_v + 1][no_xi] + m.uvxi[no_u][no_v + 3][no_xi] + m.uvxi[no_u][no_v + 1][no_xi + 1]) +
		a[3] * 0.5 * (m.uvxi[no_u + 4][no_v][no_xi] + m.uvxi[no_u][no_v + 4][no_xi] + m.uvxi[no_u][no_v][no_xi + 2] + 2 * m.uvxi[no_u + 2][no_v + 2][no_xi] + 2 * m.uvxi[no_u + 2][no_v][no_xi + 1] + 2 * m.uvxi[no_u][no_v + 2][no_xi + 1]));
}

double Get_Tau_NS(double density0, double lambda0)
{
	if (tau_type == Euler)
	{
		return 0.0;
	}
	else
	{
		// NS
		if (Mu > 0.0)
		{
			//cout << "here" << endl;
			return 2.0 * Mu * lambda0 / density0;
		}
		else if (Nu > 0.0)
		{
			return 2 * Nu * lambda0;
		}
		else
		{
			return 0.0;
		}
	}
}

double Get_Tau(double density_left, double density_right, double density0, double lambda_left, double lambda_right, double lambda0, double dt)
{
	if (tau_type == Euler)
	{
		if (c1_euler <= 0 && c2_euler <= 0)
		{
			return 0.0;
		}
		else
		{
			double C = c2_euler * abs(density_left / lambda_left - density_right / lambda_right) / abs(density_left / lambda_left + density_right / lambda_right);
			return c1_euler * dt + dt * C;
		}
	}
	else if (tau_type == NS)
	{
		double tau_n = c2_euler * abs(density_left / lambda_left - density_right / lambda_right) / abs(density_left / lambda_left + density_right / lambda_right) * dt;
		if (tau_n != tau_n)
		{
			tau_n = 0.0;
		}
		if ((Mu > 0.0 && Nu > 0.0) || (Mu < 0.0 && Nu < 0.0))
		{
			return 0.0;
		}
		else
		{
			if (Mu > 0.0)
			{
				return tau_n + 2.0 * Mu * lambda0 / density0;
			}
			else if (Nu > 0.0)
			{
				return tau_n + 2 * Nu * lambda0;
			}
			else
			{
				return 0.0;
			}
		}
	}
	else
	{
		return 0.0;
	}
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

bool negative_density_or_pressure(double* primvar)
{
	//detect whether density or pressure is negative
	if (primvar[0] < 0 || primvar[3] < 0 ||
		isnan(primvar[0])
		|| isnan(primvar[3]))
	{
		return true;
	}
	else
	{
		return false;
	}
}

void Convar_to_Primvar(Fluid2d* fluids, Block2d block)
{
	bool all_positive = true;
#pragma omp parallel  for
	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			Convar_to_primvar_2D(fluids[i * block.ny + j].primvar, fluids[i * block.ny + j].convar);
			if (negative_density_or_pressure(fluids[i * block.ny + j].primvar))
			{
				all_positive = false;
			}
		}
	}

	if (all_positive == false)
	{
		cout << "the program blows up at t=" << block.t << "!" << endl;
		exit(0);
	}
}

void Convar_to_primvar_1D(double* primvar, double convar[3])
{
	primvar[0] = convar[0];
	primvar[1] = U(convar[0], convar[1]);
	primvar[2] = Pressure(convar[0], convar[1], convar[2]);

}

void Convar_to_primvar_2D(double *primvar, double convar[4])
{
    primvar[0] = convar[0];
	primvar[1] = U(convar[0], convar[1]);
	primvar[2] = V(convar[0], convar[2]);
	primvar[3] = Pressure(convar[0], convar[1], convar[2], convar[3]);
}

void Convar_to_ULambda_1d(double* primvar, double convar[3])
{
	primvar[0] = convar[0];
	primvar[1] = U(convar[0], convar[1]);
	primvar[2] = Lambda(convar[0], primvar[1], convar[2]);
}

void Convar_to_ULambda_2d(double* primvar, double convar[4])
{
	primvar[0] = convar[0];
	primvar[1] = U(convar[0], convar[1]);
	primvar[2] = V(convar[0], convar[2]);
	//here primvar[3] refers lambda
	primvar[3] = Lambda(convar[0], primvar[1], primvar[2], convar[3]);
}

void Primvar_to_convar_1D(double* convar, double primvar[3])
{
	convar[0] = primvar[0];
	convar[1] = DensityU(primvar[0], primvar[1]);
	convar[2] = DensityE(primvar[0], primvar[1], primvar[2]);
}

void Primvar_to_convar_2D(double *convar, double primvar[4])
{
    convar[0] = primvar[0];
    convar[1] = Q_densityu(primvar[0],primvar[1]);
    convar[2] = Q_densityv(primvar[0],primvar[2]);
    convar[3] = Q_densityE(primvar[0], primvar[1], primvar[2], primvar[3]);
}

void Convar_to_char1D(double* character, double primvar[3], double convar[3])
{
	double c = sqrt(Gamma * primvar[2] / primvar[0]);
	double	alfa = (Gamma - 1.0) / (2.0 * c * c);
	double u = primvar[1];
	double s[3][3];
	s[0][0] = alfa * (0.5 * u * u + u * c / (Gamma - 1.0));
	s[0][1] = alfa * (-u - c / (Gamma - 1.0));
	s[0][2] = alfa;
	s[1][0] = alfa * (-u * u + 2.0 * c * c / (Gamma - 1.0));
	s[1][1] = alfa * 2.0 * u;
	s[1][2] = -2.0 * alfa;
	s[2][0] = alfa * (0.5 * u * u - u * c / (Gamma - 1.0));
	s[2][1] = alfa * (-u + c / (Gamma - 1.0));
	s[2][2] = alfa;

	for (int i = 0; i < 3; i++)
	{
		character[i] = 0;
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			character[i] = character[i] + s[i][j] * convar[j];
		}
	}
}

void Char_to_convar1D(double* convar, double primvar[3], double charvar[3])
{
	double c = sqrt(Gamma * primvar[2] / primvar[0]);
	double u = primvar[1];
	double	h = 0.5 * u * u + c * c / (Gamma - 1.0);
	double s[3][3];
	s[0][0] = 1.0;
	s[0][1] = 1.0;
	s[0][2] = 1.0;
	s[1][0] = u - c;
	s[1][1] = u;
	s[1][2] = u + c;
	s[2][0] = h - u * c;
	s[2][1] = u * u / 2.0;
	s[2][2] = h + u * c;

	for (int i = 0; i < 3; i++)
	{
		convar[i] = 0;
		for (int j = 0; j < 3; j++)
		{
			convar[i] = convar[i] + s[i][j] * charvar[j];
		}
	}
}

void Convar_to_char(double *character, double *primvar, double convar[4])
{
	double c = sqrt(Gamma *primvar[3] / primvar[0]);
	double overc = 1.0 / c;
	double	alfa = (Gamma - 1.0)*0.5*overc*overc;
	double u = primvar[1];
	double v = primvar[2];
	double s[4][4];
	s[0][0] = 1.0 - alfa*(u*u + v*v);
	s[0][1] = 2.0*alfa*u;
	s[0][2] = 2.0*alfa*v;
	s[0][3] = -2.0*alfa;
	s[1][0] = -v;
	s[1][1] = 0.;
	s[1][2] = 1.;
	s[1][3] = 0.;
	s[2][0] = 0.5*(-u *overc + alfa*(u*u + v*v));
	s[2][1] = 0.5 *overc - alfa*u;
	s[2][2] = -alfa*v;
	s[2][3] = alfa;
	s[3][0] = 0.5*(u *overc + alfa*(u*u + v*v));
	s[3][1] = -0.5 *overc - alfa*u;
	s[3][2] = -alfa*v;
	s[3][3] = alfa;
	for (int i = 0; i < 4; i++)
	{
		character[i] = 0;
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			character[i] = character[i] + s[i][j] * convar[j];
		}
	}
}

void Char_to_convar(double *convar, double *primvar, double character[4])
{
	double c = sqrt(Gamma *primvar[3] / primvar[0]);
	double	alfa = (Gamma - 1.0) / (2.0*c*c);
	double u = primvar[1];
	double v = primvar[2];
	double s[4][4];
	s[0][0] = 1.;
	s[0][1] = 0.0;
	s[0][2] = 1.0;
	s[0][3] = 1.0;
	s[1][0] = u;
	s[1][1] = 0.;
	s[1][2] = u + c;
	s[1][3] = u - c;
	s[2][0] = v;
	s[2][1] = 1.;
	s[2][2] = v;
	s[2][3] = v;
	s[3][0] = 0.5*(u*u + v*v);
	s[3][1] = v;
	s[3][2] = 0.5*(u*u + v*v) + c*u + c*c / (Gamma - 1.0);
	s[3][3] = 0.5*(u*u + v*v) - c*u + c*c / (Gamma - 1.0);

	for (int i = 0; i < 4; i++)
	{
		convar[i] = 0;
		for (int j = 0; j < 4; j++)
		{
			convar[i] = convar[i] + s[i][j] * character[j];
		}
	}
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

double Q_densityu(double density, double u)
{
	double q_densityu = density*u;
	return q_densityu;
}

double Q_densityv(double density, double v)
{
	double q_densityv = density*v;
	return q_densityv;
}

double Q_densityE(double density, double u, double v, double pressure)
{
	double q_densityE = density*(pressure / density / (Gamma - 1) + 0.5*(u*u + v*v));
	return q_densityE;
}

double U(double density, double q_densityu)
{
	return q_densityu / density;
}

double V(double density, double q_densityv)
{
	return q_densityv / density;
}

double Pressure(double density, double q_densityu, double q_densityv, double q_densityE)
{
	return (Gamma - 1.0)*(q_densityE - 0.5*q_densityu*q_densityu / density - 0.5*q_densityv*q_densityv / density);
}

double Lambda(double density, double u, double densityE)
{
	return (K + 1.0) * 0.25 * (density / (densityE - 0.5 * density * (u * u)));
}

double Lambda(double density, double u, double v, double densityE)
{
	return (K + 2.0) * 0.25 * (density / (densityE - 0.5 * density * (u * u + v * v)));
}

void Copy_Array(double* target, double* origin, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		target[i] = origin[i];
	}
}

void YchangetoX(double* fluidtmp, double* fluid)
{
	fluidtmp[0] = fluid[0];
	fluidtmp[1] = fluid[2];
	fluidtmp[2] = -fluid[1];
	fluidtmp[3] = fluid[3];
}

void Local_to_Global(double *change,double *normal)
{
    double temp[2];
    temp[0] = change[1];
    temp[1] = change[2];
    change[1] = temp[0] * normal[0] - temp[1] * normal[1];
    change[2] = temp[0] * normal[1] + temp[1] * normal[0];
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

Fluid2d* Setfluid(Block2d& block)
{
	Fluid2d* var = new Fluid2d[block.nx * block.ny]; // dynamic variable (since block.nx is not determined)
	if (var == 0)
	{
		cout << "fluid variable allocate fail...";
		return NULL;
	}
	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			var[i * block.ny + j].xindex = i;
			var[i * block.ny + j].yindex = j;
		}
	}
	cout << "fluid variable allocate done..." << endl;
	return var;
}

Interface2d* Setinterface_array(Block2d block)
{
	Interface2d* var = new Interface2d[(block.nx + 1) * (block.ny + 1)];  // dynamic variable (since block.nx is not determined)
	if (var == 0)
	{
		cout << "fluid variable allocate fail...";
		return NULL;
	}
	for (int i = 0; i < block.nx + 1; i++)
	{
		for (int j = 0; j < block.ny + 1; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				var[i * (block.ny + 1) + j].line.left.der1x[k] = 0.0;
				var[i * (block.ny + 1) + j].line.left.der1y[k] = 0.0;

				var[i * (block.ny + 1) + j].line.right.der1x[k] = 0.0;
				var[i * (block.ny + 1) + j].line.right.der1y[k] = 0.0;

				var[i * (block.ny + 1) + j].line.center.der1x[k] = 0.0;
				var[i * (block.ny + 1) + j].line.center.der1y[k] = 0.0;
			}
		}
	}
	cout << "interface variable allocate done..." << endl;
	return var;
}

Flux2d_gauss** Setflux_gauss_array(Block2d block)
{
	Flux2d_gauss** var = new Flux2d_gauss * [(block.nx + 1) * (block.ny + 1)];  // dynamic variable (since block.nx is not determined)

	for (int i = 0; i < block.nx + 1; i++)
	{
		for (int j = 0; j < block.ny + 1; j++)
		{
			// for m th step time marching schemes, m subflux needed
			var[i * (block.ny + 1) + j] = new Flux2d_gauss[block.stages];
		}
	}

	for (int i = 0; i < block.nx + 1; i++)
	{
		for (int j = 0; j < block.ny + 1; j++)
		{
			for (int k = 0; k < block.stages; k++)
			{
				if (gausspoint == 0)
				{
					var[i * (block.ny + 1) + j][k].gauss = new Flux2d[1];
					for (int m = 0; m < 4; m++)
					{
						var[i * (block.ny + 1) + j][k].gauss[0].f[m] = 0.0;
						var[i * (block.ny + 1) + j][k].gauss[0].derf[m] = 0.0;
					}
				}
				else
				{
					var[i * (block.ny + 1) + j][k].gauss = new Flux2d[gausspoint];
					for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
					{
						for (int m = 0; m < 4; m++)
						{
							var[i * (block.ny + 1) + j][k].gauss[num_gauss].f[m] = 0.0;
							var[i * (block.ny + 1) + j][k].gauss[num_gauss].derf[m] = 0.0;
						}
					}
				}
			}
		}
	}
	if (var == 0)
	{
		cout << "fluid variable allocate fail...";
		return NULL;
	}
	cout << "flux with gausspoint variable allocate done..." << endl;
	return var;
}

void SetUniformMesh(Block2d& block, Fluid2d* fluids, Interface2d* xinterfaces, Interface2d* yinterfaces, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes)
{
	//nodex, nodey are the real node
	//interface number = cell number + 1
	block.dx = (block.right - block.left) / block.nodex;
	block.dy = (block.up - block.down) / block.nodey;
	block.overdx = 1 / block.dx;
	block.overdy = 1 / block.dy;

	block.xcell_begin = block.ghost;
	block.xcell_end = block.ghost + block.nodex - 1;
	block.ycell_begin = block.ghost;
	block.ycell_end = block.ghost + block.nodey - 1;

	block.xinterface_begin_n = block.ghost;
	block.xinterface_end_n = block.ghost + block.nodex;
	block.xinterface_begin_t = block.ghost;
	block.xinterface_end_t = block.ghost + block.nodex - 1;

	block.yinterface_begin_n = block.ghost;
	block.yinterface_end_n = block.ghost + block.nodey;
	block.yinterface_begin_t = block.ghost;
	block.yinterface_end_t = block.ghost + block.nodey - 1;

	//cell avg information
	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			//two dimension geometry to one dimension store matrix, y direciton first and x direction second
			fluids[i * block.ny + j].dx = block.dx; //cell size
			fluids[i * block.ny + j].dy = block.dy; //cell size
			fluids[i * block.ny + j].coordx = block.left + (i + 0.5 - block.ghost) * block.dx; //cell center location
			fluids[i * block.ny + j].coordy = block.down + (j + 0.5 - block.ghost) * block.dy; //cell center location
			fluids[i * block.ny + j].area = block.dx * block.dy;
			fluids[i * block.ny + j].node[0] = block.left + (i - block.ghost) * block.dx;
			fluids[i * block.ny + j].node[1] = block.down + (j - block.ghost) * block.dy;
			fluids[i * block.ny + j].node[2] = block.left + (i + 1 - block.ghost) * block.dx;
			fluids[i * block.ny + j].node[3] = block.down + (j - block.ghost) * block.dy;
			fluids[i * block.ny + j].node[4] = block.left + (i + 1 - block.ghost) * block.dx;
			fluids[i * block.ny + j].node[5] = block.down + (j + 1 - block.ghost) * block.dy;
			fluids[i * block.ny + j].node[6] = block.left + (i - block.ghost) * block.dx;
			fluids[i * block.ny + j].node[7] = block.down + (j + 1 - block.ghost) * block.dy;
		}
	}

	// interface information
	for (int i = 0; i <= block.nx; i++)
	{
		for (int j = 0; j <= block.ny; j++)
		{
			xinterfaces[i * (block.ny + 1) + j].x = block.left + (i - block.ghost) * block.dx;
			xinterfaces[i * (block.ny + 1) + j].y = block.down + (j - block.ghost + 0.5) * block.dy;
			xinterfaces[i * (block.ny + 1) + j].length = block.dy;
			xinterfaces[i * (block.ny + 1) + j].normal[0] = 1.0;
			xinterfaces[i * (block.ny + 1) + j].normal[1] = 0.0;

			Copy_geo_from_interface_to_line(xinterfaces[i * (block.ny + 1) + j]);
			xinterfaces[i * (block.ny + 1) + j].gauss = new Recon2d[gausspoint];

			Copy_geo_from_interface_to_flux
			(xinterfaces[i * (block.ny + 1) + j], xfluxes[i * (block.ny + 1) + j], block.stages);

			yinterfaces[i * (block.ny + 1) + j].y = block.down + (j - block.ghost) * block.dy;
			yinterfaces[i * (block.ny + 1) + j].x = block.left + (i - block.ghost + 0.5) * block.dx;
			yinterfaces[i * (block.ny + 1) + j].length = block.dx;
			yinterfaces[i * (block.ny + 1) + j].normal[0] = 0.0;
			yinterfaces[i * (block.ny + 1) + j].normal[1] = 1.0;

			Copy_geo_from_interface_to_line(yinterfaces[i * (block.ny + 1) + j]);

			yinterfaces[i * (block.ny + 1) + j].gauss = new Recon2d[gausspoint];

			Copy_geo_from_interface_to_flux
			(yinterfaces[i * (block.ny + 1) + j], yfluxes[i * (block.ny + 1) + j], block.stages);

			Set_Gauss_Coordinate(xinterfaces[i * (block.ny + 1) + j], yinterfaces[i * (block.ny + 1) + j]);
		}
	}
	cout << "set uniform information done..." << endl;
}

void Copy_geo_from_interface_to_line(Interface2d& interface)
{
	interface.line.x = interface.x;
	interface.line.y = interface.y;
	interface.line.normal[0] = interface.normal[0];
	interface.line.normal[1] = interface.normal[1];
}

void Copy_geo_from_interface_to_flux(Interface2d& interface, Flux2d_gauss* flux, int stages)
{
	for (int istage = 0; istage < stages; istage++)
	{
		if (gausspoint == 0)
		{
			int igauss = 0;
			flux[istage].gauss[igauss].normal[0] = interface.normal[0];
			flux[istage].gauss[igauss].normal[1] = interface.normal[1];
			flux[istage].gauss[igauss].length = interface.length;
		}
		else
		{
			for (int igauss = 0; igauss < gausspoint; igauss++)
			{
				flux[istage].gauss[igauss].normal[0] = interface.normal[0];
				flux[istage].gauss[igauss].normal[1] = interface.normal[1];
				flux[istage].gauss[igauss].length = interface.length;
			}
		}
	}
}

void Set_Gauss_Coordinate(Interface2d& xinterface, Interface2d& yinterface)
{
	for (int num_guass = 0; num_guass < gausspoint; num_guass++)
	{
		// first is gauss parameter
		xinterface.gauss[num_guass].x = xinterface.x;
		xinterface.gauss[num_guass].y = xinterface.y + gauss_loc[num_guass] * 0.5 * xinterface.length;
		xinterface.gauss[num_guass].normal[0] = xinterface.normal[0];
		xinterface.gauss[num_guass].normal[1] = xinterface.normal[1];
		// each gauss point contain left center right point
		Set_Gauss_Intergation_Location_x(xinterface.gauss[num_guass].left, num_guass, xinterface.length);
		Set_Gauss_Intergation_Location_x(xinterface.gauss[num_guass].right, num_guass, xinterface.length);
		Set_Gauss_Intergation_Location_x(xinterface.gauss[num_guass].center, num_guass, xinterface.length);

		yinterface.gauss[num_guass].x = yinterface.x + gauss_loc[num_guass] * 0.5 * yinterface.length;
		yinterface.gauss[num_guass].y = yinterface.y;
		yinterface.gauss[num_guass].normal[0] = yinterface.normal[0];
		yinterface.gauss[num_guass].normal[1] = yinterface.normal[1];
		// each gauss point contain left center right point
		Set_Gauss_Intergation_Location_y(yinterface.gauss[num_guass].left, num_guass, yinterface.length);
		Set_Gauss_Intergation_Location_y(yinterface.gauss[num_guass].right, num_guass, yinterface.length);
		Set_Gauss_Intergation_Location_y(yinterface.gauss[num_guass].center, num_guass, yinterface.length);
	}
}

void Set_Gauss_Intergation_Location_x(Point2d& xgauss, int index, double h)
{
	xgauss.x = gauss_loc[index] * h / 2.0;
}

void Set_Gauss_Intergation_Location_y(Point2d& ygauss, int index, double h)
{
	ygauss.x = gauss_loc[index] * h / 2.0;
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

void CopyFluid_new_to_old(Fluid2d* fluids, Block2d block)
{
#pragma omp parallel  for
	for (int i = block.ghost; i < block.ghost + block.nodex; i++)
	{
		for (int j = block.ghost; j < block.ghost + block.nodey; j++)
		{
			for (int var = 0; var < 4; var++)
			{
				fluids[i * block.ny + j].convar_old[var] = fluids[i * block.ny + j].convar[var];
			}
		}
	}
}

