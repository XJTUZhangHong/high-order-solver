#include "time_advance.h"

// one-dimensional problem
TimeMarchingCoefficient timecoe_list = RK2; //initialization

void S1O1(Block1d& block)
{
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[0][0][1] = 0.0;
}

void S1O2(Block1d& block)
{
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[0][0][1] = 0.5;
}

void S2O4(Block1d& block)
{
	// two stages, so the extra stages coefficients are zero.
	block.stages = 2;
	block.timecoefficient[0][0][0] = 0.5;
	block.timecoefficient[0][0][1] = 1.0 / 8.0;
	block.timecoefficient[1][0][0] = 1.0;
	block.timecoefficient[1][1][0] = 0.0;
	block.timecoefficient[1][0][1] = 1.0 / 6.0;
	block.timecoefficient[1][1][1] = 1.0 / 3.0;
}

void RK2(Block1d& block)
{
	block.stages = 2;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[1][0][0] = 0.5;
	block.timecoefficient[1][1][0] = 0.5;

}

void RK3(Block1d& block)
{
	block.stages = 3;
}

void RK4(Block1d& block)
{
	block.stages = 4;
	block.timecoefficient[0][0][0] = 0.5;
	block.timecoefficient[1][1][0] = 0.5;
	block.timecoefficient[2][2][0] = 1.0;
	block.timecoefficient[3][0][0] = 1.0 / 6.0;
	block.timecoefficient[3][1][0] = 1.0 / 3.0;
	block.timecoefficient[3][2][0] = 1.0 / 3.0;
	block.timecoefficient[3][3][0] = 1.0 / 6.0;
}

void Initial_stages(Block1d& block)
{
	// initialize five stages (5 steps) and set all them as ZERO
	for (int i = 0; i < 5; i++) //refers the n stage
	{
		for (int j = 0; j < 5; j++) //refers the nth coefficient at n stage
		{
			for (int k = 0; k < 3; k++) //refers f, derf
			{
				block.timecoefficient[i][j][k] = 0.0;
			}
		}
	}
	// by timecoe_list set the correct stages (others still be ZERO)
	timecoe_list(block);
}

double Get_CFL(Block1d& block, Fluid1d* fluids, double tstop)
{
	double dt = block.dx;
	for (int i = 0; i < block.nodex; i++)
	{
		dt = Dtx(dt, block.dx, block.CFL, fluids[i + block.ghost].convar);
	}
	if (block.t + dt > tstop)
	{
		dt = tstop - block.t + 1e-15;
	}
	//print time step information
	cout << "step= " << block.step
		<< "time size is " << dt
		<< " time= " << block.t << endl;
	return dt;
}

double Dtx(double dtx, double dx, double CFL, double convar[3])
{
	double tmp;
	double prim[3];
	Convar_to_primvar_1D(prim, convar);
	tmp = abs(prim[1]) + sqrt(Gamma * prim[2] / prim[0]);
	if (tmp > CFL * dx / dtx) // if abs(u)+c (tmp) > abs(u)+c (old)
	{
		dtx = CFL * dx / tmp; // abs(u)+c determine one smaller time step
	}
	//consider viscous problem (tau_type == NS)
	if (dtx > 0.25 * CFL * dx * dx / Mu && tau_type == NS && Mu > 0)
	{
		// Note: if dtx (above) > CFL * dx * dx / 4
		// replace dxt by the smaller time step determined by viscous term
		dtx = 0.25 * CFL * dx * dx / Mu;
	}
	return dtx;
}

void Update(Fluid1d* fluids, Flux1d** fluxes, Block1d block, int stage)
{
	if (stage > block.stages)
	{
		cout << "wrong middle stage,pls check the time marching setting" << endl;
		exit(0);
	}

	double dt = block.dt;
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost + 1; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			double Flux = 0.0;
			for (int k = 0; k < stage + 1; ++k)
			{
				Flux = Flux
					+ block.timecoefficient[stage][k][0] * fluxes[i][k].f[j]
					+ block.timecoefficient[stage][k][1] * fluxes[i][k].derf[j];
			}
			fluxes[i][stage].F[j] = Flux;
		}
	}
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			fluids[i].convar[j] = fluids[i].convar_old[j] + 1.0 / fluids[i].dx * (fluxes[i][stage].F[j] - fluxes[i + 1][stage].F[j]);
		}
	}

}