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
	cout << "step = " << block.step
		<< " time size is " << dt
		<< " time = " << block.t << endl;
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
	if (timecoe_list == RK3)
	{
		for (int i = block.ghost; i < block.nodex + block.ghost + 1; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				double Flux = 0;
				if (stage == 0)
				{
					fluxes[i][stage].F[j] = fluxes[i][0].f[j];
				}
				if (stage == 1)
				{
					fluxes[i][stage].F[j] = 0.25 * fluxes[i][1].f[j];
				}
				if (stage == 2)
				{
					fluxes[i][stage].F[j] = 2.0 / 3.0 * fluxes[i][2].f[j];
				}
			}
		}
		for (int i = block.ghost; i < block.nodex + block.ghost; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				if (stage == 0)
				{
					fluids[i].convar[j] = fluids[i].convar_old[j] 
						+ 1.0 / fluids[i].dx * (fluxes[i][stage].F[j] - fluxes[i + 1][stage].F[j]);
					fluids[i].convar_stage1[j] = fluids[i].convar[j];
				}
				if (stage == 1)
				{
					fluids[i].convar[j] = 0.75 * fluids[i].convar_old[j] + 0.25 * fluids[i].convar_stage1[j]
						+ 1.0 / fluids[i].dx * (fluxes[i][stage].F[j] - fluxes[i + 1][stage].F[j]);
					fluids[i].convar_stage2[j] = fluids[i].convar[j];
				}
				if (stage == 2)
				{
					fluids[i].convar[j] = 1.0 / 3.0 * fluids[i].convar_old[j] + 2.0 / 3.0 * fluids[i].convar_stage2[j]
						+ 1.0 / fluids[i].dx * (fluxes[i][stage].F[j] - fluxes[i + 1][stage].F[j]);
				}
			}
		}
	}
	else
	{
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

}

// two-dimensional problem
TimeMarchingCoefficient_2d timecoe_list_2d = S1O1_2D; //initilization

void S1O1_2D(Block2d& block)
{
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[0][0][1] = 0.0;
}

void S1O2_2D(Block2d& block)
{
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[0][0][1] = 0.5;
}

void RK2_2D(Block2d& block)
{
	block.stages = 2;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[1][0][0] = 0.5;
	block.timecoefficient[1][1][0] = 0.5;
}

void S2O4_2D(Block2d& block)
{
	block.stages = 2;
	block.timecoefficient[0][0][0] = 0.5;
	block.timecoefficient[0][0][1] = 1.0 / 8.0;
	block.timecoefficient[1][0][0] = 1.0;
	block.timecoefficient[1][1][0] = 0.0;
	block.timecoefficient[1][0][1] = 1.0 / 6.0;
	block.timecoefficient[1][1][1] = 1.0 / 3.0;
}

void RK3_2D(Block2d& block)
{
	block.stages = 3;
}

void RK4_2D(Block2d& block)
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

void Initial_stages(Block2d& block)
{
	for (int i = 0; i < 5; i++) //refers the n stage
	{
		for (int j = 0; j < 5; j++) //refers the nth coefficient at n stage
		{
			for (int k = 0; k < 3; k++) //refers f, derf, der2f
			{
				block.timecoefficient[i][j][k] = 0.0;
			}
		}
	}
	timecoe_list_2d(block);
}

double Get_CFL(Block2d& block, Fluid2d* fluids, double tstop)
{
	double dt = fluids[block.ghost * block.ny + block.ghost].dx;

	for (int i = block.ghost; i < block.nx - block.ghost; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost; j++)
		{
			dt = Dtx(dt, fluids[i * block.ny + j].dx, block.CFL, fluids[i * block.ny + j].primvar[0], fluids[i * block.ny + j].primvar[1],
				fluids[i * block.ny + j].primvar[2], fluids[i * block.ny + j].primvar[3]);
			dt = Dtx(dt, fluids[i * block.ny + j].dy, block.CFL, fluids[i * block.ny + j].primvar[0], fluids[i * block.ny + j].primvar[1],
				fluids[i * block.ny + j].primvar[2], fluids[i * block.ny + j].primvar[3]);
		}
	}
	if (block.step % 10 == 0)
	{
		cout << "step = " << block.step
		<< " time size is " << dt
		<< " time = " << block.t << endl;
	}
	if (block.t + dt > tstop)
	{
		dt = tstop - block.t + 1e-16;
		cout << "last step here" << endl;
	}
	return dt;
}

double Dtx(double dtx, double dx, double CFL, double density, double u,double v, double pressure)
{
	double tmp;
	tmp = sqrt(u*u+v*v) + sqrt(Gamma *pressure / density);
	if (tmp>CFL*dx / dtx)
	{
		dtx = CFL*dx / tmp;
	}
	if (dtx > 0.25*CFL*dx*dx / Mu&&tau_type == NS&&Mu>0)
	{
		dtx = 0.25*CFL*dx*dx / Mu;
	}
	return dtx;
}

void Update(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, int stage)
{
	// Note : write the function to update the conservative variables of each cell
	if (stage > block.stages)
	{
		cout << "wrong middle stage,pls check the time marching setting" << endl;
		exit(0);
	}

	double dt = block.dt;
	if (timecoe_list_2d == RK3_2D)
	{
		#pragma omp parallel  for
		for (int i = block.ghost; i < block.nodex + block.ghost + 1; i++)
		{
			for (int j = block.ghost; j < block.nodey + block.ghost; j++)
			{
				for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
				{
					for (int var = 0; var < 4; var++)
					{
						if (stage == 0)
						{
							xfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].x[var] = 
								gauss_weight[num_gauss] * xfluxes[i * (block.ny + 1) + j][0].gauss[num_gauss].f[var];
						}
						if (stage == 1)
						{
							xfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].x[var] = 0.25 *
								gauss_weight[num_gauss] * xfluxes[i * (block.ny + 1) + j][1].gauss[num_gauss].f[var];
						}
						if (stage == 2)
						{
							xfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].x[var] = 2.0 / 3.0 *
								gauss_weight[num_gauss] * xfluxes[i * (block.ny + 1) + j][2].gauss[num_gauss].f[var];
						}
					}
				}
			}
		}
	#pragma omp parallel  for
		for (int i = block.ghost; i < block.nodex + block.ghost; i++)
		{
			for (int j = block.ghost; j < block.nodey + block.ghost + 1; j++)
			{
				for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
				{
					for (int var = 0; var < 4; var++)
					{
						if (stage == 0)
						{
							yfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].x[var] =
								gauss_weight[num_gauss] * yfluxes[i * (block.ny + 1) + j][0].gauss[num_gauss].f[var];
						}
						if (stage == 1)
						{
							yfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].x[var] = 0.25 *
								gauss_weight[num_gauss] * yfluxes[i * (block.ny + 1) + j][1].gauss[num_gauss].f[var];
						}
						if (stage == 2)
						{
							yfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].x[var] = 2.0 / 3.0 *
								gauss_weight[num_gauss] * yfluxes[i * (block.ny + 1) + j][2].gauss[num_gauss].f[var];
						}
					}
				}
			}
		}

		Update_with_gauss_RK3(fluids, xfluxes, yfluxes, block, stage);
	}
	else
	{
		#pragma omp parallel  for
		for (int i = block.ghost; i < block.nodex + block.ghost + 1; i++)
		{
			for (int j = block.ghost; j < block.nodey + block.ghost; j++)
			{
				for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
				{
					for (int var = 0; var < 4; var++)
					{
						double Flux = 0.0;
						for (int k = 0; k < stage + 1; k++)
						{
		
							Flux = Flux
								+ gauss_weight[num_gauss] *
								(block.timecoefficient[stage][k][0] * xfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].f[var]
									+ block.timecoefficient[stage][k][1] * xfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].derf[var]);

						}
						xfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].x[var] = Flux;
					}
				}
			}
		}

	#pragma omp parallel  for
		for (int i = block.ghost; i < block.nodex + block.ghost; i++)
		{
			for (int j = block.ghost; j < block.nodey + block.ghost + 1; j++)
			{
				for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
				{
					for (int var = 0; var < 4; var++)
					{
						double Flux = 0.0;
						for (int k = 0; k < stage + 1; k++)
						{
							Flux = Flux
								+ gauss_weight[num_gauss] *
								(block.timecoefficient[stage][k][0] * yfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].f[var]
									+ block.timecoefficient[stage][k][1] * yfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].derf[var]);
						}
						yfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].x[var] = Flux;
					}
				}
			}
		}
		Update_with_gauss(fluids, xfluxes, yfluxes, block, stage);
		}
	}

void Update_with_gauss(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, int stage)
{
	//Note : calculate the final flux of the cell, in fluids, by the obtained final flux of interface, in xfluxes and yfluxes
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost + 1; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost + 1; j++)
		{
			int face = i * (block.ny + 1) + j;
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				Local_to_Global(yfluxes[face][stage].gauss[num_gauss].x, yfluxes[face][stage].gauss[num_gauss].normal);
				Local_to_Global(xfluxes[face][stage].gauss[num_gauss].x, xfluxes[face][stage].gauss[num_gauss].normal);
			}
		}
	}

#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			int cell = i * (block.ny) + j;
			int face = i * (block.ny + 1) + j;
			
			double total_flux;
			for (int var = 0; var < 4; var++)
			{
				fluids[cell].convar[var] = fluids[cell].convar_old[var]; //get the Wn from convar_old
				total_flux = 0.0;
				for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
				{
					total_flux += yfluxes[face][stage].gauss[num_gauss].length * yfluxes[face][stage].gauss[num_gauss].x[var];
					total_flux += -yfluxes[face + 1][stage].gauss[num_gauss].length * yfluxes[face + 1][stage].gauss[num_gauss].x[var];
					total_flux += xfluxes[face][stage].gauss[num_gauss].length * xfluxes[face][stage].gauss[num_gauss].x[var];
					total_flux += -xfluxes[face + block.ny + 1][stage].gauss[num_gauss].length * xfluxes[face + block.ny + 1][stage].gauss[num_gauss].x[var];
				}
				fluids[cell].convar[var] += total_flux / fluids[cell].area;
				// calculate the final flux of the cell, in fluids, by the obtained final flux of interface, in xfluxes and yfluxes
				// in 2d, the total flux be updated by the line-averaged flux of four interface
			}
		}
	}
}

void Update_with_gauss_RK3(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, int stage)
{
	//Note : calculate the final flux of the cell, in fluids, by the obtained final flux of interface, in xfluxes and yfluxes
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost + 1; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost + 1; j++)
		{
			int face = i * (block.ny + 1) + j;
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				Local_to_Global(yfluxes[face][stage].gauss[num_gauss].x, yfluxes[face][stage].gauss[num_gauss].normal);
				Local_to_Global(xfluxes[face][stage].gauss[num_gauss].x, xfluxes[face][stage].gauss[num_gauss].normal);
			}
		}
	}

#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			int cell = i * (block.ny) + j;
			int face = i * (block.ny + 1) + j;

			for (int var = 0; var < 4; var++)
			{
				double total_flux = 0.0;
				for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
				{
					total_flux += yfluxes[face][stage].gauss[num_gauss].length * yfluxes[face][stage].gauss[num_gauss].x[var];
					total_flux += -yfluxes[face + 1][stage].gauss[num_gauss].length * yfluxes[face + 1][stage].gauss[num_gauss].x[var];
					total_flux += xfluxes[face][stage].gauss[num_gauss].length * xfluxes[face][stage].gauss[num_gauss].x[var];
					total_flux += -xfluxes[face + block.ny + 1][stage].gauss[num_gauss].length * xfluxes[face + block.ny + 1][stage].gauss[num_gauss].x[var];
				}
				if (stage == 0)
				{
					fluids[cell].convar[var] = fluids[cell].convar_old[var] + total_flux / fluids[cell].area;
					fluids[cell].convar_stage1[var] = fluids[cell].convar[var];
				}
				if (stage == 1)
				{
					fluids[cell].convar[var] = 0.75 * fluids[cell].convar_old[var] + 0.25 * fluids[cell].convar_stage1[var] + 
						total_flux / fluids[cell].area;
					fluids[cell].convar_stage2[var] = fluids[cell].convar[var];
				}
				if (stage == 2)
				{
					fluids[cell].convar[var] = 1.0 / 3.0 * fluids[cell].convar_old[var] + 2.0 / 3.0 * fluids[cell].convar_stage2[var] +
						total_flux / fluids[cell].area;
				}
				// calculate the final flux of the cell, in fluids, by the obtained final flux of interface, in xfluxes and yfluxes
				// in 2d, the total flux be updated by the line-averaged flux of four interface
			}
		}
	}
}

void Update_with_source(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, int stage)
{
	double dt = block.dt;
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost + 1; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost + 1; j++)
		{
			int face = i * (block.ny + 1) + j;
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				Local_to_Global(yfluxes[face][stage].gauss[num_gauss].f, yfluxes[face][stage].gauss[num_gauss].normal);
				Local_to_Global(xfluxes[face][stage].gauss[num_gauss].f, xfluxes[face][stage].gauss[num_gauss].normal);
			

				Local_to_Global(yfluxes[face][stage].gauss[num_gauss].derf, yfluxes[face][stage].gauss[num_gauss].normal);
				Local_to_Global(xfluxes[face][stage].gauss[num_gauss].derf, xfluxes[face][stage].gauss[num_gauss].normal);

			}
		}
	}
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost + 1; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				for (int var = 0; var < 4; var++)
				{
					double Flux1 = 0.0;
					for (int k = 0; k < stage + 1; k++)
					{
						Flux1 = Flux1
							+ gauss_weight[num_gauss] *
							(block.timecoefficient[stage][k][0] * xfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].f[var]
								+ block.timecoefficient[stage][k][1] * xfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].derf[var]);
					}
					xfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].x[var] = Flux1;
				}
			}
		}
	}
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost + 1; j++)
		{
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				for (int var = 0; var < 4; var++)
				{
					double Flux1 = 0.0;
					for (int k = 0; k < stage + 1; k++)
					{
						Flux1 = Flux1
							+ gauss_weight[num_gauss] *
							(block.timecoefficient[stage][k][0] * yfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].f[var]
								+ block.timecoefficient[stage][k][1] * yfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].derf[var]);
					}
					yfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].x[var] = Flux1;
				}
			}
		}
	}
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			int cell = i * (block.ny) + j;
			int face = i * (block.ny + 1) + j;
			double Lw1[4]{ 0, 0, 0, 0 };
			double derLw1[4]{ 0, 0, 0, 0 };

			for (int var = 0; var < 4; var++)
			{
				for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
				{
					Lw1[var] += gauss_weight[num_gauss] * ((xfluxes[face][stage].gauss[num_gauss].f[var] - xfluxes[face + block.ny + 1][stage].gauss[num_gauss].f[var]) / block.dx / block.dt
						+ (yfluxes[face][stage].gauss[num_gauss].f[var] - yfluxes[face + 1][stage].gauss[num_gauss].f[var]) / block.dy / block.dt);
				
					derLw1[var] += gauss_weight[num_gauss] * ((xfluxes[face][stage].gauss[num_gauss].derf[var] - xfluxes[face + block.ny + 1][stage].gauss[num_gauss].derf[var]) / block.dx / block.dt / block.dt
						+ (yfluxes[face][stage].gauss[num_gauss].derf[var] - yfluxes[face + 1][stage].gauss[num_gauss].derf[var]) / block.dy / block.dt / block.dt);
				}
			}
			fluids[cell].Sw1[stage][0] = 0.0; fluids[cell].Sw1[stage][1] = 0.0;
			fluids[cell].Sw1[stage][2] = fluids[cell].convar[0];
			fluids[cell].Sw1[stage][3] = fluids[cell].convar[2];
			for (int var = 0; var < 4; var++)
			{
				fluids[cell].Lw1[stage][var] = Lw1[var] + fluids[cell].Sw1[stage][var];
			}
			fluids[cell].derSw1[stage][0] = 0.0; fluids[cell].derSw1[stage][1] = 0.0;
			fluids[cell].derSw1[stage][2] = fluids[cell].Lw1[stage][0];
			fluids[cell].derSw1[stage][3] = fluids[cell].Lw1[stage][2];
			for (int var = 0; var < 4; var++)
			{
				fluids[cell].derLw1[stage][var] = derLw1[var] + fluids[cell].derSw1[stage][var];
			}

			// RK2
			if (timecoe_list_2d == RK2_2D)
			{
				if (stage == 0)
				{
					for (int k = 0; k < 4; k++)
					{
						fluids[cell].convar[k] = fluids[cell].convar_old[k] + block.dt * fluids[cell].Lw1[stage][k];
					}
				}
				if (stage == 1)
				{
					for (int k = 0; k < 4; k++)
					{
						fluids[cell].convar[k] = fluids[cell].convar_old[k] + 0.5 * block.dt * fluids[cell].Lw1[0][k]
							+ 0.5 * block.dt * fluids[cell].Lw1[1][k];
					}
				}
			}
			// RK3
			if (timecoe_list_2d == RK3_2D)
			{
				if (stage == 0)
				{
					for (int k = 0; k < 4; k++)
					{
						fluids[cell].convar[k] = fluids[cell].convar_old[k] + block.dt * fluids[cell].Lw1[stage][k];
					}
				}
				if (stage == 1)
				{
					for (int k = 0; k < 4; k++)
					{
						fluids[cell].convar[k] = 0.75 * fluids[cell].convar_old[k] + 0.25 * fluids[cell].convar[k] + 0.25 * block.dt * fluids[cell].Lw1[stage][k];
					}
				}
				if (stage == 2)
				{
					for (int k = 0; k < 4; k++)
					{
						fluids[cell].convar[k] = fluids[cell].convar_old[k] / 3 + 2 * fluids[cell].convar[k] / 3 + 2 * block.dt * fluids[cell].Lw1[stage][k] / 3;
					}
				}
			}
			// RK4
			if(timecoe_list_2d == RK4_2D)
			{
				if (stage == 0)
			{
				for (int k = 0; k < 4; k++)
				{
					fluids[cell].convar[k] = fluids[cell].convar_old[k] + 0.5 * block.dt * fluids[cell].Lw1[stage][k];
				}
			}
				if (stage == 1)
			{
				for (int k = 0; k < 4; k++)
				{
					fluids[cell].convar[k] = fluids[cell].convar_old[k] + 0.5 * block.dt * fluids[cell].Lw1[stage][k];
				}
			}
				if (stage == 2)
			{
				for (int k = 0; k < 4; k++)
				{
					fluids[cell].convar[k] = fluids[cell].convar_old[k] + block.dt * fluids[cell].Lw1[stage][k];
				}
			}
				if (stage == 3)
			{
				for (int k = 0; k < 4; k++)
				{
					fluids[cell].convar[k] = fluids[cell].convar_old[k] 
						+ block.dt * fluids[cell].Lw1[0][k] / 6.0
						+ block.dt * fluids[cell].Lw1[1][k] / 3.0
						+ block.dt * fluids[cell].Lw1[2][k] / 3.0
						+ block.dt * fluids[cell].Lw1[3][k] / 6.0;
				}
			}
			}
		}
	}
}

// three-dimensional problem
TimeMarchingCoefficient_3d timecoe_list_3d = S1O1_3D;

void S1O1_3D(Block3d &block)
{
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient_hweno[0][0][0] = 1.0;
}

void S1O2_3D(Block3d &block)
{
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[0][0][1] = 0.5;
	block.timecoefficient_hweno[0][0][0] = 1.0;
}

void S1O3_3D(Block3d &block)
{
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[0][0][1] = 0.5;
	block.timecoefficient[0][0][2] = 1.0 / 6.0;
}

void S2O4_3D(Block3d &block)
{
	block.stages = 2;
	block.timecoefficient[0][0][0] = 0.5;
	block.timecoefficient[0][0][1] = 1.0 / 8.0;
	block.timecoefficient[1][0][0] = 1.0;
	block.timecoefficient[1][1][0] = 0.0;
	block.timecoefficient[1][0][1] = 1.0 / 6.0;
	block.timecoefficient[1][1][1] = 1.0 / 3.0;

	block.timecoefficient_hweno[0][0][0] = 0.5;
	block.timecoefficient_hweno[1][0][0] = 0.0;
	block.timecoefficient_hweno[1][1][0] = 1.0;
}

void RK2_3D(Block3d &block)
{
	block.stages = 2;
	block.timecoefficient[0][0][0] = 0.5;
	block.timecoefficient[1][1][0] = 1.0;
}

void RK3_3D(Block3d& block)
{
	block.stages = 3;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[1][0][0] = 0.25;
	block.timecoefficient[1][1][0] = 0.25;
	block.timecoefficient[2][0][0] = 1.0 / 6.0;
	block.timecoefficient[2][1][0] = 1.0 / 6.0;
	block.timecoefficient[2][2][0] = 2.0 / 3.0;
}

void RK4_3D(Block3d &block)
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

void Initial_stages(Block3d &block)
{
	for (int i = 0; i < 5; i++) //refers the n stage
	{
		for (int j = 0; j < 5; j++) //refers the nth coefficient at n stage 
		{
			for (int k = 0; k < 3; k++) //refers f, derf, der2f
			{
				block.timecoefficient[i][j][k] = 0.0;
				block.timecoefficient_hweno[i][j][k] = 0.0;
			}
		}
	}
	timecoe_list_3d(block);
}

double Get_CFL(Block3d &block, Fluid3d *fluids, double tstop)
{
	//double dt = block.dx;
	double dt = fluids[block.ghost * block.ny * block.nz
		+ block.ghost * block.nz + block.ghost].R_c;

#pragma omp parallel
		{
			double dtPerThread = dt;
#pragma omp for 
			for (int i = block.ghost; i < block.nx - block.ghost; i++)
			{
				for (int j = block.ghost; j < block.ny - block.ghost; j++)
				{
					for (int k = block.ghost; k < block.nz - block.ghost; k++)
					{
						int index = i * (block.ny * block.nz) + j * block.nz + k;
						dtPerThread = Dtx_3d(dtPerThread, fluids[index].R_c, block.CFL, fluids[index].primvar[0], fluids[index].primvar[1], fluids[index].primvar[2], fluids[index].primvar[3], fluids[index].primvar[4]);
					}
				}
			}
#pragma omp critical 
			{
				if (dtPerThread < dt) {
					dt = dtPerThread;
				}
			}
		}

	if (block.t + dt>tstop)
	{
		dt = tstop - block.t + 1e-16;
	}
	//print time step information
	cout << "step = " << block.step
		<< "time size is " << dt
		<< " time = " << block.t << endl;
	return dt;
}

double Dtx_3d(double dtx, double dx, double CFL, double density, double u, double v,double w, double pressure)
{
	double overden =1.0/ density;
	double tmp;
	tmp = sqrt(u*u + v*v+w*w) + sqrt(Gamma*pressure* overden);
	if (tmp>CFL*dx / dtx)
	{
		dtx = CFL*dx / tmp;
	}

	if (dtx > 0.33*CFL*dx*dx * density / Mu &&tau_type !=Euler&&Mu>0)
	{
		dtx = 0.33*CFL*dx*dx * density / Mu;
	}
	return dtx;
}

void Update_gauss(Fluid3d *fluids, Flux3d_gauss** xfluxes, Flux3d_gauss** yfluxes, Flux3d_gauss** zfluxes, Block3d block, int stage)
{
	if (stage > block.stages)
	{
		cout << "wrong middle stage ,pls check the time marching setting" << endl;
		exit(0);
	}

	double dt = block.dt;
	bool is_gausspoint = true;

	if (gausspoint == 0)
	{
		gausspoint = 1; is_gausspoint = false;
	}

#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost + 1; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			for (int k = block.ghost; k < block.nodez + block.ghost; k++)
			{
				int index = i * (block.ny + 1)*(block.nz + 1) + j * (block.nz + 1) + k;
				for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
				{
					
					for (int var = 0; var < 5; var++)
					{
						double Flux = 0.0;
						for (int current_stage = 0; current_stage < stage + 1; current_stage++)
						{
							Flux = Flux +
								(block.timecoefficient[stage][current_stage][0] * xfluxes[index][current_stage].gauss[num_gauss].f[var]
									+ block.timecoefficient[stage][current_stage][1] * xfluxes[index][current_stage].gauss[num_gauss].derf[var]);

						}
						xfluxes[index][stage].gauss[num_gauss].x[var] = Flux;
					}
				}
			}
		}
	}


#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost + 1; j++)
		{
			for (int k = block.ghost; k < block.nodez + block.ghost; k++)
			{
				int index = i * (block.ny + 1)*(block.nz + 1) + j * (block.nz + 1) + k;
				for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
				{
					for (int var = 0; var < 5; var++)
					{
						double Flux = 0.0;
						for (int current_stage = 0; current_stage < stage + 1; current_stage++)
						{
							Flux = Flux +
								(block.timecoefficient[stage][current_stage][0] * yfluxes[index][current_stage].gauss[num_gauss].f[var]
									+ block.timecoefficient[stage][current_stage][1] * yfluxes[index][current_stage].gauss[num_gauss].derf[var]);

						}
						yfluxes[index][stage].gauss[num_gauss].x[var] = Flux;
					}
				}
			}
		}
	}

#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			for (int k = block.ghost; k < block.nodez + block.ghost + 1; k++)
			{
				int index = i * (block.ny + 1)*(block.nz + 1) + j * (block.nz + 1) + k;
				for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
				{
					for (int var = 0; var < 5; var++)
					{
						double Flux = 0.0;
						for (int current_stage = 0; current_stage < stage + 1; current_stage++)
						{
							Flux = Flux +
								(block.timecoefficient[stage][current_stage][0] * zfluxes[index][current_stage].gauss[num_gauss].f[var]
									+ block.timecoefficient[stage][current_stage][1] * zfluxes[index][current_stage].gauss[num_gauss].derf[var]);
						}
						zfluxes[index][stage].gauss[num_gauss].x[var] = Flux;
					}
				}
			}
		}
	}



#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			for (int k = block.ghost; k < block.nodez + block.ghost; k++)
			{
				int cell = i * (block.ny*block.nz) + j * block.nz + k;
				int face = i * (block.ny + 1)*(block.nz + 1) + j * (block.nz + 1) + k;
				int facex = face + (block.ny + 1)*(block.nz + 1);
				int facey = face + (block.nz + 1);
				int facez = face + 1;

				double total_flux[5];
				for (int var = 0; var < 5; var++)
				{
					fluids[cell].convar[var] = fluids[cell].convar_old[var];
					total_flux[var] = 0.0;
				}

				for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
				{
					double weight = gauss_weight[num_gauss];
					total_flux[0] += weight * block.overdz*(zfluxes[face][stage].gauss[num_gauss].x[0] - zfluxes[facez][stage].gauss[num_gauss].x[0]);
					total_flux[1] += weight * block.overdz*(-zfluxes[face][stage].gauss[num_gauss].x[3] + zfluxes[facez][stage].gauss[num_gauss].x[3]);
					total_flux[2] += weight * block.overdz*(zfluxes[face][stage].gauss[num_gauss].x[2] - zfluxes[facez][stage].gauss[num_gauss].x[2]);
					total_flux[3] += weight * block.overdz*(zfluxes[face][stage].gauss[num_gauss].x[1] - zfluxes[facez][stage].gauss[num_gauss].x[1]);
					total_flux[4] += weight * block.overdz*(zfluxes[face][stage].gauss[num_gauss].x[4] - zfluxes[facez][stage].gauss[num_gauss].x[4]);

					////int overdy = gauss_weight[0] / block.dy;
					total_flux[0] += weight * block.overdy*(yfluxes[face][stage].gauss[num_gauss].x[0] - yfluxes[facey][stage].gauss[num_gauss].x[0]);
					total_flux[1] += weight * block.overdy*(-yfluxes[face][stage].gauss[num_gauss].x[2] + yfluxes[facey][stage].gauss[num_gauss].x[2]);
					total_flux[2] += weight * block.overdy*(yfluxes[face][stage].gauss[num_gauss].x[1] - yfluxes[facey][stage].gauss[num_gauss].x[1]);
					total_flux[3] += weight * block.overdy*(yfluxes[face][stage].gauss[num_gauss].x[3] - yfluxes[facey][stage].gauss[num_gauss].x[3]);
					total_flux[4] += weight * block.overdy*(yfluxes[face][stage].gauss[num_gauss].x[4] - yfluxes[facey][stage].gauss[num_gauss].x[4]);

					//int overdx = gauss_weight[0] / block.dx;
					for (int var = 0; var < 5; var++)
					{
						total_flux[var] +=
							weight * block.overdx*(xfluxes[face][stage].gauss[num_gauss].x[var] - xfluxes[facex][stage].gauss[num_gauss].x[var]);
					}
				}
				for (int var = 0; var < 5; var++)
				{
					fluids[cell].convar[var] += total_flux[var];
				}
			}
		}
	}
}