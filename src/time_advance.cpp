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
					// calculate the final flux of the interface, in x in xfluxes, by the obtained the flux and its derivative (f, def, der2f) at guass points, in xfluxes, and the corresponding weight factors
					// calculate by several stages according to the time marching method. same for yfluxes
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

void Update_RT(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, int stage)
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
				for (int num_gauss = 0; num_gauss < 2; num_gauss++)
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

