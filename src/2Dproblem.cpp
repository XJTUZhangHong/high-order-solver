#include "2Dproblem.h"

void PlanarShock()
{
	Runtime runtime;
	runtime.start_initial = clock();
	Block2d block;
	block.uniform = true;
	block.nodex = 100;
	block.nodey = 100;
	block.ghost = 3;

	block.CFL = 0.5;
	Fluid2d* bcvalue = new Fluid2d[4];

	K = 3;
	Gamma = 1.4;

	//prepare the boundary condtion function
	BoundaryCondition2d leftboundary(0);
	BoundaryCondition2d rightboundary(0);
	BoundaryCondition2d downboundary(0);
	BoundaryCondition2d upboundary(0);

	leftboundary = free_boundary_left;
	rightboundary = free_boundary_right;
	downboundary = free_boundary_down;
	upboundary = free_boundary_up;

	//prepare the reconstruction function
	gausspoint = 2;
	SetGuassPoint();

	reconstruction_variable = characteristic;
	wenotype = wenoz;

	cellreconstruction_2D_normal = WENO5_AO_normal;
	cellreconstruction_2D_tangent = WENO5_AO_tangent;
	g0reconstruction_2D_normal = Center_do_nothing_normal;
	g0reconstruction_2D_tangent = Center_all_collision_multi;

	is_reduce_order_warning = true; 

	//prepare the flux function
	gks2dsolver = gks2nd_2d;
	tau_type = Euler;
	c1_euler = 0.05;
	c2_euler = 1.0;
	flux_function_2d = GKS2D;

	//prepare time marching stratedgy
	//time coe list must be 2d
	timecoe_list_2d = S2O4_2D;
	Initial_stages(block);


	// allocate memory for 2-D fluid field
	// in a standard finite element method, we have 
	// first the cell average value, N*N
	block.nx = block.nodex + 2 * block.ghost;
	block.ny = block.nodey + 2 * block.ghost;

	Fluid2d* fluids = Setfluid(block);
	// then the interfaces reconstructioned value, (N+1)*(N+1)
	Interface2d* xinterfaces = Setinterface_array(block);
	Interface2d* yinterfaces = Setinterface_array(block);
	// then the flux, which the number is identical to interfaces
	Flux2d_gauss** xfluxes = Setflux_gauss_array(block);
	Flux2d_gauss** yfluxes = Setflux_gauss_array(block);
	//end the allocate memory part

	block.left = 0.0;
	block.right = 1.0;
	block.down = 0.0;
	block.up = 1.0;
	block.dx = (block.right - block.left) / block.nodex;
	block.dy = (block.up - block.down) / block.nodey;
	block.overdx = 1 / block.dx;
	block.overdy = 1 / block.dy;
	//set the uniform geometry information
	SetUniformMesh(block, fluids, xinterfaces, yinterfaces, xfluxes, yfluxes);
	//ended mesh part

	//RM 2 T=0.6 with x,y = 0.7
	double tstop[]{ 0.6 };
	double zone1[]{ 0.138, 1.206, 1.206, 0.029 };
	double zone2[]{ 0.5323, 0, 1.206, 0.3 };
	double zone3[]{ 1.5, 0, 0, 1.5 };
	double zone4[]{ 0.5323, 1.206, 0, 0.3 };

    double discon[]{0.7, 0.7};
	IC_for_riemann_2d(fluids, zone1, zone2, zone3, zone4, block, discon);

	runtime.finish_initial = clock();
	block.t = 0;//the current simulation time
	block.step = 0; //the current step
	int tstop_dim = sizeof(tstop) / sizeof(double);

	int inputstep = 1;//input a certain step,
					  //initialize inputstep=1, to avoid a 0 result

	for (int instant = 0; instant < tstop_dim; ++instant)
	{

		if (inputstep == 0)
		{
			break;
		}
		while (block.t < tstop[instant])
		{

			if (block.step % inputstep == 0)
			{
				cout << "pls cin interation step, if input is 0, then the program will exit " << endl;
				cin >> inputstep;
				if (inputstep == 0)
				{
					output2d(fluids, block);
					break;
				}
			}
			if (runtime.start_compute == 0.0)
			{
				runtime.start_compute = clock();
				cout << "runtime-start " << endl;
			}

			CopyFluid_new_to_old(fluids, block);

			block.dt = Get_CFL(block, fluids, tstop[instant]);

			for (int i = 0; i < block.stages; i++)
			{
				leftboundary(fluids, block, bcvalue[0]);
				rightboundary(fluids, block, bcvalue[1]);
				downboundary(fluids, block, bcvalue[2]);
				upboundary(fluids, block, bcvalue[3]);

				Convar_to_Primvar(fluids, block);

				Reconstruction_within_cell(xinterfaces, yinterfaces, fluids, block);

				Reconstruction_forg0(xinterfaces, yinterfaces, fluids, block);

				Calculate_flux(xfluxes, yfluxes, xinterfaces, yinterfaces, block, i);

				Update(fluids, xfluxes, yfluxes, block, i);
			}
			block.step++;
			//cout << "The step is " << block.step << endl;
			block.t = block.t + block.dt;
		}
		output2d(fluids, block);
	}
	runtime.finish_compute = clock();
	cout << "the total run time is " << (double)(runtime.finish_compute - runtime.start_initial) / CLOCKS_PER_SEC << " second !" << endl;
	cout << "initializing time is" << (double)(runtime.finish_initial - runtime.start_initial) / CLOCKS_PER_SEC << " second !" << endl;
	cout << "the pure computational time is" << (double)(runtime.finish_compute - runtime.start_compute) / CLOCKS_PER_SEC << " second !" << endl;
}

void IC_for_riemann_2d(Fluid2d* fluid, double* zone1, double* zone2, double* zone3, double* zone4, Block2d block, double* discon)
{	
    for (int i = block.ghost; i < block.nx - block.ghost; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost; j++)
		{
			if (i - block.ghost + 1 < discon[0] / block.dx && j - block.ghost + 1 < discon[1] / block.dx)
			{
				Copy_Array(fluid[i * block.ny + j].primvar, zone1, 4);
			}
			else
			{
				if (i - block.ghost + 1 >= discon[0] / block.dx && j - block.ghost + 1 < discon[1] / block.dx)
				{
					Copy_Array(fluid[i * block.ny + j].primvar, zone2, 4);
				}
				else
				{
					if (i - block.ghost + 1 >= discon[0] / block.dx && j - block.ghost + 1 >= discon[1] / block.dx)
					{
						Copy_Array(fluid[i * block.ny + j].primvar, zone3, 4);
					}
					else
					{
						Copy_Array(fluid[i * block.ny + j].primvar, zone4, 4);
					}
				}
			}

		}
	}

	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			Primvar_to_convar_2D(fluid[i * block.ny + j].convar, fluid[i * block.ny + j].primvar);

		}
	}
}

