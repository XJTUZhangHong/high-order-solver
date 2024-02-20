#include "2Dproblem.h"

void PlanarShock()
{
	Runtime runtime;
	runtime.start_initial = clock();
	Block2d block;
	block.uniform = true;
	block.nodex = 500;
	block.nodey = 500;
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

			if (block.step > 0 && is_using_df_factor)
			{
				cellreconstruction_2D_normal = WENO5_AO_with_df_normal;
				cellreconstruction_2D_tangent = WENO5_AO_with_df_tangent;
			}
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

				if (is_using_df_factor)
				{
					Update_alpha(xinterfaces, yinterfaces, fluids, block);
				}
			}
			block.step++;
			block.t = block.t + block.dt;
			if (block.step % 100 == 0) { output2d(fluids, block); }
		}
		output2d(fluids, block);
	}
	runtime.finish_compute = clock();
	cout << "the total run time is " << (double)(runtime.finish_compute - runtime.start_initial) / CLOCKS_PER_SEC << " second !" << endl;
	cout << "initializing time is" << (double)(runtime.finish_initial - runtime.start_initial) / CLOCKS_PER_SEC << " second !" << endl;
	cout << "the pure computational time is" << (double)(runtime.finish_compute - runtime.start_compute) / CLOCKS_PER_SEC << " second !" << endl;
}

void PlanarSheer()
{
	Runtime runtime;
	runtime.start_initial = clock();
	Block2d block;
	block.uniform = true;
	block.nodex = 500;
	block.nodey = 500;
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
	flux_function_2d = LF2D;

	//prepare time marching stratedgy
	//time coe list must be 2d
	timecoe_list_2d = RK3_2D;
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
	block.right = 2.0;
	block.down = 0.0;
	block.up = 2.0;
	block.dx = (block.right - block.left) / block.nodex;
	block.dy = (block.up - block.down) / block.nodey;
	block.overdx = 1 / block.dx;
	block.overdy = 1 / block.dy;
	//set the uniform geometry information
	SetUniformMesh(block, fluids, xinterfaces, yinterfaces, xfluxes, yfluxes);
	//ended mesh part

	double tstop[] = { 1.6 };
	double p0 = 1.0; // this p0 is for adjust the speed of shear layer
	double zone1[]{ 1, -0.75, 0.5, p0 };
	double zone2[]{ 3.0, -0.75, -0.5, p0 };
	double zone3[]{ 1.0, 0.75, -0.5, p0 };
	double zone4[]{ 2.0, 0.75, 0.5, p0 };

    double discon[]{1.0, 1.0};
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

			if (block.step > 0 && is_using_df_factor)
			{
				cellreconstruction_2D_normal = WENO5_AO_with_df_normal;
				cellreconstruction_2D_tangent = WENO5_AO_with_df_tangent;
			}
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

				if (is_using_df_factor)
				{
					Update_alpha(xinterfaces, yinterfaces, fluids, block);
				}
			}
			block.step++;
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

void High_mach_astrophusical_jet()
{
	Runtime runtime;
	runtime.start_initial = clock();
	Block2d block;
	block.uniform = true;
	block.nodex = 256;
	block.nodey = 256;
	block.ghost = 3;

	block.CFL = 0.5;
	Fluid2d* bcvalue = new Fluid2d[4];

	Gamma = 5.0 / 3.0;
	K = (int)((5.0 - 3.0 * Gamma) / (Gamma - 1.0) + 1.0);

	//prepare the boundary condtion function
	BoundaryCondition2d leftboundary(0);
	BoundaryCondition2d rightboundary(0);
	BoundaryCondition2d downboundary(0);
	BoundaryCondition2d upboundary(0);

	//leftboundary = free_boundary_left;
	rightboundary = free_boundary_right;
	downboundary = periodic_boundary_down;
	upboundary = periodic_boundary_up;

	//prepare the reconstruction function
	gausspoint = 2;
	SetGuassPoint();

	reconstruction_variable = characteristic;
	wenotype = wenoz;

	// First order tangent reconstruction corresponding to the positive-preserving gks
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
	double tstop[]{ 0.0002 };
	IC_for_astrophusical_jet(fluids, block);

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
			block.dt = 2e-7;
			if (block.step > 0 && is_using_df_factor)
			{
				cellreconstruction_2D_normal = WENO5_AO_with_df_normal;
				cellreconstruction_2D_tangent = WENO5_AO_with_df_tangent;
			}
			for (int i = 0; i < block.stages; i++)
			{
				//leftboundary(fluids, block, bcvalue[0]);
				inflow_boundary_left(fluids, block);
				rightboundary(fluids, block, bcvalue[1]);
				downboundary(fluids, block, bcvalue[2]);
				upboundary(fluids, block, bcvalue[3]);

				Convar_to_Primvar(fluids, block);

				Reconstruction_within_cell(xinterfaces, yinterfaces, fluids, block);

				Reconstruction_forg0(xinterfaces, yinterfaces, fluids, block);

				Calculate_flux(xfluxes, yfluxes, xinterfaces, yinterfaces, block, i);

				Update(fluids, xfluxes, yfluxes, block, i);

				if (is_using_df_factor)
				{
					Update_alpha(xinterfaces, yinterfaces, fluids, block);
				}
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

void IC_for_astrophusical_jet(Fluid2d* fluid, Block2d block)
{	
    int n = block.nx;
	for (int i = 0; i < block.nx; i++)
	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.nx; j++)
		{
			fluid[i * block.nx + j].primvar[0] = 0.5;
			fluid[i * block.nx + j].primvar[1] = 0.0;
			fluid[i * block.nx + j].primvar[2] = 0.0;
			fluid[i * block.nx + j].primvar[3] = 0.4127;
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

void inflow_boundary_left(Fluid2d* fluids, Block2d block)
{
	for (int i = block.ghost - 1; i >= 0; i--)
	{
		for (int j = 0; j < block.ny; j++)
		{
			if (j * block.dy > 0.45 && j * block.dy < 0.55)
			{
				fluids[i * block.nx + j].primvar[0] = 5.0;
				fluids[i * block.nx + j].primvar[1] = 4000.0;
				fluids[i * block.nx + j].primvar[2] = 0.0;
				fluids[i * block.nx + j].primvar[3] = 0.4127;
			}
			else
			{
				fluids[i * block.nx + j].primvar[0] = fluids[(i + 1) * block.nx + j].primvar[0];
				fluids[i * block.nx + j].primvar[1] = fluids[(i + 1) * block.nx + j].primvar[1];
				fluids[i * block.nx + j].primvar[2] = fluids[(i + 1) * block.nx + j].primvar[2];
				fluids[i * block.nx + j].primvar[3] = fluids[(i + 1) * block.nx + j].primvar[3];
			}

			Primvar_to_convar_2D(fluids[i * block.ny + j].convar, fluids[i * block.ny + j].primvar);
		}
	}
}

void RT_instability()
{
	Runtime runtime;
	runtime.start_initial = clock();
	Block2d block;
	block.uniform = true;
	block.nodex = 64;
	block.nodey = 256;
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
	c2_euler = 1;
	flux_function_2d = LF2D;

	//prepare time marching stratedgy


	//time coe list must be 2d
	timecoe_list_2d = RK2_2D;
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
	block.right = 0.25;
	block.down = 0.0;
	block.up = 1.0;
	block.dx = (block.right - block.left) / block.nodex;
	block.dy = (block.up - block.down) / block.nodey;
	block.overdx = 1 / block.dx;
	block.overdy = 1 / block.dy;
	//set the uniform geometry information
	SetUniformMesh(block, fluids, xinterfaces, yinterfaces, xfluxes, yfluxes);
	//ended mesh part

	IC_for_RT_instability(fluids, block);
	double tstop[] = { 1.95 };

	//IC_for_hurricane_problem(fluids, block);
	//double tstop[] = { 0.045 };

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

			if (block.step > 00 && is_using_df_factor)
			{
				cellreconstruction_2D_normal = WENO5_AO_with_df_normal;
				cellreconstruction_2D_tangent = WENO5_AO_with_df_tangent;
			}
			for (int i = 0; i < block.stages; i++)
			{
				RT_boundary(fluids, block);

				Convar_to_Primvar(fluids, block);

				Reconstruction_within_cell(xinterfaces, yinterfaces, fluids, block);

				Reconstruction_forg0(xinterfaces, yinterfaces, fluids, block);

				Calculate_flux(xfluxes, yfluxes, xinterfaces, yinterfaces, block, i);

				Update_RT(fluids, xfluxes, yfluxes, block, i);

				if (is_using_df_factor)
				{
					Update_alpha(xinterfaces, yinterfaces, fluids, block);
				}
			}
			if (block.step % 100 == 0)
			{
				output2d(fluids, block);
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

void IC_for_RT_instability(Fluid2d* fluid, Block2d block)
{
	double c, p, v, pi = acos(-1);
	for (int i = block.ghost; i < block.nx - block.ghost; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost; j++)
		{
			if (j - block.ghost + 1 < 0.5 / block.dy)
			{
				// p = 2y + 1
				p = 2.0 * (j - block.ghost + 0.5) * block.dy + 1.0;
				c = sqrt(Gamma * p / 2.0);
				v = -0.025 * c * cos(8.0 * pi * (i - block.ghost + 0.5) * block.dx);
				fluid[i * block.ny + j].primvar[0] = 2.0;
				fluid[i * block.ny + j].primvar[1] = 0.0;
				fluid[i * block.ny + j].primvar[2] = v;
				fluid[i * block.ny + j].primvar[3] = p;
			}
			else
			{
				p = (j - block.ghost + 0.5) * block.dy + 1.5;
				c = sqrt(Gamma * p / 1.0);
				v = -0.025 * c * cos(8.0 * pi * (i - block.ghost + 0.5) * block.dx);
				fluid[i * block.ny + j].primvar[0] = 1.0;
				fluid[i * block.ny + j].primvar[1] = 0.0;
				fluid[i * block.ny + j].primvar[2] = v;
				fluid[i * block.ny + j].primvar[3] = p;
			}
		}
	}
	int i, j;
	for (i = 0; i < block.nx; i++)
	{
		for (j = 0; j < block.ny; j++)
		{
			fluid[i * block.ny + j].convar[0] = fluid[i * block.ny + j].primvar[0];
			fluid[i * block.ny + j].convar[1] = Q_densityu(fluid[i * block.ny + j].primvar[0], fluid[i * block.ny + j].primvar[1]);
			fluid[i * block.ny + j].convar[2] = Q_densityv(fluid[i * block.ny + j].primvar[0], fluid[i * block.ny + j].primvar[2]);
			fluid[i * block.ny + j].convar[3] = Q_densityE(fluid[i * block.ny + j].primvar[0], fluid[i * block.ny + j].primvar[1], fluid[i * block.ny + j].primvar[2], fluid[i * block.ny + j].primvar[3]);
		}
	}
}

void doubleMach()
{
	Runtime runtime;
	runtime.start_initial = omp_get_wtime();
	Block2d block;
	block.uniform = true;
	block.nodex = 960;
	block.nodey = 240;
	block.ghost = 3;

	double tstop = 0.2;
	block.CFL = 0.5;
	Fluid2d* bcvalue = new Fluid2d[1];

	K = 3;
	Gamma = 1.4;

	//this part should rewritten ad gks2dsolver blabla
	gks2dsolver = gks2nd_2d;
	tau_type = Euler;
	c1_euler = 0.2;
	c2_euler = 1;


	//prepare the boundary condtion function

	//prepare the reconstruction function

	gausspoint = 2;
	SetGuassPoint();

	reconstruction_variable = characteristic;
	wenotype = wenoz;

	cellreconstruction_2D_normal = WENO5_AO_normal;
	cellreconstruction_2D_tangent = WENO5_AO_tangent;
	g0reconstruction_2D_normal = Center_do_nothing_normal;
	g0reconstruction_2D_tangent = Center_all_collision_multi;

	is_reduce_order_warning = false;
	//prepare the flux function
	flux_function_2d = GKS2D;
	//flux_function_2d = GKS2D_with_avg_der;
	//prepare time marching stratedgy
	//block.stages = 1;
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
	//bulid or read mesh,
	//since the mesh is all structured from left to right, down to up
	//there is no need for buliding the complex topology between cells and interfaces
	//just using the index for address searching
	if (block.uniform == true)
	{
		block.left = 0.0;
		block.right = 4.0;
		block.down = 0.0;
		block.up = 1.0;
		block.dx = (block.right - block.left) / block.nodex;
		block.dy = (block.up - block.down) / block.nodey;
		block.overdx = 1 / block.dx;
		block.overdy = 1 / block.dy;
		//set the uniform geometry information
		SetUniformMesh(block, fluids, xinterfaces, yinterfaces, xfluxes, yfluxes);
	}
	else
	{
		//set non-uniform geomertry information, should be added later
	}
	//ended mesh part

	ICforDoubleMach(fluids, block);
	//initializing part end

	runtime.finish_initial = omp_get_wtime();
	//then we are about to do the loop
	block.t = 0;//the current simulation time
	block.step = 0; //the current step

	int inputstep = 1;//input a certain step,
	//initialize inputstep=1, to avoid a 0 result
	while (block.t < tstop)
	{
		// assume you are using command window,
		// you can specify a running step conveniently
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
			runtime.start_compute = omp_get_wtime();
			cout << "runtime-start " << endl;
		}
		//Copy the fluid vales to fluid old
		//and install primvar for it
		CopyFluid_new_to_old(fluids, block);
		//determine the cfl condtion
		block.dt = Get_CFL(block, fluids, tstop);
		if (block.step > 00 && is_using_df_factor)
		{
			cellreconstruction_2D_normal = WENO5_AO_with_df_normal;
			cellreconstruction_2D_tangent = WENO5_AO_with_df_tangent;
		}
		for (int i = 0; i < block.stages; i++)
		{

			boundaryforDoubleMach(fluids, block, bcvalue[0]);

			Convar_to_Primvar(fluids, block);

			Reconstruction_within_cell(xinterfaces, yinterfaces, fluids, block);

			Reconstruction_forg0(xinterfaces, yinterfaces, fluids, block);

			Calculate_flux(xfluxes, yfluxes, xinterfaces, yinterfaces, block, i);

			Update(fluids, xfluxes, yfluxes, block, i);
			if (is_using_df_factor)
			{
				Update_alpha(xinterfaces, yinterfaces, fluids, block);
			}
		}
		block.step++;
		block.t = block.t + block.dt;
	}
	output2d(fluids, block);
	runtime.finish_compute = omp_get_wtime();

	cout << "the total run time is " << (double)(runtime.finish_compute - runtime.start_initial) / CLOCKS_PER_SEC << " second !" << endl;
	cout << "initializing time is" << (double)(runtime.finish_initial - runtime.start_initial) / CLOCKS_PER_SEC << " second !" << endl;
	cout << "the pure computational time is" << (double)(runtime.finish_compute - runtime.start_compute) / CLOCKS_PER_SEC << " second !" << endl;
}

void ICforDoubleMach(Fluid2d* fluid, Block2d block)
{
	for (int i = block.ghost; i < block.nx - block.ghost; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost; j++)
		{
			int nodex = block.nodex;
			int index = i * block.ny + j;

			if (fluid[index].coordx < 1.0 / 6.0)
			{
				fluid[index].primvar[0] = 8.0;
				fluid[index].primvar[1] = 4.125 * sqrt(3);
				fluid[index].primvar[2] = -4.125;
				fluid[index].primvar[3] = 116.5;
			}
			else
			{
				if (fluid[index].coordy > sqrt(3) * (fluid[index].coordx - 1.0 / 6.0))
				{
					fluid[index].primvar[0] = 8.0;
					fluid[index].primvar[1] = 4.125 * sqrt(3);
					fluid[index].primvar[2] = -4.125;
					fluid[index].primvar[3] = 116.5;
				}
				else
				{
					fluid[index].primvar[0] = 1.4;
					fluid[index].primvar[1] = 0;
					fluid[index].primvar[2] = 0;
					fluid[index].primvar[3] = 1;
				}
			}
		}
	}
	// the derivative information shall be zero, in default
	for (int i = block.ghost; i < block.nx - block.ghost; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost; j++)
		{
			int index = i * block.ny + j;
			Primvar_to_convar_2D(fluid[index].convar, fluid[index].primvar);
		}
	}
}

void boundaryforDoubleMach(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	Fluid2d zone1;
	zone1.primvar[0] = 8; zone1.primvar[1] = 4.125 * sqrt(3); zone1.primvar[2] = -4.125; zone1.primvar[3] = 116.5;
	Fluid2d zone2;
	zone2.primvar[0] = 1.4; zone2.primvar[1] = 0, zone2.primvar[2] = 0, zone2.primvar[3] = 1.0;

	inflow_boundary_left(fluids, block, zone1);
	free_boundary_right(fluids, block, bcvalue);

	//boundary up
	double convar1[4];
	Primvar_to_convar_2D(convar1, zone1.primvar);
	double convar2[4];
	Primvar_to_convar_2D(convar2, zone2.primvar);
	int order = block.ghost;
	int nodex = block.nx - 2 * order;
	double height = block.up - block.down;
	for (int j = block.ny - order; j < block.ny; j++)
	{
		for (int i = order; i < block.nx; i++)
		{
			int index = i * block.ny + j;
			if (fluids[index].coordx < (1.0 / 6.0 + height / sqrt(3.0) + 20 / sqrt(3.0) * block.t))
			{
				for (int var = 0; var < 4; var++)
				{
					fluids[index].convar[var] = convar1[var];
				}
			}
			else
			{
				for (int var = 0; var < 4; var++)
				{
					fluids[index].convar[var] = convar2[var];
				}
			}
		}
	}
	////boundary down
	for (int j = order - 1; j >= 0; j--)
	{
		for (int i = order; i < block.nx; i++)
		{
			int index = i * block.ny + j;
			if (fluids[index].coordx < 1.0 / 6.0)
			{
				//cout << i <<" "<< j << endl;
				for (int var = 0; var < 4; var++)
				{
					fluids[index].convar[var] = convar1[var];
				}
			}
			else
			{
				int ref = i * block.ny + 2 * order - 1 - j;
				fluids[index].convar[0] = fluids[ref].convar[0];
				fluids[index].convar[1] = fluids[ref].convar[1];
				fluids[index].convar[2] = -fluids[ref].convar[2];
				fluids[index].convar[3] = fluids[ref].convar[3];
			}
		}

	}
}




