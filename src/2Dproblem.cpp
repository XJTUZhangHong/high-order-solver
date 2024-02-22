#include "2Dproblem.h"

void PlanarShock()
{
	Runtime runtime;
	runtime.start_initial = clock();
	Block2d block;
	block.uniform = true;
	block.nodex = 500;
	block.nodey = 500;
	block.ghost = 4;

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
				cellreconstruction_2D_normal = WENO7_AO_with_df_normal;
				cellreconstruction_2D_tangent = WENO7_AO_with_df_tangent;
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
	block.nodex = 200;
	block.nodey = 200;
	block.ghost = 4;

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
	gausspoint = 4;
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

			if (is_using_df_factor)
			{
				if (block.step == 0)
				{
					cellreconstruction_2D_normal = First_order_normal;
					cellreconstruction_2D_tangent = First_order_tangent;
				}
				else
				{
					cellreconstruction_2D_normal = WENO7_AO_with_df_normal;
					cellreconstruction_2D_tangent = WENO7_AO_with_df_tangent;
				}
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
	block.ghost = 4;



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
	flux_function_2d = GKS2D;

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
				cellreconstruction_2D_normal = WENO7_AO_with_df_normal;
				cellreconstruction_2D_tangent = WENO7_AO_with_df_tangent;
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
	block.ghost = 4;

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
			cellreconstruction_2D_normal = WENO7_AO_with_df_normal;
			cellreconstruction_2D_tangent = WENO7_AO_with_df_tangent;
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

void viscous_sod_shock_problem()
{
	Runtime runtime;
	runtime.start_initial = omp_get_wtime();
	Block2d block;
	block.uniform = true;
	block.nodex = 1000;
	block.nodey = 500;
	block.ghost = 3;



	block.CFL = 0.5;
	Fluid2d* bcvalue = new Fluid2d[4];

	K = 3;
	Gamma = 1.4;
	double Renum = 1000;
	Mu = 1.0 / Renum;


	//this part should rewritten ad gks2dsolver blabla
	gks2dsolver = gks2nd_2d;
	tau_type = NS;
	c1_euler = 0.01;
	c2_euler = 5;


	//prepare the boundary condtion function
	BoundaryCondition2d leftboundary(0);
	BoundaryCondition2d rightboundary(0);
	BoundaryCondition2d downboundary(0);
	BoundaryCondition2d upboundary(0);

	for (int i = 0; i < 3; i++)
	{
		bcvalue[0].primvar[i] = 0.0;
		bcvalue[1].primvar[i] = 0.0;
		bcvalue[2].primvar[i] = 0.0;
		bcvalue[3].primvar[i] = 0.0;
	}
	leftboundary = noslip_adiabatic_boundary_left;
	rightboundary = noslip_adiabatic_boundary_right;
	downboundary = noslip_adiabatic_boundary_down;
	upboundary = reflection_boundary_up;


	gausspoint = 2;
	SetGuassPoint();

	reconstruction_variable = characteristic;
	wenotype = wenoz;

	cellreconstruction_2D_normal = WENO5_AO_normal;
	cellreconstruction_2D_tangent = WENO5_AO_tangent;
	g0reconstruction_2D_normal = Center_do_nothing_normal;
	g0reconstruction_2D_tangent = Center_all_collision_multi;

	is_reduce_order_warning = true;
	flux_function_2d = GKS2D;


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
	block.left = 0.0;
	block.right = 1.0;
	block.down = 0.0;
	block.up = 0.5;
	block.dx = (block.right - block.left) / block.nodex;
	block.dy = (block.up - block.down) / block.nodey;
	block.overdx = 1 / block.dx;
	block.overdy = 1 / block.dy;
	//set the uniform geometry information
	SetUniformMesh(block, fluids, xinterfaces, yinterfaces, xfluxes, yfluxes);
	//ended mesh part
	double tstop[] = { 1.0 };

	IC_for_vishocktube(block.ghost, fluids, block);

	runtime.finish_initial = omp_get_wtime();
	//then we are about to do the loop
	block.t = 0;//the current simulation time
	block.step = 0; //the current step
	int tstop_dim = 1;

	int inputstep = 1;//input a certain step,
	//initialize inputstep=1, to avoid a 0 result
	for (int instant = 0; instant < tstop_dim; instant++)
	{
		if (inputstep == 0)
		{
			break;
		}

		while (block.t < tstop[instant])
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
			block.dt = Get_CFL(block, fluids, tstop[instant]);

			if (block.step > 00 && is_using_df_factor)
			{
				cellreconstruction_2D_normal = WENO5_AO_with_df_normal;
				cellreconstruction_2D_tangent = WENO5_AO_with_df_tangent;
			}
			for (int i = 0; i < block.stages; i++)
			{
				//after determine the cfl condition, let's implement boundary condtion
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
	}
	output2d(fluids, block);
	runtime.finish_compute = omp_get_wtime();
	cout << "the total run time is " << (double)(runtime.finish_compute - runtime.start_initial) / 1.0 << " second !" << endl;
	cout << "initializing time is" << (double)(runtime.finish_initial - runtime.start_initial) / 1.0 << " second !" << endl;
	cout << "the pure computational time is" << (double)(runtime.finish_compute - runtime.start_compute) / 1.0 << " second !" << endl;

}

void IC_for_vishocktube(int order, Fluid2d* fluid, Block2d block)
{
	int i, j;
	for (i = order; i < block.nx - order; i++)
	{
		for (j = order; j < block.ny - order; j++)
		{
			if (i < 0.5 * block.nx)
			{
				fluid[i * block.ny + j].primvar[0] = 120.0;
				fluid[i * block.ny + j].primvar[1] = 0.0;
				fluid[i * block.ny + j].primvar[2] = 0.0;
				fluid[i * block.ny + j].primvar[3] = 120.0 / Gamma;

			}
			else
			{
				fluid[i * block.ny + j].primvar[0] = 1.2;
				fluid[i * block.ny + j].primvar[1] = 0.0;
				fluid[i * block.ny + j].primvar[2] = 0.0;
				fluid[i * block.ny + j].primvar[3] = 1.2 / Gamma;
			}

		}

	}
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

void accuracy_sinwave_2d()
{
	int mesh_set = 4;
	int mesh_number_start = 10;
	double length = 2.0;
	double CFL = 0.05;
	double dt_ratio = 1.0;

	double mesh_size_start = length / mesh_number_start;
	int* mesh_number = new int[mesh_set];
	double* mesh_size = new double[mesh_set];
	double** error = new double* [mesh_set];

	for (int i = 0; i < mesh_set; ++i)
	{
		error[i] = new double[3];
		mesh_size[i] = mesh_size_start / pow(2, i);
		mesh_number[i] = mesh_number_start * pow(2, i);
		sinwave_2d(CFL, dt_ratio, mesh_number[i], error[i]);
	}

	output_error_form(CFL, dt_ratio, mesh_set, mesh_number, error);
}

void sinwave_2d(double& CFL, double& dt_ratio, int& mesh_number, double* error)
{
	Runtime runtime;
	runtime.start_initial = clock();
	Block2d block; // Class, geometry variables
	block.uniform = true;
	block.nodex = mesh_number;
	block.nodey = mesh_number;
	block.ghost = 4; // 5th-order reconstruction, should 3 ghost cell

	double tstop = 2;
	block.CFL = CFL;


	K = 3; // 1d K=4, 2d K=3, 3d K=2;
	Gamma = 1.4; // diatomic gas, r=1.4

	//this part should rewritten ad gks2dsolver blabla
	gks2dsolver = gks2nd_2d; // Emumeration, choose solver type
	// 2nd means spatical 2nd-order, the traditional GKS; 2d means two dimensions

	tau_type = Euler; // Emumeration, choose collision time tau type
	// for accuracy test case, the tau type should be Euler and more accurately ZERO
	// for inviscid case, use Euler type tau; for viscous case, use NS type tau

	c1_euler = 0.0; // first coefficient in Euler type tau calculation
	c2_euler = 0.0; // second coefficient in Euler type tau calculation
	// when Euler type, tau is always zero;
	// if the c1 c2 were given, (c1, c2, should be given together), the tau_num is determined by c1 c2, (c1*dt + c2*deltaP)
	// if c1, c2 are both zero, (or not given value either), tau_num is also zero, which might be useful for accuracy test.
	// when NS type, tau is always physical tau, relating to Mu or Nu;
	// besides, tau_num would be determined by tau and c2 part, (tau + c2*deltaP),
	// so, c2 should be given (c1 is useless now), and Mu or Nu should be given ONLY ONE;
	// for NS type, specially, if Smooth is true, tau_num is determined ONLY by tau, (tau_num = tau), which means c2 is zero

	Mu = 0.0;
	Pr = 0.73;
	R_gas = 1;

	//prepare the boundary condtion function
	Fluid2d* bcvalue = new Fluid2d[4]; // Class, primitive variables only for boundary
	BoundaryCondition2d leftboundary(0); 
	BoundaryCondition2d rightboundary(0); 
	BoundaryCondition2d downboundary(0);  
	BoundaryCondition2d upboundary(0); 

	leftboundary = periodic_boundary_left; 
	rightboundary = periodic_boundary_right;  
	downboundary = periodic_boundary_down; 
	upboundary = periodic_boundary_up;  

	//prepare the reconstruction
	gausspoint = 4; // fifth-order or sixth-order use THREE gauss points
	// WENO5 has the function relating to arbitrary gausspoints
	// WENO5_AO supports 2 gausspoint now, so fourth-order at most for spacial reconstruction (enough for two step fourth-order GKS)
	SetGuassPoint(); // Function, set Gauss points coordinates and weight factor

	reconstruction_variable = conservative; // Emumeration, choose the variables used for reconstruction type
	wenotype = wenoz; // Emumeration, choose reconstruction type

	cellreconstruction_2D_normal = WENO7_AO_with_df_normal;  // reconstruction in normal directon
	cellreconstruction_2D_tangent = WENO7_AO_with_df_tangent;  // reconstruction in tangential directon
	g0reconstruction_2D_normal = Center_do_nothing_normal;  // reconstruction for g0 in normal directon
	g0reconstruction_2D_tangent = Center_all_collision_multi;  // reconstruction for g0 in tangential directon


	//prepare the flux function
	flux_function_2d = GKS2D; // 给函数指针 赋值 flux calculation type
	// solver is GKS
	// GKS2D means full solver (used for shock), GKS2D_smooth means smooth solver (used for non_shock)
	// for accuracy test case, tau is ZERO, so GKS2D is equal to GKS2D_smooth for results, but the later is faster
	// for inviscid flow, Euler type tau equals artifical viscosity;
	// for viscous flow, NS type tau equals the modified collision time;
	//end

	//prepare time marching stratedgy
	//time coe list must be 2d, end by _2D
	timecoe_list_2d = S2O4_2D; 
	Initial_stages(block);


	// allocate memory for 2-D fluid field
	// in a standard finite element method, we have
	// first the cell average value, N*N
	block.nx = block.nodex + 2 * block.ghost;
	block.ny = block.nodey + 2 * block.ghost;

	Fluid2d* fluids = Setfluid(block); // Function, input a class (geometry), output the pointer of one class (conservative variables)
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
	SetUniformMesh(block, fluids, xinterfaces, yinterfaces, xfluxes, yfluxes); // Function, set mesh

	//end

	ICfor_sinwave_2d(fluids, block);

	runtime.finish_initial = clock();
	block.t = 0; //the current simulation time
	block.step = 0; //the current step

	int inputstep = 1;//input a certain step

	while (block.t < tstop)
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

		CopyFluid_new_to_old(fluids, block); // Function, copy variables
		//determine the cfl condtion
		block.dt = Get_CFL(block, fluids, tstop); // Function, get real CFL number

		if ((block.t + block.dt - tstop) > 0)
		{
			block.dt = tstop - block.t + 1e-16;
		}

		for (int i = 0; i < block.stages; i++)
		{
			//after determine the cfl condition, let's implement boundary condtion
			leftboundary(fluids, block, bcvalue[0]);
			rightboundary(fluids, block, bcvalue[1]);
			downboundary(fluids, block, bcvalue[2]);
			upboundary(fluids, block, bcvalue[3]);
			//cout << "after bc " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;

			Convar_to_Primvar(fluids, block); // Function
			//then is reconstruction part, which we separate the left or right reconstrction
			//and the center reconstruction
			Reconstruction_within_cell(xinterfaces, yinterfaces, fluids, block); // Function
			//cout << "after cell recon " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;

			Reconstruction_forg0(xinterfaces, yinterfaces, fluids, block); // Function
			//cout << "after g0 recon " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;

			//then is solver part
			Calculate_flux(xfluxes, yfluxes, xinterfaces, yinterfaces, block, i); // Function
			//cout << "after flux calcu " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;

			//then is update flux part
			Update(fluids, xfluxes, yfluxes, block, i); // Function

			if (is_using_df_factor)
			{
				Update_alpha(xinterfaces, yinterfaces, fluids, block);
			}
		}
		block.step++;
		block.t = block.t + block.dt;
	}

	error_for_sinwave_2d(fluids, block, tstop, error); // Function, get error
	runtime.finish_compute = clock();
}

void ICfor_sinwave_2d(Fluid2d* fluids, Block2d block)
{
	for (int i = block.ghost; i < block.nx - block.ghost; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost; j++)
		{
			double pi = 3.14159265358979323846;
			int index = i * (block.ny) + j;
			double xleft = (i - block.ghost) * block.dx;
			double xright = (i + 1 - block.ghost) * block.dx;
			double yleft = (j - block.ghost) * block.dy;
			double yright = (j + 1 - block.ghost) * block.dy;

			//case in two dimensional x-y-plane
			double k1 = sin(pi * (xright + yright));
			double k2 = sin(pi * (xright + yleft));
			double k3 = sin(pi * (xleft + yright));
			double k4 = sin(pi * (xleft + yleft));
			fluids[index].primvar[0] = 1.0 - 0.2 / pi / pi / block.dx / block.dy * ((k1 - k2) - (k3 - k4));
			fluids[index].exact = fluids[index].primvar[0];
			fluids[index].primvar[1] = 1;
			fluids[index].primvar[2] = 1;
			fluids[index].primvar[3] = 1;
		}
	}
#pragma omp parallel for
	for (int i = block.ghost; i < block.nx - block.ghost; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost; j++)
		{
			int index = i * (block.ny) + j;
			Primvar_to_convar_2D(fluids[index].convar, fluids[index].primvar);
		}
	}
}

void error_for_sinwave_2d(Fluid2d* fluids, Block2d block, double tstop, double* error)
{
	if (abs(block.t - tstop) <= (1e-10))
	{
		cout << "Accuracy-residual-file-output" << endl;
		double error1 = 0;
		double error2 = 0;
		double error3 = 0;
		for (int i = block.ghost; i < block.nx - block.ghost; i++)
		{
			for (int j = block.ghost; j < block.ny - block.ghost; j++)
			{
				double pi = 3.14159265358979323846;
				int index = i * block.ny + j;
				double xleft = (i - block.ghost) * block.dx;
				double xright = (i + 1 - block.ghost) * block.dx;
				double yleft = (j - block.ghost) * block.dy;
				double yright = (j + 1 - block.ghost) * block.dy;

				//case in two dimensional x-y-plane
				double k1 = sin(pi * (xright + yright));
				double k2 = sin(pi * (xright + yleft));
				double k3 = sin(pi * (xleft + yright));
				double k4 = sin(pi * (xleft + yleft));
				double primvar0 = 1.0 - 0.2 / pi / pi / block.dx / block.dy * ((k1 - k2) - (k3 - k4));


				error1 = error1 + abs(fluids[index].convar[0] - primvar0);
				error2 = error2 + pow(fluids[index].convar[0] - primvar0, 2);
				error3 = (error3 > abs(fluids[index].convar[0] - primvar0)) ? error3 : abs(fluids[index].convar[0] - primvar0);
			}
		}
		error1 /= (block.nodex * block.nodey);
		error2 = sqrt(error2 / (block.nodex * block.nodey));
		error[0] = error1; error[1] = error2; error[2] = error3;
		cout << error1 << endl;
	}
}

