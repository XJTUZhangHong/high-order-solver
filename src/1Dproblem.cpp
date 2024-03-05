#include "1Dproblem.h"

void SodTubeProblem()
{
	Runtime runtime;//statement for recording the running time
	runtime.start_initial = clock();

	Block1d block;
	block.nodex = 100;
	block.ghost = 5;

	double tstop = 0.2;
	block.CFL = 0.5;
	Fluid1d* bcvalue = new Fluid1d[2];
	K = 4;
	Gamma = 1.4;
	R_gas = 1.0;

	gks1dsolver = gks2nd;
	//g0type = collisionn;
	tau_type = Euler;
	//Smooth = false;
	c1_euler = 0.05;
	c2_euler = 1.0;

	//prepare the boundary condtion function
	BoundaryCondition leftboundary(0);
	BoundaryCondition rightboundary(0);
	leftboundary = free_boundary_left;
	rightboundary = free_boundary_right;
	//prepare the reconstruction function

	cellreconstruction = WENO5_AO;
	wenotype = wenoz;
	reconstruction_variable = characteristic;
	g0reconstruction = Center_collision;
	is_reduce_order_warning = true;
	//prepare the flux function
	flux_function = GKS;
	//prepare time marching stratedgy
	timecoe_list = S2O4;
	Initial_stages(block);

	// allocate memory for 1-D fluid field
	// in a standard finite element method, we have 
	// first the cell average value, N
	block.nx = block.nodex + 2 * block.ghost;
	block.nodex_begin = block.ghost;
	block.nodex_end = block.nodex + block.ghost - 1;
	Fluid1d* fluids = new Fluid1d[block.nx];
	// then the interfaces reconstructioned value, N+1
	Interface1d* interfaces = new Interface1d[block.nx + 1];
	// then the flux, which the number is identical to interfaces
	Flux1d** fluxes = Setflux_array(block);
	//end the allocate memory part

	//bulid or read mesh,
	//since the mesh is all structured from left and right,
	//there is no need for buliding the complex topology between cells and interfaces
	//just using the index for address searching

	block.left = 0.0;
	block.right = 1.0;
	block.dx = (block.right - block.left) / block.nodex;
	//set the uniform geometry information
	SetUniformMesh(block, fluids, interfaces, fluxes);
	//ended mesh part

	ICforSodTube(fluids, block);
	//initializing part end

	//then we are about to do the loop
	block.t = 0;//the current simulation time
	block.step = 0; //the current step

	int inputstep = 1;//input a certain step,
					  //initialize inputstep=1, to avoid a 0 result
	runtime.finish_initial = clock();
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
				output1d(fluids, block);
				break;
			}
		}
		if (runtime.start_compute == 0.0)
		{
			runtime.start_compute = clock();
			cout << "runtime-start " << endl;
		}
		//Copy the fluid vales to fluid old
		CopyFluid_new_to_old(fluids, block);
		//determine the cfl condtion
		block.dt = Get_CFL(block, fluids, tstop);

		if (block.step > 0 && is_using_df_factor)
		{
			cellreconstruction = WENO9_AO_with_DF;
		}

		for (int i = 0; i < block.stages; i++)
		{
			//after determine the cfl condition, let's implement boundary condtion
			leftboundary(fluids, block, bcvalue[0]);
			rightboundary(fluids, block, bcvalue[1]);
			// here the boudary type, you shall go above the search the key words"BoundaryCondition leftboundary;"
			// to see the pointer to the corresponding function

			//then is reconstruction part, which we separate the left or right reconstrction
			//and the center reconstruction
			Reconstruction_within_cell(interfaces, fluids, block);

			Reconstruction_forg0(interfaces, fluids, block);
			
			//then is solver part
			Calculate_flux(fluxes, interfaces, block, i);
			//then is update flux part
			Update(fluids, fluxes, block, i);

			if (is_using_df_factor)
			{
				Update_alpha(interfaces, fluids, block);
			}
		}
		// update the compression factor
		//for (int j = 3; j < 402; j++) { cout << fluids[j].convar[0] << endl; }
		block.step++;
		block.t = block.t + block.dt;
		if (block.step % 1000 == 0)
		{
			cout << "step 1000 time is " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;
		}
		if ((block.t - tstop) > 0)
		{
			runtime.finish_compute = clock();
			cout << "initializiation time is " << (float)(runtime.finish_initial - runtime.start_initial) / CLOCKS_PER_SEC << "seconds" << endl;
			cout << "computational time is " << (float)(runtime.finish_compute - runtime.start_compute) / CLOCKS_PER_SEC << "seconds" << endl;
			output1d(fluids, block);
			//output1d_checking(fluids, interfaces, fluxes, block);
		}
	}
}

void ICforSodTube(Fluid1d* fluids, Block1d block)
{
	for (int i = block.ghost; i <= block.nodex + 2 * block.ghost - 1; i++)
	{
		if ((double)(i - block.ghost + 1) * block.dx <= 0.5)
		{
			fluids[i].primvar[0] = 1.0;
			fluids[i].primvar[1] = 0.0;
			fluids[i].primvar[2] = 1.0;
		}
		else
		{
			fluids[i].primvar[0] = 0.125;
			fluids[i].primvar[1] = 0.0;
			fluids[i].primvar[2] = 0.1;
		}
		//cout << fluids[i].primvar[0] << endl;
	}

	for (int i = 0; i < block.nx; i++)
	{
		Primvar_to_convar_1D(fluids[i].convar, fluids[i].primvar);
	}
}

void Blastwave()
{
	Runtime runtime;//statement for recording the running time
	runtime.start_initial = clock();

	Block1d block;
	block.nodex = 400;
	block.ghost = 5;

	double tstop = 0.038;
	block.CFL = 0.5;
	Fluid1d* bcvalue = new Fluid1d[2];
	K = 4;
	Gamma = 1.4;
	R_gas = 1.0;

	gks1dsolver = gks2nd;
	//g0type = collisionn;
	tau_type = Euler;
	//Smooth = false;
	c1_euler = 0.01;
	c2_euler = 5;

	//prepare the boundary condtion function
	BoundaryCondition leftboundary(0);
	BoundaryCondition rightboundary(0);
	leftboundary = reflection_boundary_left;
	rightboundary = reflection_boundary_right;
	//prepare the reconstruction function

	cellreconstruction = WENO5_AO;
	wenotype = wenoz;
	reconstruction_variable = characteristic;
	g0reconstruction = Center_collision;
	is_reduce_order_warning = true;
	//prepare the flux function
	flux_function = GKS;
	//prepare time marching stratedgy
	timecoe_list = S2O4;
	Initial_stages(block);

	// allocate memory for 1-D fluid field
	// in a standard finite element method, we have 
	// first the cell average value, N
	block.nx = block.nodex + 2 * block.ghost;
	block.nodex_begin = block.ghost;
	block.nodex_end = block.nodex + block.ghost - 1;
	Fluid1d* fluids = new Fluid1d[block.nx];
	// then the interfaces reconstructioned value, N+1
	Interface1d* interfaces = new Interface1d[block.nx + 1];
	// then the flux, which the number is identical to interfaces
	Flux1d** fluxes = Setflux_array(block);
	//end the allocate memory part

	//bulid or read mesh,
	//since the mesh is all structured from left and right,
	//there is no need for buliding the complex topology between cells and interfaces
	//just using the index for address searching

	block.left = 0.0;
	block.right = 1.0;
	block.dx = (block.right - block.left) / block.nodex;
	//set the uniform geometry information
	SetUniformMesh(block, fluids, interfaces, fluxes);
	//ended mesh part

	ICforBlastwave(fluids, block);
	//initializing part end


	//then we are about to do the loop
	block.t = 0;//the current simulation time
	block.step = 0; //the current step

	int inputstep = 1;//input a certain step,
					  //initialize inputstep=1, to avoid a 0 result
	runtime.finish_initial = clock();
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
				output1d(fluids, block);
				break;
			}
		}
		if (runtime.start_compute == 0.0)
		{
			runtime.start_compute = clock();
			cout << "runtime-start " << endl;
		}
		//Copy the fluid vales to fluid old
		CopyFluid_new_to_old(fluids, block);
		//determine the cfl condtion
		block.dt = Get_CFL(block, fluids, tstop);

		if (block.step > 0 && is_using_df_factor)
		{
			cellreconstruction = WENO5_AO_with_DF;
		}

		for (int i = 0; i < block.stages; i++)
		{
			//after determine the cfl condition, let's implement boundary condtion
			leftboundary(fluids, block, bcvalue[0]);
			rightboundary(fluids, block, bcvalue[1]);
			// here the boudary type, you shall go above the search the key words"BoundaryCondition leftboundary;"
			// to see the pointer to the corresponding function

			//then is reconstruction part, which we separate the left or right reconstrction
			//and the center reconstruction
			Reconstruction_within_cell(interfaces, fluids, block);

			Reconstruction_forg0(interfaces, fluids, block);
			
			//then is solver part
			Calculate_flux(fluxes, interfaces, block, i);
			//then is update flux part
			Update(fluids, fluxes, block, i);

			if (is_using_df_factor)
			{
				Update_alpha(interfaces, fluids, block);
			}
		}
		block.step++;
		block.t = block.t + block.dt;
		if (block.step % 1000 == 0)
		{
			cout << "step 1000 time is " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;
		}
		if ((block.t - tstop) > 0)
		{
			runtime.finish_compute = clock();
			cout << "initializiation time is " << (float)(runtime.finish_initial - runtime.start_initial) / CLOCKS_PER_SEC << "seconds" << endl;
			cout << "computational time is " << (float)(runtime.finish_compute - runtime.start_compute) / CLOCKS_PER_SEC << "seconds" << endl;
			output1d(fluids, block);
		}
	}

}

void ICforBlastwave(Fluid1d* fluids, Block1d block)
{
	for (int i = block.ghost; i <= block.nodex + block.ghost - 1; i++)
	{
		if ((double)(i - block.ghost + 1) * block.dx <= 0.1)
		{
			fluids[i].primvar[0] = 1.0; fluids[i].primvar[1] = 0.0; fluids[i].primvar[2] = 1000;
		}
		else if ((double)(i - block.ghost + 1) * block.dx >= 0.9)
		{
			fluids[i].primvar[0] = 1.0; fluids[i].primvar[1] = 0.0; fluids[i].primvar[2] = 100;
		}
		else
		{
			fluids[i].primvar[0] = 1.0; fluids[i].primvar[1] = 0.0; fluids[i].primvar[2] = 0.01;
		}
	}

	for (int i = 0; i < block.nx; i++)
	{
		Primvar_to_convar_1D(fluids[i].convar, fluids[i].primvar);
	}
}

void ShuOsher()
{
	Runtime runtime;//statement for recording the running time
	runtime.start_initial = clock();

	Block1d block;
	block.nodex = 200;
	block.ghost = 5;

	double tstop = 1.8;
	block.CFL = 0.5;
	Fluid1d* bcvalue = new Fluid1d[2];
	K = 4;
	Gamma = 1.4;
	R_gas = 1.0;

	gks1dsolver = gks2nd;
	//g0type = collisionn;
	tau_type = Euler;
	//Smooth = false;
	c1_euler = 0.05;
	c2_euler = 1.0;

	//prepare the boundary condtion function
	BoundaryCondition leftboundary(0);
	BoundaryCondition rightboundary(0);
	leftboundary = free_boundary_left;
	rightboundary = free_boundary_right;
	//prepare the reconstruction function

	cellreconstruction = WENO5_AO;
	wenotype = wenoz;
	reconstruction_variable = characteristic;
	g0reconstruction = Center_collision;
	is_reduce_order_warning = true;
	//prepare the flux function
	flux_function = GKS;
	//prepare time marching stratedgy
	timecoe_list = S2O4;
	Initial_stages(block);

	// allocate memory for 1-D fluid field
	// in a standard finite element method, we have 
	// first the cell average value, N
	block.nx = block.nodex + 2 * block.ghost;
	block.nodex_begin = block.ghost;
	block.nodex_end = block.nodex + block.ghost - 1;
	Fluid1d* fluids = new Fluid1d[block.nx];
	// then the interfaces reconstructioned value, N+1
	Interface1d* interfaces = new Interface1d[block.nx + 1];
	// then the flux, which the number is identical to interfaces
	Flux1d** fluxes = Setflux_array(block);
	//end the allocate memory part

	//bulid or read mesh,
	//since the mesh is all structured from left and right,
	//there is no need for buliding the complex topology between cells and interfaces
	//just using the index for address searching

	block.left = 0.0;
	block.right = 10.0;
	block.dx = (block.right - block.left) / block.nodex;
	//set the uniform geometry information
	SetUniformMesh(block, fluids, interfaces, fluxes);
	//ended mesh part

	ICforShuOsher(fluids, block);
	//initializing part end

	//then we are about to do the loop
	block.t = 0;//the current simulation time
	block.step = 0; //the current step

	int inputstep = 1;//input a certain step,
					  //initialize inputstep=1, to avoid a 0 result
	runtime.finish_initial = clock();
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
				output1d(fluids, block);
				break;
			}
		}
		if (runtime.start_compute == 0.0)
		{
			runtime.start_compute = clock();
			cout << "runtime-start " << endl;
		}
		//Copy the fluid vales to fluid old
		CopyFluid_new_to_old(fluids, block);
		//determine the cfl condtion
		block.dt = Get_CFL(block, fluids, tstop);

		
		if (block.step > 0 && is_using_df_factor)
		{
			cellreconstruction = WENO9_AO_with_DF;
		}

		for (int i = 0; i < block.stages; i++)
		{
			//after determine the cfl condition, let's implement boundary condtion
			leftboundary(fluids, block, bcvalue[0]);
			//rightboundary(fluids, block, bcvalue[1]);
			// here the boudary type, you shall go above the search the key words"BoundaryCondition leftboundary;"
			// to see the pointer to the corresponding function

			//then is reconstruction part, which we separate the left or right reconstrction
			//and the center reconstruction
			Reconstruction_within_cell(interfaces, fluids, block);

			Reconstruction_forg0(interfaces, fluids, block);
			
			//then is solver part
			Calculate_flux(fluxes, interfaces, block, i);
			//then is update flux part
			Update(fluids, fluxes, block, i);

			if (is_using_df_factor)
			{
				Update_alpha(interfaces, fluids, block);
			}
		}
		block.step++;
		block.t = block.t + block.dt;
		if (block.step % 1000 == 0)
		{
			cout << "step 1000 time is " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;
		}
		if ((block.t - tstop) > 0)
		{
			runtime.finish_compute = clock();
			cout << "initializiation time is " << (float)(runtime.finish_initial - runtime.start_initial) / CLOCKS_PER_SEC << "seconds" << endl;
			cout << "computational time is " << (float)(runtime.finish_compute - runtime.start_compute) / CLOCKS_PER_SEC << "seconds" << endl;
			output1d(fluids, block);
		}
	}
}

void ICforShuOsher(Fluid1d* fluids, Block1d block)
{
	for (int i = block.ghost; i <= block.nodex + 2 * block.ghost - 1; i++)
	{
		double coordx = (i - block.ghost + 0.5) * block.dx;
		if (coordx <= 1.0)
		{
			fluids[i].primvar[0] = 3.85714;
			fluids[i].primvar[1] = 2.629369;
			fluids[i].primvar[2] = 10.33333;
		}
		else
		{
			fluids[i].primvar[0] = 1 - (0.2 / (5 * block.dx)) * (cos(5 * (coordx + 0.5 * block.dx)) - cos(5 * (coordx - 0.5 * block.dx)));
			fluids[i].primvar[1] = 0;
			fluids[i].primvar[2] = 1;
		}
		//cout << fluids[i].primvar[0] << endl;
	}

	for (int i = 0; i < block.nx; i++)
	{
		Primvar_to_convar_1D(fluids[i].convar, fluids[i].primvar);
	}
}

void accuracy_sinwave_1d()
{
	int mesh_set = 4; 
	int mesh_number_start = 10; 
	double length = 2.0; 
	double CFL = 0.02;

	double dt_ratio = 1.0;
	//end

	double mesh_size_start = length / mesh_number_start; 

	int* mesh_number = new int[mesh_set];
	double* mesh_size = new double[mesh_set];
	double** error = new double* [mesh_set];
	//end

	for (int i = 0; i < mesh_set; i++)
	{
		error[i] = new double[3]; 
		mesh_size[i] = mesh_size_start / pow(2, i); 
		mesh_number[i] = mesh_number_start * pow(2, i);

		accuracy_sinwave_1d(CFL, dt_ratio, mesh_number[i], error[i]);
		//end
	}

	output_error_form(CFL, dt_ratio, mesh_set, mesh_number, error);
}

void accuracy_sinwave_1d(double& CFL, double& dt_ratio, int& mesh_number, double* error)
{
	Runtime runtime;
	runtime.start_initial = clock();

	Block1d block; 

	block.nodex = mesh_number;
	block.ghost = 5; 

	double tstop = 2.0; 

	block.CFL = CFL; 



	K = 4; 
	Gamma = 1.4; 

	tau_type = Euler; 


	flux_function = GKS; 

	gks1dsolver = gks2nd; 
	c1_euler = 0.0; 
	c2_euler = 0.0; 
	//end

	Fluid1d* bcvalue = new Fluid1d[2];
	BoundaryCondition leftboundary(0);
	BoundaryCondition rightboundary(0);
	leftboundary = periodic_boundary_left;
	rightboundary = periodic_boundary_right;
	//end


	cellreconstruction = WENO9_AO_with_DF;
	wenotype = wenoz;
	reconstruction_variable = conservative;
	g0reconstruction = Center_collision;
	//end


	timecoe_list = S2O4;
	Initial_stages(block); 
	//end

	// allocate memory for 1-D fluid field
	// in a standard finite element method, we have 
	// first the cell average value, N
	block.nx = block.nodex + 2 * block.ghost;
	block.nodex_begin = block.ghost;
	block.nodex_end = block.nodex + block.ghost - 1;
	Fluid1d* fluids = new Fluid1d[block.nx];

	// then the interfaces reconstructioned value, N+1
	Interface1d* interfaces = new Interface1d[block.nx + 1];
	// then the flux, which the number is identical to interfaces
	Flux1d** fluxes = Setflux_array(block);
	//end the allocate memory part

	//bulid or read mesh,
	//since the mesh is all structured from left and right,
	//there is no need for buliding the complex topology between cells and interfaces
	//just using the index for address searching

	block.left = 0.0;
	block.right = 2.0;
	block.dx = (block.right - block.left) / block.nodex;
	//set the uniform geometry information
	SetUniformMesh(block, fluids, interfaces, fluxes);



	ICfor_sinwave(fluids, block);


	block.t = 0;//the current simulation time
	block.step = 0; //the current step
	int inputstep = 1;//输入一个运行的时间步数

	runtime.finish_initial = clock(); //记录一下初始化需要的cpu时间....

	while (block.t < tstop) //显式时间步循环，from t^n to t^n+1
	{

		if (block.step % inputstep == 0)
		{
			cout << "pls cin interation step, if input is 0, then the program will exit " << endl;
			cin >> inputstep;
			if (inputstep == 0)
			{
				output1d(fluids, block); //这个是输出一维的流场数据
				break;
			}
		}

		if (runtime.start_compute == 0.0)
		{
			runtime.start_compute = clock();
			cout << "runtime-start " << endl;
		}
		//end
		//Copy the fluid vales to fluid old
		CopyFluid_new_to_old(fluids, block);

		//determine the cfl condtion
		block.dt = Get_CFL(block, fluids, tstop);

		for (int i = 0; i < block.stages; ++i)
		{

			leftboundary(fluids, block, bcvalue[0]);
			rightboundary(fluids, block, bcvalue[1]);

			Reconstruction_within_cell(interfaces, fluids, block);

			Reconstruction_forg0(interfaces, fluids, block);
			
			Calculate_flux(fluxes, interfaces, block, i);

			Update(fluids, fluxes, block, i);
			
			// if (is_using_df_factor)
			// {
			// 	Update_alpha(interfaces, fluids, block);
			// }
		}
		block.step++;
		block.t = block.t + block.dt;
		if ((block.t - tstop) > 0)
		{
			error_for_sinwave(fluids, block, tstop, error);
		}
	}
	runtime.finish_compute = clock();
	cout << "initializiation time is " << (float)(runtime.finish_initial - runtime.start_initial) / CLOCKS_PER_SEC << "seconds" << endl;
	cout << "computational time is " << (float)(runtime.finish_compute - runtime.start_compute) / CLOCKS_PER_SEC << "seconds" << endl;
	output1d(fluids, block);
}

void ICfor_sinwave(Fluid1d* fluids, Block1d block)
{
//#pragma omp parallel for 
	for (int i = block.nodex_begin; i <= block.nodex_end; i++)
	{
		double pi = 3.14159265358979323846;
		fluids[i].primvar[0]
			= 1.0 - 0.2 / pi / block.dx *
			(cos(pi * (i + 1 - block.ghost) * block.dx) - cos(pi * (i - block.ghost) * block.dx));
		fluids[i].primvar[1] = 1;
		fluids[i].primvar[2] = 1;
		Primvar_to_convar_1D(fluids[i].convar, fluids[i].primvar);

	}
}

void error_for_sinwave(Fluid1d* fluids, Block1d block, double tstop, double* error)
{

	if (abs(block.t - tstop) <= (1e-10))
	{

		double error1 = 0;
		double error2 = 0;
		double error3 = 0;

		for (int i = block.ghost; i < block.nx - block.ghost; i++)
		{
			double pi = 3.14159265358979323846;
			int index = i;
			double primvar0 = 1 - 0.2 / pi / block.dx * (cos(pi * (i + 1 - block.ghost) * block.dx) - cos(pi * (i - block.ghost) * block.dx));
			error1 = error1 + abs(fluids[index].convar[0] - primvar0);
			error2 = error2 + pow(fluids[index].convar[0] - primvar0, 2);
			error3 = (error3 > abs(fluids[index].convar[0] - primvar0)) ? error3 : abs(fluids[index].convar[0] - primvar0);
		}
		error1 /= block.nodex;
		error2 = sqrt(error2 / block.nodex);

		error[0] = error1; error[1] = error2; error[2] = error3;

		cout << scientific << "L1 error=" << error1 << " L2 error=" << error2 << " Linf error=" << error3 << endl;

	}
}
