#include "3Dproblem.h"

void CubicTube()
{
    Runtime runtime;
    runtime.start_initial = omp_get_wtime();

    Block3d block;
    block.uniform = true;
	block.nodex = 50;
	block.nodey = 50;
	block.nodez = 50;
	block.ghost = 3;

	double tstop = 0.2;
	block.CFL = 0.5;
	Fluid3d *bcvalue = new Fluid3d[6];

    K = 2;
	Gamma = 1.4;    

    //this part should rewritten ad gks2dsolver blabla
	gks3dsolver = gks2nd_3d;
	tau_type = Euler;
	c1_euler = 0.05;
	c2_euler = 1;
	//prepare the boundary condtion function
	BoundaryCondition3d leftboundary(0);
	BoundaryCondition3d rightboundary(0);
	BoundaryCondition3d downboundary(0);
	BoundaryCondition3d upboundary(0);
	BoundaryCondition3d backboundary(0);
	BoundaryCondition3d frontboundary(0);

    leftboundary = free_boundary_xleft;
	rightboundary = free_boundary_xright;
	downboundary = free_boundary_yleft;
	upboundary = free_boundary_yright;
	backboundary = free_boundary_zleft;
	frontboundary = free_boundary_zright;

	//prepare the reconstruction function
	// for 3-D reconstruction there is no need to set gauss point 
    gausspoint = 4;
    SetGuassPoint_3D();

    reconstruction_variable = characteristic;
	wenotype = wenoz;

    cellreconstruction_3D_of_face_value = WENO5_AO_splitting_3d;
	cellreconstruction_3D_of_line_value = WENO5_AO_of_line_value;
	cellreconstruction_3D_of_point_value = WENO5_AO_of_point_value;
	g0reconstruction_3D_of_face_value = Do_nothing_splitting_3d;
	g0reconstruction_3D_of_line_value = Do_nothing_reconstruction_of_line_value;
	g0reconstruction_3D_of_point_value = Do_nothing_reconstruction_of_point_value;

	gks3dsolver = gks2nd_3d;
	g0type = all_collisionn;
	flux_function_3d = GKS3D;

	timecoe_list_3d = S2O4_3D;
	Initial_stages(block);

	// allocate memory for 2-D fluid field
	// in a standard finite element method, we have 
	// first the cell average value, N*N
	block.nx = block.nodex + 2 * block.ghost;
	block.ny = block.nodey + 2 * block.ghost;
	block.nz = block.nodez + 2 * block.ghost;
	Fluid3d *fluids = Setfluid_array(block);

	// then the interfaces reconstructioned value, (N+1)*(N+1)
	Interface3d *xinterfaces = Setinterface_array(block);
	Interface3d *yinterfaces = Setinterface_array(block);
	Interface3d *zinterfaces = Setinterface_array(block);
	// then the flux, which the number is identical to interfaces
	Flux3d_gauss **xfluxes = Setflux_gauss_array(block);
	Flux3d_gauss **yfluxes = Setflux_gauss_array(block);
	Flux3d_gauss **zfluxes = Setflux_gauss_array(block);
	//end the allocate memory part

	//bulid or read mesh,
	//since the mesh is all structured from left to right, down to up
	//there is no need for buliding the complex topology between cells and interfaces
	//just using the index for address searching
	block.left = 0.0; block.right = 1.0;
	block.down = 0.0; block.up = 1.0;
	block.back = 0.0; block.front = 1.0;
	block.dx = (block.right - block.left) / block.nodex;
	block.dy = (block.up - block.down) / block.nodey;
	block.dz = (block.front - block.back) / block.nodez;
	block.xarea = block.dy*block.dz; block.yarea = block.dx*block.dz; block.zarea = block.dx*block.dy;
	block.overdx = 1 / block.dx; block.overdy = 1 / block.dy; block.overdz = 1 / block.dz;
	//set the uniform geometry information
	SetUniformMesh(block, fluids, xinterfaces, yinterfaces, zinterfaces);
	//ended mesh part

	double zone1[5] = { 1.0, 0,0, 0, 1.0 };
	double zone2[5] = { 0.125, 0,0, 0, 0.1 };

	ICfor_cubic_sod(fluids, zone1, zone2, block);
	runtime.finish_initial = omp_get_wtime();
	//then we are about to do the loop
	block.t = 0;//the current simulation time
	block.step = 0; //the current step

	int inputstep = 1;//input a certain step,
					  //initialize inputstep=1, to avoid a 0 result
	
}

void ICfor_cubic_sod(Fluid3d *fluids, double * zone1, double * zone2, Block3d block)
{

	for (int i = block.ghost; i < block.nx - block.ghost; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost; j++)
		{
			for (int k = block.ghost; k < block.nz - block.ghost; k++)
			{
				int index = i*(block.ny*block.nz) + j*block.nz + k;

				double xr = std::abs(fluids[index].loc[0] - 0.5*(block.right + block.left));
				double yr = std::abs(fluids[index].loc[1] - 0.5*(block.up + block.down));
				double zr = std::abs(fluids[index].loc[2] - 0.5*(block.front + block.back));
				double r = sqrt((xr*xr) + (yr*yr) + (zr*zr));
				if (xr<0.2*(block.right - block.left)&&yr<0.2*(block.up - block.down) &&zr<0.2*(block.front - block.back))
				{
					for (int var = 0; var < 5; var++)
					{

						fluids[index].primvar[var] = zone1[var];
					}
				}
				else
				{
					for (int var = 0; var < 5; var++)
					{

						fluids[index].primvar[var] = zone2[var];
					}
				}
			}
		}
	}
#pragma omp parallel for 
	for (int i = block.ghost; i < block.nx - block.ghost; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost; j++)
		{
			for (int k = block.ghost; k < block.nz - block.ghost; k++)
			{
				int index = i*(block.ny*block.nz) + j*block.nz + k;
				Primvar_to_convar_3D(fluids[index].convar, fluids[index].primvar);
			}
		}
	}
}