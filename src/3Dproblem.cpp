#include "3Dproblem.h"

void CubicTube()
{
    Runtime runtime;
    runtime.start_initial = clock();

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
}