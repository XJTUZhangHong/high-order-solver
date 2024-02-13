#include "flux_function.h"

// one-dimensional problem
GKS1d_type gks1dsolver = nothing; //initialization
Flux_function flux_function = LF; //initialization

void Calculate_flux(Flux1d** fluxes, Interface1d* interfaces, Block1d& block, int stage)
{
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost + 1; ++i)
	{
		flux_function(fluxes[i][stage], interfaces[i], block.dt);
	}
}

void GKS(Flux1d& flux, Interface1d& interface, double dt)
{
	if (gks1dsolver == nothing)
	{
		cout << "no gks solver specify" << endl;
		exit(0);
	}
	double Flux[2][3];
	//change conservative variables to rho u lambda
	double convar_left[3], convar_right[3], convar0[3];
	for (int i = 0; i < 3; i++)
	{
		convar_left[i] = interface.left.convar[i];
		convar_right[i] = interface.right.convar[i];
		convar0[i] = interface.center.convar[i];
	}

	double prim_left[3], prim_right[3], prim0[3];
	Convar_to_ULambda_1d(prim_left, convar_left);
	Convar_to_ULambda_1d(prim_right, convar_right);
	Convar_to_ULambda_1d(prim0, convar0);
	//then lets get the coefficient of time intergation factors

	double tau;
	double tau_num;
	tau = Get_Tau_NS(prim0[0], prim0[2]);
	tau_num = Get_Tau(prim_left[0], prim_right[0], prim0[0], prim_left[2], prim_right[2], prim0[2], dt);
	double eta = exp(-dt / tau_num);
	double t[10];
	// non equ part time coefficient for gks_2nd algorithm (f0)
	t[0] = tau_num * (1 - eta); // this refers glu, gru part
	t[1] = tau_num * (eta * (dt + tau_num) - tau_num) + tau * tau_num * (eta - 1); //this refers aluu, aruu part
	t[2] = tau * tau_num * (eta - 1); //this refers Alu, Aru part
	// then, equ part time coefficient for gks 2nd (g0)
	t[3] = tau_num * eta + dt - tau_num; //this refers g0u part
	t[4] = tau_num * (tau_num - eta * (dt + tau_num) - tau * (eta - 1)) - dt * tau; //this refers a0uu part
	t[5] = 0.5 * dt * dt - tau * tau_num * (eta - 1) - tau * dt; //this refers A0u part

	if (gks1dsolver == kfvs1st)
	{
		t[0] = dt;
		for (int i = 1; i < 6; i++)
		{
			t[i] = 0.0;
		}
		//do nothing, kfvs1st only use t[0]=dt part;
	}
	else if (gks1dsolver == kfvs2nd)
	{
		t[0] = dt;
		t[1] = -dt * dt / 2.0;
		for (int i = 2; i < 6; i++)
		{
			t[i] = 0.0;
		}
	}
	MMDF1d ml(prim_left[1], prim_left[2]);
	MMDF1d mr(prim_right[1], prim_right[2]);

	double unit[3] = { 1, 0.0, 0.0 };

	double glu[3], gru[3];
	GL(1, 0, glu, unit, ml);
	GR(1, 0, gru, unit, mr);

	//only one part, the kfvs1st part
	for (int i = 0; i < 3; i++)
	{
		Flux[0][i] = prim_left[0] * t[0] * glu[i] + prim_right[0] * t[0] * gru[i];
	}

	if (gks1dsolver == kfvs1st)
	{
		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = Flux[0][i];
		}
		return;
	}
	// kfvs1st part ended

	//now the equ part added, m0 term added, gks1st part begin
	MMDF1d m0(prim0[1], prim0[2]);

	double g0u[3];
	G(1, 0, g0u, unit, m0);

	//the equ g0u part, the gks1st part
	for (int i = 0; i < 3; i++)
	{
		Flux[0][i] = Flux[0][i] + prim0[0] * t[3] * g0u[i];
	}

	if (gks1dsolver == gks1st)
	{
		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = Flux[0][i];
		}
		return;
	}
	// gks1d solver ended

	//for kfvs2nd part
	double der1left[3], der1right[3];
	for (int i = 0; i < 3; i++)
	{
		der1left[i] = interface.left.der1[i];
		der1right[i] = interface.right.der1[i];
	}

	double alx[3];
	Microslope(alx, der1left, prim_left);

	double alxuul[3];
	GL(2, 0, alxuul, alx, ml);

	double arx[3];
	Microslope(arx, der1right, prim_right);
	double arxuur[3];
	GR(2, 0, arxuur, arx, mr);

	for (int i = 0; i < 3; i++)
	{	// t1 part
		Flux[0][i] = Flux[0][i] + prim_left[0] * t[1] * (alxuul[i]) + prim_right[0] * t[1] * (arxuur[i]);
	}
	if (gks1dsolver == kfvs2nd)
	{
		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = Flux[0][i];
		}
		return;
	}
	// the kfvs2nd part ended

	// then we still need t[2], t[4] t[5] part for gks 2nd
	//for t[2] Aru,Alu part
	double alxu[3];
	double arxu[3];

	//take <u> moment for al, ar
	G(1, 0, alxu, alx, ml);
	G(1, 0, arxu, arx, mr);

	double Al[3], Ar[3];
	double der_AL[3], der_AR[3];

	//using compatability condition to get the time derivative
	for (int i = 0; i < 3; i++)
	{
		der_AL[i] = -prim_left[0] * (alxu[i]);
		der_AR[i] = -prim_right[0] * (arxu[i]);
	}
	// solve the coefficient martix b=ma
	Microslope(Al, der_AL, prim_left);
	Microslope(Ar, der_AR, prim_right);

	//to obtain the Alu and Aru
	double Alul[3];
	double Arur[3];
	GL(1, 0, Alul, Al, ml);
	GR(1, 0, Arur, Ar, mr);

	for (int i = 0; i < 3; i++)
	{	// t2 part
		Flux[0][i] = Flux[0][i] + prim_left[0] * t[2] * (Alul[i]) + prim_right[0] * t[2] * (Arur[i]);
	}

	// for t[4] a0xuu part

	double a0x[3];
	double der1[3];

	for (int i = 0; i < 3; i++)
	{
		der1[i] = interface.center.der1[i];
	}

	//solve the microslope
	Microslope(a0x, der1, prim0);
	//a0x <u> moment
	double a0xu[3];
	G(1, 0, a0xu, a0x, m0);
	//a0x <u^2> moment
	double a0xuu[3];
	G(2, 0, a0xuu, a0x, m0);

	for (int i = 0; i < 3; i++)
	{	// t4 part
		Flux[0][i] = Flux[0][i] + prim0[0] * t[4] * (a0xuu[i]);
	}


	// for t[5] A0u part
	double derA0[3];

	for (int i = 0; i < 3; i++)
	{
		derA0[i] = -prim0[0] * (a0xu[i]);
	}
	double A0[3];
	Microslope(A0, derA0, prim0);
	double A0u[3];
	G(1, 0, A0u, A0, m0);
	for (int i = 0; i < 3; i++)
	{	// t5 part
		Flux[0][i] = Flux[0][i] + prim0[0] * t[5] * (A0u[i]);
	}
	if (gks1dsolver == gks2nd && timecoe_list == S1O1)
	{
		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = Flux[0][i];
		}
		return;
	}
	if (gks1dsolver == gks2nd)
	{
		double dt2 = 0.5 * dt;
		tau_num = Get_Tau(prim_left[0], prim_right[0], prim0[0], prim_left[2], prim_right[2], prim0[2], dt);
		eta = exp(-dt2 / tau_num);
		// non equ part time coefficient for gks_2nd algorithm
		t[0] = tau_num * (1 - eta); // this refers glu, gru part
		t[1] = tau_num * (eta * (dt2 + tau_num) - tau_num) + tau * tau_num * (eta - 1); //this refers aluu, aruu part
		t[2] = tau * tau_num * (eta - 1); //this refers Alu, Aru part
		// then, equ part time coefficient for gks 2nd
		t[3] = tau_num * eta + dt2 - tau_num; //this refers g0u part
		t[4] = tau_num * (tau_num - eta * (dt2 + tau_num) - tau * (eta - 1)) - dt2 * tau; //this refers a0uu part
		t[5] = 0.5 * dt2 * dt2 - tau * tau_num * (eta - 1) - tau * dt2; //this refers A0u part

		for (int i = 0; i < 3; i++)
		{
			// t0 part
			Flux[1][i] = prim_left[0] * t[0] * glu[i] + prim_right[0] * t[0] * gru[i];
			// t1 part
			Flux[1][i] = Flux[1][i] + prim_left[0] * t[1] * (alxuul[i]) + prim_right[0] * t[1] * (arxuur[i]);
			// t2 part
			Flux[1][i] = Flux[1][i] + prim_left[0] * t[2] * (Alul[i]) + prim_right[0] * t[2] * (Arur[i]);
			// t3 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[3] * g0u[i];
			// t4 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[4] * (a0xuu[i]);
			// t5 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[5] * (A0u[i]);
		}

		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = (4.0 * Flux[1][i] - Flux[0][i]);
			flux.derf[i] = 4.0 * (Flux[0][i] - 2.0 * Flux[1][i]);
		}

		return;
	}
	else
	{
		cout << "no valid solver specify" << endl;
		exit(0);
	}
}

void LF(Flux1d& flux, Interface1d& interface, double dt)
{
	double pl[3], pr[3];
	Convar_to_primvar_1D(pl, interface.left.convar);
	Convar_to_primvar_1D(pr, interface.right.convar);

	//Reimann invariants, u +/- 2c/(r-1)
	//Sound speed, c = sqrt(r*p/density)
	double k[2];
	k[0] = abs(pl[1]) + sqrt(Gamma * pl[2] / pl[0]); //abs(u)+c, left
	k[1] = abs(pr[1]) + sqrt(Gamma * pr[2] / pr[0]); //abs(u)+c, right

	double beta = k[0]; // beta, means, abs(partialF/partialX)
	if (k[1] > k[0]) { beta = k[1]; }
	double flux_l[3], flux_r[3];
	get_Euler_flux(pl, flux_l);
	get_Euler_flux(pr, flux_r);

	for (int m = 0; m < 3; m++)
	{
		flux.f[m] = 0.5 * ((flux_l[m] + flux_r[m]) - beta * (interface.right.convar[m] - interface.left.convar[m]));
		flux.f[m] *= dt;
	}
}

void get_Euler_flux(double p[3], double* flux)
{
	//the flux of Euler equation, density*u, density*u*u+p, (p+0.5*density*u*u+p/(r-1))*u
	//p[0 1 2], means, density, u, p
	flux[0] = p[0] * p[1];
	flux[1] = p[0] * p[1] * p[1] + p[2];
	double ENERGS = 0.5 * (p[1] * p[1]) * p[0] + p[2] / (Gamma - 1.0);
	flux[2] = p[1] * (ENERGS + p[2]);
}

// two-dimensional problem
GKS2d_type gks2dsolver = kfvs1st_2d; //initilization
Flux_function_2d flux_function_2d = GKS2D; //initilization

void Calculate_flux(Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Interface2d* xinterfaces, Interface2d* yinterfaces, Block2d block, int stage)
{
	// Note : calculate the final flux of gauss points, in xfluxes,  by the reconstructed variables in xinterfaces; the same for y
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost + 1; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				flux_function_2d(xfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss], xinterfaces[i * (block.ny + 1) + j].gauss[num_gauss], block.dt);
				// calculate the final flux in xfluxes,  by the reconstructed variables in xinterfaces; the same for y
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
				flux_function_2d(yfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss], yinterfaces[i * (block.ny + 1) + j].gauss[num_gauss], block.dt);
			}
		}
	}
}

void GKS2D(Flux2d& flux, Recon2d& interface, double dt)
{
	// Note : // calculate the final flux in flux,  by the reconstructed variables in interface
	if (gks2dsolver == nothing_2d)
	{
		cout << "no gks solver specify" << endl;
		exit(0);
	}
	double Flux[2][4];
	double convar_left[4], convar_right[4], convar0[4];
	for (int i = 0; i < 4; ++i)
	{
		convar_left[i] = interface.left.convar[i];
		convar_right[i] = interface.right.convar[i];
	}
	//change conservative variables to rho u lambda
	double prim_left[4], prim_right[4], prim0[4];
	Convar_to_ULambda_2d(prim_left, convar_left);
	Convar_to_ULambda_2d(prim_right, convar_right);

	//prim_left[1], prim_left[2], prim_left[3], means U, V, Lambda
	MMDF ml(prim_left[1], prim_left[2], prim_left[3]);
	MMDF mr(prim_right[1], prim_right[2], prim_right[3]);
	
	if (g0reconstruction_2D_tangent == Center_all_collision_multi)
	{
		Collision(interface.center.convar, prim_left[0], prim_right[0], ml, mr); // ml, mr, pass the whole class into the function
		double a0[4] = { 1.0, 0.0, 0.0, 0.0 };
		double alx_t[4], arx_t[4];
		//w_x
		A(alx_t, interface.left.der1x, prim_left); // input the slope of macroscopic variables, output the a coefficient
		A(arx_t, interface.right.der1x, prim_right);
		double al0x_t[4];
		double ar0x_t[4];
		GL_address(0, 0, 0, al0x_t, alx_t, ml); // ml, mr, pass the whole class into the function
		GR_address(0, 0, 0, ar0x_t, arx_t, mr); // ml, mr, pass the whole class into the function
		//GL_address, GR_address, G_address, are calculating the macroscopic variables by the microscopic variables based on moment calculation
		//Note: when the input of moment calculation is a coefficient, the result is the derivative of W; when the input is (0 0 0 0), the result is W
		//when the input is (1 0 0 0), the result is the relative flux, without considering the integral of time

		//w_y
		double aly_t[4], ary_t[4];
		A(aly_t, interface.left.der1y, prim_left);
		A(ary_t, interface.right.der1y, prim_right);
		//The difference of A_point with A is just the input form, matrix or pointer
		//The content, input variables, output variables are the same

		double al0y_t[4];
		double ar0y_t[4];
		GL_address(0, 0, 0, al0y_t, aly_t, ml);
		GR_address(0, 0, 0, ar0y_t, ary_t, mr);
		for (int var = 0; var < 4; ++var)
		{
			interface.center.der1x[var] = prim_left[0] * al0x_t[var] + prim_right[0] * ar0x_t[var];
			interface.center.der1y[var] = prim_left[0] * al0y_t[var] + prim_right[0] * ar0y_t[var];
		}
	}
	
	for (int i = 0; i < 4; ++i)
	{
		convar0[i] = interface.center.convar[i];
	}

	Convar_to_ULambda_2d(prim0, convar0);
	//then lets get the coefficient of time intergation factors
	double tau;
	double tau_num;
	tau = Get_Tau_NS(prim0[0], prim0[3]);
	tau_num = Get_Tau(prim_left[0], prim_right[0], prim0[0], prim_left[3], prim_right[3], prim0[3], dt);
	double eta = exp(-dt / tau_num);
	double t[10];
	// non equ part time coefficient for gks_2nd algorithm
	t[0] = tau_num * (1 - eta); // this refers glu, gru part
	t[1] = tau_num * (eta * (dt + tau_num) - tau_num) + tau * tau_num * (eta - 1); //this refers aluu, aruu part

	t[2] = tau * tau_num * (eta - 1); //this refers Alu, Aru part
	// then, equ part time coefficient for gks 2nd
	t[3] = tau_num * eta + dt - tau_num; //this refers g0u part
	t[4] = tau_num * (tau_num - eta * (dt + tau_num) - tau * (eta - 1)) - dt * tau; //this refers a0uu part
	t[5] = 0.5 * dt * dt - tau * tau_num * (eta - 1) - tau * dt; //this refers A0u part
	//reset the t[0~5], used ONLY for kfvs method
	//in kfvs, only use time step deltat, but NOT collision time
	if (gks2dsolver == kfvs1st_2d)
	{
		t[0] = dt;
		for (int i = 1; i < 6; i++)
		{
			t[i] = 0.0;
		}
		//do nothing, kfvs1st only use t[0]=dt part;
	}
	else if (gks2dsolver == kfvs2nd_2d)
	{
		t[0] = dt;
		t[1] = -0.5 * dt * dt;
		for (int i = 2; i < 6; i++)
		{
			t[i] = 0.0;
		}
	}

	double unit[4] = { 1, 0.0, 0.0, 0.0 };
	// alx_t = a(0)*1 + a(1)*u + a(2)*v + a(3)*e; arx_t, aly_t, ary_t are the same
	// the slope of the (macroscopic) conservative variables value (/density) can be obtained by (in al0x_t), GL_address(0, 0, 0, al0x_t, alx_t, ml);
	// thus, when replace a by unit (=[1 0 0 0]), GL_address can get the relative (macroscopic) conservative variables (/density); but pay attention to the (1 0 0) but not (0 0 0)
	// G_address, GR_address, are the same with GL_address

	double glu[4], gru[4];

	GL_address(1, 0, 0, glu, unit, ml); // (1 0 0) get GLu
	GR_address(1, 0, 0, gru, unit, mr); // (1 0 0) get GRu
	//only one part, the kfvs1st part
	for (int i = 0; i < 4; i++)
	{
		Flux[0][i] = prim_left[0] * t[0] * glu[i] + prim_right[0] * t[0] * gru[i];
	}
	//cout<< endl;
	if (gks2dsolver == kfvs1st_2d)
	{
		for (int i = 0; i < 4; i++)
		{
			flux.f[i] = Flux[0][i];
		}

		return;
	}
	// kfvs1st part ended

	//now the equ part added, m0 term added, gks1st part begin
	MMDF m0(prim0[1], prim0[2], prim0[3]);

	double g0u[4];
	G_address(1, 0, 0, g0u, unit, m0);

	if (gks2dsolver != kfvs2nd_2d)
	{
		//the equ g0u part, the gks1st part
		for (int i = 0; i < 4; i++)
		{
			Flux[0][i] = Flux[0][i] + prim0[0] * t[3] * g0u[i];
		}
	}

	if (gks2dsolver == gks1st_2d)
	{
		for (int i = 0; i < 4; i++)
		{
			flux.f[i] = Flux[0][i];
		}
		return;
	}
	// gks1d solver ended

	//for kfvs2nd part
	double der1xleft[4], der1xright[4];
	double der1yleft[4], der1yright[4];
	for (int i = 0; i < 4; i++)
	{
		der1xleft[i] = interface.left.der1x[i];
		der1xright[i] = interface.right.der1x[i];
		der1yleft[i] = interface.left.der1y[i];
		der1yright[i] = interface.right.der1y[i];
	}

	double alx[4];
	A(alx, der1xleft, prim_left);
	double alxuul[4];
	GL_address(2, 0, 0, alxuul, alx, ml);

	double arx[4];
	A(arx, der1xright, prim_right);
	double arxuur[4];
	GR_address(2, 0, 0, arxuur, arx, mr);

	double aly[4];
	A(aly, der1yleft, prim_left);
	double alyuvl[4];
	GL_address(1, 1, 0, alyuvl, aly, ml);

	double ary[4];
	A(ary, der1yright, prim_right);
	double aryuvr[4];
	GR_address(1, 1, 0, aryuvr, ary, mr);

	for (int i = 0; i < 4; i++)
	{	// t1 part
		Flux[0][i] = Flux[0][i] + t[1] * (prim_left[0] * (alxuul[i] + alyuvl[i]) + prim_right[0] * (arxuur[i] + aryuvr[i]));
	}
	if (gks2dsolver == kfvs2nd_2d)
	{
		for (int i = 0; i < 4; i++)
		{
			flux.f[i] = Flux[0][i];
		}
		return;
	}
	// the kfvs2nd part ended

	// then we still need t[2], t[4] t[5] part for gks 2nd
	//for t[2] Aru,Alu part
	double alxu[4];		double alyv[4];
	double arxu[4];		double aryv[4];

	//take <u> moment for al, ar
	G_address(1, 0, 0, alxu, alx, ml);
	G_address(1, 0, 0, arxu, arx, mr);
	G_address(0, 1, 0, alyv, aly, ml);
	G_address(0, 1, 0, aryv, ary, mr);

	double Al[4], Ar[4];
	double der_AL[4], der_AR[4];

	//using compatability condition to get the time derivative
	for (int i = 0; i < 4; i++)
	{
		der_AL[i] = -prim_left[0] * (alxu[i] + alyv[i]);
		der_AR[i] = -prim_right[0] * (arxu[i] + aryv[i]);
	}
	// solve the coefficient martix b=ma
	A(Al, der_AL, prim_left);
	A(Ar, der_AR, prim_right);

	//to obtain the Alu and Aru
	double Alul[4];
	double Arur[4];
	GL_address(1, 0, 0, Alul, Al, ml);
	GR_address(1, 0, 0, Arur, Ar, mr);

	for (int i = 0; i < 4; i++)
	{	// t2 part
		Flux[0][i] = Flux[0][i] + t[2] * (prim_left[0] * Alul[i] + prim_right[0] * Arur[i]);
		
	}


	// for t[4] a0xuu part
	double a0x[4];
	double der1[4];
	for (int i = 0; i < 4; i++)
	{
		der1[i] = interface.center.der1x[i]; //Only the averaged value for g0
	}
	//solve the microslope
	A(a0x, der1, prim0);
	
	//a0x <u> moment
	double a0xu[4];
	G_address(1, 0, 0, a0xu, a0x, m0); //get a0xu, used for the following determination of derA0, and then A0
	//a0x <u^2> moment
	double a0xuu[4];
	G_address(2, 0, 0, a0xuu, a0x, m0);

	double a0y[4];
	double dery[4];
	for (int i = 0; i < 4; i++)
	{
		dery[i] = interface.center.der1y[i];
	}
	A(a0y, dery, prim0);
	double a0yv[4];
	G_address(0, 1, 0, a0yv, a0y, m0); //get a0yv, used for the following determination of derA0, and then A0
	//a0x <u^2> moment
	double a0yuv[4];
	G_address(1, 1, 0, a0yuv, a0y, m0);

	for (int i = 0; i < 4; i++)
	{	// t4 part
		Flux[0][i] = Flux[0][i] + prim0[0] * t[4] * (a0xuu[i] + a0yuv[i]);
		
	}

	// for t[5] A0u part
	double derA0[4];
	for (int i = 0; i < 4; i++)
	{
		derA0[i] = -prim0[0] * (a0xu[i] + a0yv[i]);
	}
	double A0[4];
	A(A0, derA0, prim0);
	double A0u[4];
	G_address(1, 0, 0, A0u, A0, m0);
	for (int i = 0; i < 4; i++)
	{	// t5 part
		Flux[0][i] = Flux[0][i] + prim0[0] * t[5] * (A0u[i]);
	}
	if (is_Prandtl_fix == true)
	{
		double qs;
		qs = a0xuu[3] + A0u[3] - prim0[1] * (a0xuu[1] + A0u[1]) - prim0[2] * (a0xuu[2] + A0u[2]);
		qs *= -prim0[0] * tau * dt;
		Flux[0][3] += (1.0 / Pr - 1.0) * qs;
	}
	if (gks2dsolver == gks2nd_2d && timecoe_list_2d == S1O1_2D)
	{
		for (int i = 0; i < 4; i++)
		{
			flux.f[i] = Flux[0][i];
			flux.derf[i] = 0.0;
		}
		return;
	}

	if (gks2dsolver == gks2nd_2d && timecoe_list_2d != S1O1_2D)
	{
		double dt2 = 0.5 * dt; // the following is dt2
		tau_num = Get_Tau(prim_left[0], prim_right[0], prim0[0], prim_left[3], prim_right[3], prim0[3], dt2);
		
		eta = exp(-dt2 / tau_num);
		// non equ part time coefficient for gks_2nd algorithm
		t[0] = tau_num * (1 - eta); // this refers glu, gru part
		t[1] = tau_num * (eta * (dt2 + tau_num) - tau_num) + tau * tau_num * (eta - 1); //this refers aluu, aruu part
		t[2] = tau * tau_num * (eta - 1); //this refers Alu, Aru part
		// then, equ part time coefficient for gks 2nd
		t[3] = tau_num * eta + dt2 - tau_num; //this refers g0u part
		t[4] = tau_num * (tau_num - eta * (dt2 + tau_num) - tau * (eta - 1)) - dt2 * tau; //this refers a0uu part
		t[5] = 0.5 * dt2 * dt2 - tau * tau_num * (eta - 1) - tau * dt2; //this refers A0u part

		for (int i = 0; i < 4; i++)
		{
			// t0 part
			Flux[1][i] = t[0] * (prim_left[0] * glu[i] + prim_right[0] * gru[i]);
			// t1 part
			Flux[1][i] = Flux[1][i] + t[1] * (prim_left[0] * (alxuul[i] + alyuvl[i]) + prim_right[0] * (arxuur[i] + aryuvr[i]));
			// t2 part
			Flux[1][i] = Flux[1][i] + t[2] * (prim_left[0] * Alul[i] + prim_right[0] * Arur[i]);
			// t3 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[3] * g0u[i];
			// t4 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[4] * (a0xuu[i] + a0yuv[i]);
			// t5 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[5] * (A0u[i]);
		}

		if (is_Prandtl_fix == true)
		{
			double qs;
			qs = a0xuu[3] + A0u[3] - prim0[1] * (a0xuu[1] + A0u[1]) - prim0[2] * (a0xuu[2] + A0u[2]);
			qs *= -prim0[0] * tau * 0.5 * dt;
			Flux[1][3] += (1.0 / Pr - 1.0) * qs;
		}

		for (int i = 0; i < 4; i++)
		{
			flux.f[i] = (4.0 * Flux[1][i] - Flux[0][i]);
			flux.derf[i] = 4.0 * (Flux[0][i] - 2.0 * Flux[1][i]);
		}
		return;
	}
	else
	{
		cout << "no valid solver specify" << endl;
		exit(0);
	}
}