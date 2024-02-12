#include "reconstruction.h"

Reconstruction_within_Cell cellreconstruction = Vanleer; //initialization
Reconstruction_forG0 g0reconstruction = Center_collision; //initialization
Reconstruction_variable reconstruction_variable = conservative; //initialization
WENOtype wenotype = wenojs; //initialization
bool is_reduce_order_warning = false; //initialization

void Check_Order_Reduce(Point1d& left, Point1d& right, Fluid1d& fluid)
{
	bool order_reduce[2];
	Check_Order_Reduce_by_Lambda_1D(order_reduce[0], right.convar);
	Check_Order_Reduce_by_Lambda_1D(order_reduce[1], left.convar);
	//if lambda <0, then reduce to the first order
	if (order_reduce[0] == true || order_reduce[1] == true)
	{
		//cout << "warning" << endl;
		//if (is_reduce_order_warning == true)
			//cout << "reconstruction reduce order in location  " << fluid.cx << endl;
		for (int m = 0; m < 3; m++)
		{
			left.convar[m] = fluid.convar[m];
			right.convar[m] = fluid.convar[m];
			left.der1[m] = 0.0;
			right.der1[m] = 0.0;
		}
	}
}

void Check_Order_Reduce_by_Lambda_1D(bool& order_reduce, double* convar)
{
	order_reduce = false;
	double lambda;
	lambda = Lambda(convar[0], convar[1] / convar[0], convar[2]);
	//if lambda <0, then reduce to the first order
	if (lambda <= 0.0 || (lambda == lambda) == false)
	{
		order_reduce = true;
	}
}

void Reconstruction_within_cell(Interface1d* interfaces, Fluid1d* fluids, Block1d block)
{
#pragma omp parallel  for
	for (int i = block.ghost - 1; i < block.nx - block.ghost + 1; i++)
	{
		cellreconstruction(interfaces[i].right, interfaces[i + 1].left, &fluids[i], block);
	}
}

void Vanleer(Point1d& left, Point1d& right, Fluid1d* fluids, Block1d block)
{
	//Note: function by vanleer reconstruction
	//the final updated variables is conservative variables
	Fluid1d wn1 = fluids[-1];
	Fluid1d w = fluids[0];
	Fluid1d wp1 = fluids[1];
	double splus[3], sminus[3];

	for (int i = 0; i < 3; i++)
	{
		splus[i] = (wp1.convar[i] - w.convar[i]) / ((wp1.dx + w.dx) / 2.0);
		sminus[i] = (w.convar[i] - wn1.convar[i]) / ((wn1.dx + w.dx) / 2.0);

		if ((splus[i] * sminus[i]) > 0)
		{
			left.der1[i] = 2 * splus[i] * sminus[i] / (splus[i] + sminus[i]);
			right.der1[i] = left.der1[i];
		}
		else
		{
			left.der1[i] = 0.0;
			right.der1[i] = 0.0;
		}
		left.convar[i] = w.convar[i] - 0.5 * w.dx * left.der1[i];
		right.convar[i] = w.convar[i] + 0.5 * w.dx * right.der1[i];
	}
	Check_Order_Reduce(left, right, fluids[0]);
}

void WENO5_AO(Point1d& left, Point1d& right, Fluid1d* fluids, Block1d block)
{
	//Note: function by WENO5_AO reconstruction

	double wn2[3]; Copy_Array(wn2, fluids[-2].convar, 3);
	double wn1[3]; Copy_Array(wn1, fluids[-1].convar, 3);
	double w0[3];  Copy_Array(w0, fluids[0].convar, 3);
	double wp1[3]; Copy_Array(wp1, fluids[1].convar, 3);
	double wp2[3]; Copy_Array(wp2, fluids[2].convar, 3);
	double tmp;
	//non-uniform grid was treated as uniform grid
	double h = fluids[0].dx;
	double beta[4];

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 3; i++)
		{
			beta[0] = 13.0 / 12.0 * pow((wn2[i] - 2.0 * wn1[i] + w0[i]), 2) + 0.25 * pow((wn2[i] - 4.0 * wn1[i] + 3.0 * w0[i]), 2);
			beta[1] = 13.0 / 12.0 * pow((wn1[i] - 2.0 * w0[i] + wp1[i]), 2) + 0.25 * pow((wn1[i] - wp1[i]), 2);
			beta[2] = 13.0 / 12.0 * pow((w0[i] - 2.0 * wp1[i] + wp2[i]), 2) + 0.25 * pow((3.0 * w0[i] - 4.0 * wp1[i] + wp2[i]), 2);
			beta[3] = (1.0 / 5040.0) * (231153.0 * w0[i] * w0[i] + 104963.0 * wn1[i] * wn1[i] + 6908.0 * wn2[i] * wn2[i] -
			38947.0 * wn2[i] * wp1[i] + 104963.0 * wp1[i] * wp1[i] +
			wn1[i] * (-51001.0 * wn2[i] + 179098.0 * wp1[i] - 38947.0 * wp2[i]) -
			3.0 * w0[i] * (99692.0 * wn1[i] - 22641.0 * wn2[i] + 99692.0 * wp1[i] - 22641.0 * wp2[i]) +
			8209.0 * wn2[i] * wp2[i] - 51001.0 * wp1[i] * wp2[i] + 6908.0 * wp2[i] * wp2[i]);
			
			weno_5th_ao_right(right.convar[i], right.der1[i], tmp, wn2[i], wn1[i], w0[i], wp1[i], wp2[i], beta, h);
			weno_5th_ao_left(left.convar[i], left.der1[i], tmp, wn2[i], wn1[i], w0[i], wp1[i], wp2[i], beta, h);
		}
	}

	if (reconstruction_variable == characteristic)
	{

		double l_ren2[3], l_ren1[3], l_re0[3], l_rep1[3], l_rep2[3];
		double r_ren2[3], r_ren1[3], r_re0[3], r_rep1[3], r_rep2[3];
		double l_var[3], l_der1[3], l_der2[3];
		double r_var[3], r_der1[3], r_der2[3];

		double base_left[3], base_right[3];

		double wn1_primvar[3], w_primvar[3], wp1_primvar[3];
		Convar_to_primvar_1D(wn1_primvar, wn1);
		Convar_to_primvar_1D(w_primvar, w0);
		Convar_to_primvar_1D(wp1_primvar, wp1);

		for (int i = 0; i < 3; i++)
		{
			base_left[i] = 0.5 * (wn1_primvar[i] + w_primvar[i]);
			base_right[i] = 0.5 * (wp1_primvar[i] + w_primvar[i]);
		}
		// left side
		Convar_to_char1D(l_ren2, base_left, wn2);
		Convar_to_char1D(l_ren1, base_left, wn1);
		Convar_to_char1D(l_re0, base_left, w0);
		Convar_to_char1D(l_rep1, base_left, wp1);
		Convar_to_char1D(l_rep2, base_left, wp2);
		// right side
		Convar_to_char1D(r_ren2, base_right, wn2);
		Convar_to_char1D(r_ren1, base_right, wn1);
		Convar_to_char1D(r_re0, base_right, w0);
		Convar_to_char1D(r_rep1, base_right, wp1);
		Convar_to_char1D(r_rep2, base_right, wp2);

		for (int i = 0; i < 3; i++)
		{
			beta[0] = 13.0 / 12.0 * pow((wn2[i] - 2.0 * wn1[i] + w0[i]), 2) + 0.25 * pow((wn2[i] - 4.0 * wn1[i] + 3.0 * w0[i]), 2);
			beta[1] = 13.0 / 12.0 * pow((wn1[i] - 2.0 * w0[i] + wp1[i]), 2) + 0.25 * pow((wn1[i] - wp1[i]), 2);
			beta[2] = 13.0 / 12.0 * pow((w0[i] - 2.0 * wp1[i] + wp2[i]), 2) + 0.25 * pow((3.0 * w0[i] - 4.0 * wp1[i] + wp2[i]), 2);
			beta[3] = (1.0 / 5040.0) * (231153.0 * w0[i] * w0[i] + 104963.0 * wn1[i] * wn1[i] + 6908.0 * wn2[i] * wn2[i] -
			38947.0 * wn2[i] * wp1[i] + 104963.0 * wp1[i] * wp1[i] +
			wn1[i] * (-51001.0 * wn2[i] + 179098.0 * wp1[i] - 38947.0 * wp2[i]) -
			3.0 * w0[i] * (99692.0 * wn1[i] - 22641.0 * wn2[i] + 99692.0 * wp1[i] - 22641.0 * wp2[i]) +
			8209.0 * wn2[i] * wp2[i] - 51001.0 * wp1[i] * wp2[i] + 6908.0 * wp2[i] * wp2[i]);
			
			weno_5th_ao_left(l_var[i], l_der1[i], l_der2[i], l_ren2[i], l_ren1[i], l_re0[i], l_rep1[i], l_rep2[i], beta, h);
			weno_5th_ao_right(r_var[i], r_der1[i], r_der2[i], r_ren2[i], r_ren1[i], r_re0[i], r_rep1[i], r_rep2[i], beta, h);
		}
		Char_to_convar1D(left.convar, base_left, l_var);
		Char_to_convar1D(left.der1, base_left, l_der1);
		Char_to_convar1D(right.convar, base_right, r_var);
		Char_to_convar1D(right.der1, base_right, r_der1);
	}
	Check_Order_Reduce(left, right, fluids[0]);
}

void weno_5th_ao_left(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double* beta, double h)
{
	double dhi = 0.85;
	double dlo = 0.85;
	//-- - parameter of WENO-- -
	double d[4], ww[4], alpha[4];
	double epsilonW = 1e-10;
	//-- - intermediate parameter-- -
	double p[4], px[4], pxx[4], tempvar;
	double sum_alpha;

	//three small stencil
	d[0] = (1.0 - dhi) * (1.0 - dlo) / 2.0;
	d[1] = (1.0 - dhi) * dlo;
	d[2] = (1.0 - dhi) * (1.0 - dlo) / 2.0;
	//one big stencil
	d[3] = dhi;

	if (wenotype == linear)
	{
		for (int k = 0; k < 4; k++)
		{
			ww[k] = d[k];
		}
	}
	else
	{
		double tau5 = 1.0 / 3.0 * (abs(beta[3] - beta[0]) + abs(beta[3] - beta[1]) + abs(beta[3] - beta[2]));

		if (wenotype == wenojs)
		{
			sum_alpha = 0.0;
			for (int k = 0; k < 4; k++)
			{
				alpha[k] = d[k] / ((epsilonW + beta[k]) * (epsilonW + beta[k]));
				sum_alpha += alpha[k];
			}
		}
		else if (wenotype == wenoz)
		{
			sum_alpha = 0.0;
			for (int i = 0; i < 4; i++)
			{
				double global_div = tau5 / (beta[i] + epsilonW);
				alpha[i] = d[i] * (1 + global_div * global_div);
				sum_alpha += alpha[i];
			}
		}

		for (int k = 0; k < 4; k++)
		{
			ww[k] = alpha[k] / sum_alpha;
		}
	}
	//-- - candidate polynomial-- -
	p[0] = -1.0 / 6.0 * wn2 + 5.0 / 6.0 * wn1 + 1.0 / 3.0 * w0;
	p[1] = 1.0 / 3.0 * wn1 + 5.0 / 6.0 * w0 - 1.0 / 6.0 * wp1;
	p[2] = 11.0 / 6.0 * w0 - 7.0 / 6.0 * wp1 + 1.0 / 3.0 * wp2;
	p[3] = (1.0 / 60.0) * (47.0 * w0 + 27.0 * wn1 - 3.0 * wn2 - 13.0 * wp1 + 2.0 * wp2);

	px[0] = (w0 - wn1) / h;
	px[1] = (w0 - wn1) / h;
	px[2] = -((2.0 * w0 - 3.0 * wp1 + wp2) / h);
	px[3] = (15.0 * w0 - 15.0 * wn1 + wn2 - wp1) / (12.0 * h);

	pxx[0] = (w0 - 2.0 * wn1 + wn2) / h / h;
	pxx[1] = (-2.0 * w0 + wn1 + wp1) / h / h;
	pxx[2] = (w0 - 2.0 * wp1 + wp2) / h / h;
	pxx[3] = ((-8.0 * w0 + 2.0 * wn1 + wn2 + 6.0 * wp1 - wp2) / (4.0 * h * h));

	//-- - combination-- -
	var = 0.0;
	der1 = 0.0;
	der2 = 0.0;
	double final_weight[4];
	final_weight[3] = ww[3] / d[3];
	for (int k = 0; k < 3; k++)
	{
		final_weight[k] = ww[k] - ww[3] / d[3] * d[k];
	}

	for (int k = 0; k < 4; k++)
	{
		var += final_weight[k] * p[k];
		der1 += final_weight[k] * px[k];
		der2 += final_weight[k] * pxx[k];
	}
}

void weno_5th_ao_right(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double* beta, double h)
{
	double dhi = 0.85;
	double dlo = 0.85;
	//-- - parameter of WENO-- -
	double d[4], ww[4], alpha[4];
	double epsilonW = 1e-10;

	//-- - intermediate parameter-- -
	double p[4], px[4], pxx[4], tempvar;
	double sum_alpha;

	//three small stencil
	d[0] = (1 - dhi) * (1 - dlo) / 2.0;
	d[1] = (1 - dhi) * dlo;
	d[2] = (1 - dhi) * (1 - dlo) / 2.0;
	//one big stencil
	d[3] = dhi;

	if (wenotype == linear)
	{
		for (int k = 0; k < 4; k++)
		{
			ww[k] = d[k];
		}
	}
	else
	{
		double tau5 = 1.0 / 3.0 * (abs(beta[3] - beta[0]) + abs(beta[3] - beta[1]) + abs(beta[3] - beta[2]));

		if (wenotype == wenojs)
		{
			sum_alpha = 0.0;
			for (int k = 0; k < 4; k++)
			{
				alpha[k] = d[k] / ((epsilonW + beta[k]) * (epsilonW + beta[k]));
				sum_alpha += alpha[k];
			}
		}
		else if (wenotype == wenoz)
		{
			sum_alpha = 0.0;
			for (int i = 0; i < 4; i++)
			{
				double global_div = tau5 / (beta[i] + epsilonW);
				alpha[i] = d[i] * (1.0 + global_div * global_div);
				sum_alpha += alpha[i];
			}
		}

		for (int k = 0; k < 4; k++)
		{
			ww[k] = alpha[k] / sum_alpha;
		}

	}
	//-- - candidate polynomial-- -

	p[0] = 1.0 / 3.0 * wn2 - 7.0 / 6.0 * wn1 + 11.0 / 6.0 * w0;
	p[1] = -1.0 / 6.0 * wn1 + 5.0 / 6.0 * w0 + 1.0 / 3.0 * wp1;
	p[2] = 1.0 / 3.0 * w0 + 5.0 / 6.0 * wp1 - 1.0 / 6.0 * wp2;
	p[3] = (1.0 / 60.0) * (47.0 * w0 - 13.0 * wn1 + 2.0 * wn2 + 27.0 * wp1 - 3.0 * wp2);

	px[0] = (2.0 * w0 - 3.0 * wn1 + wn2) / h;
	px[1] = (-w0 + wp1) / h;
	px[2] = (-w0 + wp1) / h;
	px[3] = (-15.0 * w0 + wn1 + 15.0 * wp1 - wp2) / (12.0 * h);

	pxx[0] = (w0 - 2.0 * wn1 + wn2) / h / h;
	pxx[1] = (-2.0 * w0 + wn1 + wp1) / h / h;
	pxx[2] = (w0 - 2.0 * wp1 + wp2) / h / h;
	pxx[3] = (-8.0 * w0 + 6.0 * wn1 - wn2 + 2.0 * wp1 + wp2) / (4.0 * h * h);

	//-- - combination-- -
	var = 0.0;
	der1 = 0.0;
	der2 = 0.0;
	double final_weight[4];
	final_weight[3] = ww[3] / d[3];
	for (int k = 0; k < 3; k++)
	{
		final_weight[k] = ww[k] - ww[3] / d[3] * d[k];
	}

	for (int k = 0; k < 4; k++)
	{
		var += final_weight[k] * p[k];
		der1 += final_weight[k] * px[k];
		der2 += final_weight[k] * pxx[k];
	}
}

void Reconstruction_forg0(Interface1d* interfaces, Fluid1d* fluids, Block1d block)
{
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nx - block.ghost + 1; ++i)
	{
		g0reconstruction(interfaces[i], &fluids[i - 1]);
	}
}

void Center_collision(Interface1d& interfaces, Fluid1d* fluids)
{
	double convar_left[3], convar_right[3];
	for (int i = 0; i < 3; i++)
	{
		convar_left[i] = interfaces.left.convar[i];
		convar_right[i] = interfaces.right.convar[i];
	}

	double prim_left[3], prim_right[3]; //rho, U, lambda
	Convar_to_ULambda_1d(prim_left, convar_left);
	Convar_to_ULambda_1d(prim_right, convar_right);

	MMDF1d ml(prim_left[1], prim_left[2]);
	MMDF1d mr(prim_right[1], prim_right[2]);
	
	double unit[3]{ 1.0, 0.0, 0.0 };

	double gl[3], gr[3];
	GL(0, 0, gl, unit, ml); // gl, means Wl(u>0), by input uint
	GR(0, 0, gr, unit, mr); // gr, means Wr(u<0), by input uint
	double axl[3], axr[3];
	Microslope(axl, interfaces.left.der1, prim_left); // axl, means a coefficient indicating slope
	Microslope(axr, interfaces.right.der1, prim_right); // axr, means a coefficient indicating slope
	double ax0l[3], ax0r[3];
	GL(0, 0, ax0l, axl, ml);  // ax0l, means Wlx(u>0), by input axl
	GR(0, 0, ax0r, axr, mr); // ax0r, means Wrx(u<0), by input axr
	for (int i = 0; i < 3; i++)
	{
		interfaces.center.convar[i] = convar_left[0] * gl[i] + convar_right[0] * gr[i];
		interfaces.center.der1[i] = convar_left[0] * ax0l[i] + convar_right[0] * ax0r[i];
	}
}