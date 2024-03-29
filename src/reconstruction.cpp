#include "reconstruction.h"

// one-dimensional problem
Reconstruction_within_Cell cellreconstruction = Vanleer; //initialization
Reconstruction_forG0 g0reconstruction = Center_collision; //initialization
Reconstruction_variable reconstruction_variable = conservative; //initialization
WENOtype wenotype = wenojs; //initialization
double df_thres = 3.0;
bool is_reduce_order_warning = false; //initialization
bool is_using_df_factor = false;
bool smooth = false;

void Check_Order_Reduce(Point1d& left, Point1d& right, Fluid1d& fluid)
{
	bool order_reduce[2];
	Check_Order_Reduce_by_Lambda_1D(order_reduce[0], right.convar);
	Check_Order_Reduce_by_Lambda_1D(order_reduce[1], left.convar);
	//if lambda <0, then reduce to the first order
	if (order_reduce[0] == true || order_reduce[1] == true)
	{
		//cout << "warning" << endl;
		if (is_reduce_order_warning == true)
			cout << "reconstruction reduce order in location  " << fluid.cx << endl;
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
	double  h = fluids[0].dx;

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 3; i++)
		{
			weno_5th_ao_right(right.convar[i], right.der1[i], tmp, wn2[i], wn1[i], w0[i], wp1[i], wp2[i], h);
			weno_5th_ao_left(left.convar[i], left.der1[i], tmp, wn2[i], wn1[i], w0[i], wp1[i], wp2[i], h);
		}
	}

	if (reconstruction_variable == characteristic)
	{

		double ren3[3], ren2[3], ren1[3], re0[3], rep1[3], rep2[3], rep3[3];
		double var[3], der1[3], der2[3];

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

		Convar_to_char1D(ren2, base_left, wn2);
		Convar_to_char1D(ren1, base_left, wn1);
		Convar_to_char1D(re0, base_left, w0);
		Convar_to_char1D(rep1, base_left, wp1);
		Convar_to_char1D(rep2, base_left, wp2);

		// left_reconstruction
		for (int i = 0; i < 3; i++)
		{
			weno_5th_ao_left(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h);
		}
		Char_to_convar1D(left.convar, base_left, var);
		Char_to_convar1D(left.der1, base_left, der1);
	
		// right reconstruction

		Convar_to_char1D(ren2, base_right, wn2);
		Convar_to_char1D(ren1, base_right, wn1);
		Convar_to_char1D(re0, base_right, w0);
		Convar_to_char1D(rep1, base_right, wp1);
		Convar_to_char1D(rep2, base_right, wp2);

		for (int i = 0; i < 3; i++)
		{
			weno_5th_ao_right(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h);
		}
		Char_to_convar1D(right.convar, base_right, var);
		Char_to_convar1D(right.der1, base_right, der1);

	}

	Check_Order_Reduce(left, right, fluids[0]);

}

void weno_5th_ao_left(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double h)
{
	double dhi = 0.85;
	double dlo = 0.85;
	//-- - parameter of WENO-- -
	double beta[4], d[4], ww[4], alpha[4];
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
		//cout << "here" << endl;
		beta[0] = 13.0 / 12.0 * pow((wn2 - 2.0 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4.0 * wn1 + 3.0 * w0), 2);
		beta[1] = 13.0 / 12.0 * pow((wn1 - 2.0 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
		beta[2] = 13.0 / 12.0 * pow((w0 - 2.0 * wp1 + wp2), 2) + 0.25 * pow((3.0 * w0 - 4.0 * wp1 + wp2), 2);
		beta[3] = (1.0 / 5040.0) * (231153.0 * w0 * w0 + 104963.0 * wn1 * wn1 + 6908.0 * wn2 * wn2 -
			38947.0 * wn2 * wp1 + 104963.0 * wp1 * wp1 +
			wn1 * (-51001.0 * wn2 + 179098.0 * wp1 - 38947.0 * wp2) -
			3.0 * w0 * (99692.0 * wn1 - 22641.0 * wn2 + 99692.0 * wp1 - 22641.0 * wp2) +
			8209.0 * wn2 * wp2 - 51001.0 * wp1 * wp2 + 6908.0 * wp2 * wp2);

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

void weno_5th_ao_right(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double h)
{
	double dhi = 0.85;
	double dlo = 0.85;
	//-- - parameter of WENO-- -
	double beta[4], d[4], ww[4], alpha[4];
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

		beta[0] = 13.0 / 12.0 * pow((wn2 - 2.0 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4.0 * wn1 + 3.0 * w0), 2);
		beta[1] = 13.0 / 12.0 * pow((wn1 - 2.0 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
		beta[2] = 13.0 / 12.0 * pow((w0 - 2.0 * wp1 + wp2), 2) + 0.25 * pow((3.0 * w0 - 4.0 * wp1 + wp2), 2);

		beta[3] = (1.0 / 5040.0) * (231153.0 * w0 * w0 + 104963.0 * wn1 * wn1 + 6908.0 * wn2 * wn2 -
			38947.0 * wn2 * wp1 + 104963.0 * wp1 * wp1 +
			wn1 * (-51001.0 * wn2 + 179098.0 * wp1 - 38947.0 * wp2) -
			3.0 * w0 * (99692.0 * wn1 - 22641.0 * wn2 + 99692.0 * wp1 - 22641.0 * wp2) +
			8209.0 * wn2 * wp2 - 51001.0 * wp1 * wp2 + 6908.0 * wp2 * wp2);
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

double Calculate_alpha_k_1d(double* prim_left, double* prim_right)
{
	double A, Ma_left_normal, Ma_right_normal;

	Ma_left_normal = prim_left[1] / sqrt(Gamma * prim_left[2] / prim_left[0]);
	Ma_right_normal = prim_right[1] / sqrt(Gamma * prim_right[2] / prim_right[0]);

	A = abs(prim_left[2] - prim_right[2]) / prim_left[2] + abs(prim_left[2] - prim_right[2]) / prim_right[2];
	A = A + pow(Ma_left_normal - Ma_right_normal, 2);
	return A;
}

void Update_alpha(Interface1d* interfaces, Fluid1d* fluids, Block1d block)
{
	double prim_left_left[4], prim_left_right[4], prim_right_left[4], prim_right_right[4];
	double alpha_left, alpha_right;
	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		// here two gaussian points are used by default
		Convar_to_primvar_1D(prim_left_left, interfaces[i].left.convar);
		Convar_to_primvar_1D(prim_left_right, interfaces[i].right.convar);
		Convar_to_primvar_1D(prim_right_left, interfaces[i + 1].left.convar);
		Convar_to_primvar_1D(prim_right_right, interfaces[i + 1].right.convar);
		
		alpha_left = Calculate_alpha_k_1d(prim_left_left, prim_left_right);
		alpha_right = Calculate_alpha_k_1d(prim_right_left, prim_right_right);
		alpha_left = alpha_left < 2.0 ? 0.0 : alpha_left;
		alpha_right = alpha_right < 2.0 ? 0.0 : alpha_right;
		fluids[i].alpha = alpha_left + alpha_right;
	}
}

void WENO5_AO_with_DF(Point1d& left, Point1d& right, Fluid1d* fluids, Block1d block)
{
	// discontiunity feedback factor
	double df[4];
	double sum_df[4];
	sum_df[0] = fluids[-2].alpha + fluids[0].alpha;
	sum_df[1] = fluids[-1].alpha + fluids[1].alpha;
	sum_df[2] = fluids[0].alpha + fluids[2].alpha;
	sum_df[3] = fluids[-2].alpha + fluids[0].alpha + fluids[2].alpha;
	for (int i = 0; i < 4; i++)
	{
		df[i] = sum_df[i] < 2.0 ? 1.0 : 2.0 / sum_df[i];
	}
	//Note: function by WENO5_AO reconstruction
	double wn2[3]; Copy_Array(wn2, fluids[-2].convar, 3);
	double wn1[3]; Copy_Array(wn1, fluids[-1].convar, 3);
	double w0[3];  Copy_Array(w0, fluids[0].convar, 3);
	double wp1[3]; Copy_Array(wp1, fluids[1].convar, 3);
	double wp2[3]; Copy_Array(wp2, fluids[2].convar, 3);
	double tmp;
	//non-uniform grid was treated as uniform grid
	double  h = fluids[0].dx;

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 3; i++)
		{
			weno_5th_ao_with_df_right(right.convar[i], right.der1[i], tmp, wn2[i], wn1[i], w0[i], wp1[i], wp2[i], df, h);
			weno_5th_ao_with_df_left(left.convar[i], left.der1[i], tmp, wn2[i], wn1[i], w0[i], wp1[i], wp2[i], df, h);
		}
	}

	if (reconstruction_variable == characteristic)
	{

		double ren3[3], ren2[3], ren1[3], re0[3], rep1[3], rep2[3], rep3[3];
		double var[3], der1[3], der2[3];

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

		Convar_to_char1D(ren2, base_left, wn2);
		Convar_to_char1D(ren1, base_left, wn1);
		Convar_to_char1D(re0, base_left, w0);
		Convar_to_char1D(rep1, base_left, wp1);
		Convar_to_char1D(rep2, base_left, wp2);

		// left_reconstruction
		for (int i = 0; i < 3; i++)
		{
			weno_5th_ao_with_df_left(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], df, h);
		}
		Char_to_convar1D(left.convar, base_left, var);
		Char_to_convar1D(left.der1, base_left, der1);
	
		// right reconstruction

		Convar_to_char1D(ren2, base_right, wn2);
		Convar_to_char1D(ren1, base_right, wn1);
		Convar_to_char1D(re0, base_right, w0);
		Convar_to_char1D(rep1, base_right, wp1);
		Convar_to_char1D(rep2, base_right, wp2);

		for (int i = 0; i < 3; i++)
		{
			weno_5th_ao_with_df_right(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], df, h);
		}
		Char_to_convar1D(right.convar, base_right, var);
		Char_to_convar1D(right.der1, base_right, der1);

	}

	Check_Order_Reduce(left, right, fluids[0]);

}

void weno_5th_ao_with_df_left(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double* df, double h)
{
	double dhi = 0.85;
	double dlo = 0.85;
	//-- - parameter of WENO-- -
	double beta[4], d[4], ww[4], alpha[4];
	double epsilonW = 1e-6;
	//-- - intermediate parameter-- -
	double p[4], px[4], pxx[4], tempvar;
	double sum_alpha;

	//three small stencil
	d[0] = (1.0 - dhi) * (1.0 - dlo) / 2.0;
	d[1] = (1.0 - dhi) * dlo;
	d[2] = (1.0 - dhi) * (1.0 - dlo) / 2.0;
	//one big stencil
	d[3] = dhi;

	beta[0] = 13.0 / 12.0 * pow((wn2 - 2.0 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4.0 * wn1 + 3.0 * w0), 2);
	beta[1] = 13.0 / 12.0 * pow((wn1 - 2.0 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
	beta[2] = 13.0 / 12.0 * pow((w0 - 2.0 * wp1 + wp2), 2) + 0.25 * pow((3.0 * w0 - 4.0 * wp1 + wp2), 2);
	beta[3] = 1.0 / 6.0 * (beta[0] + 4.0 * beta[1] + beta[2]) + abs(beta[0] - beta[2]);
	double tau5 = 1.0 / 3.0 * (abs(beta[3] - beta[0]) + abs(beta[3] - beta[1]) + abs(beta[3] - beta[2]));

	sum_alpha = 0.0;
	for (int i = 0; i < 4; i++)
	{
		double global_div = tau5 / (beta[i] + epsilonW);
		alpha[i] = d[i] * (1.0 + global_div * global_div);
		sum_alpha += alpha[i];
	}

	for (int k = 0; k < 4; k++)
	{
		ww[k] = alpha[k] / sum_alpha;
	}
	//-- - candidate polynomial-- -
	p[0] = w0 + df[0] * (-1.0 / 6.0 * wn2 + 5.0 / 6.0 * wn1 + 1.0 / 3.0 * w0 - w0);
	p[1] = w0 + df[1] * (1.0 / 3.0 * wn1 + 5.0 / 6.0 * w0 - 1.0 / 6.0 * wp1 - w0);
	p[2] = w0 + df[2] * (11.0 / 6.0 * w0 - 7.0 / 6.0 * wp1 + 1.0 / 3.0 * wp2 - w0);
	p[3] = w0 + df[3] * ((1.0 / 60.0) * (47.0 * w0 + 27.0 * wn1 - 3.0 * wn2 - 13.0 * wp1 + 2.0 * wp2) - w0);

	px[0] = df[0] * (w0 - wn1) / h;
	px[1] = df[1] * (w0 - wn1) / h;
	px[2] = -df[2] * ((2.0 * w0 - 3.0 * wp1 + wp2) / h);
	px[3] = df[3] * (15.0 * w0 - 15.0 * wn1 + wn2 - wp1) / (12.0 * h);

	// double b2, c2, b4, c4, d4, e4, x = -1.0;
	// b2 = wn2 - 3.0 * wn1 + 2.0 * w0;
	// c2 = 0.5 * wn2 - wn1 + 0.5 * w0;
	// p[0] = w0 + df[0] * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	// px[0] = df[0] * (b2 + 2.0 * c2 * x) / h;
	// // 3th order polynomial 2
	// b2 = wp1 - w0;
	// c2 = 0.5 * wn1 - w0 + 0.5 * wp1;
	// p[1] = w0 + df[1] * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	// px[1] = df[1] * (b2 + 2.0 * c2 * x) / h;
	// // 3th order polynomial 3
	// b2 = wp1 - w0;
	// c2 = 0.5 * w0 - wp1 + 0.5 * wp2;
	// p[2] = w0 + df[2] * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	// px[2] = df[2] * (b2 + 2.0 * c2 * x) / h;
	// // 5th order polynomial
	// // a4 = (2 * wn2 - 13 * wn1 + 47 * w0 + 27 * wp1 - 3 * wp2) / 60;
	// b4 = (wn1 - 15 * w0 + 15 * wp1 - wp2) / 12.0;
	// c4 = (-wn2 + 6 * wn1 - 8 * w0 + 2 * wp1 + wp2) / 8.0;
	// d4 = (-wn1 + 3 * w0 - 3 * wp1 + wp2) / 6.0;
	// e4 = (wn2 - 4 * wn1 + 6 * w0 - 4 * wp1 + wp2) / 24.0;
	// p[3] = w0 +
	// 	df[3] * (b4 * (x + 1.0 / 2.0) + c4 * (x * x - 1.0 / 3.0) + d4 * (x * x * x + 1.0 / 4.0) + e4 * (x * x * x * x - 1.0 / 5.0));
	// px[3] = df[3] * (b4 + 2.0 * c4 * x + 3.0 * d4 * x * x + 4.0 * e4 * x * x * x) / h;

	//-- - combination-- -
	var = 0.0;
	der1 = 0.0;
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
	}
}

void weno_5th_ao_with_df_right(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double* df, double h)
{
	double dhi = 0.85;
	double dlo = 0.85;
	//-- - parameter of WENO-- -
	double beta[4], d[4], ww[4], alpha[4];
	double epsilonW = 1e-6;

	//-- - intermediate parameter-- -
	double p[4], px[4], pxx[4], tempvar;
	double sum_alpha;

	//three small stencil
	d[0] = (1 - dhi) * (1 - dlo) / 2.0;
	d[1] = (1 - dhi) * dlo;
	d[2] = (1 - dhi) * (1 - dlo) / 2.0;
	//one big stencil
	d[3] = dhi;

	beta[0] = 13.0 / 12.0 * pow((wn2 - 2.0 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4.0 * wn1 + 3.0 * w0), 2);
	beta[1] = 13.0 / 12.0 * pow((wn1 - 2.0 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
	beta[2] = 13.0 / 12.0 * pow((w0 - 2.0 * wp1 + wp2), 2) + 0.25 * pow((3.0 * w0 - 4.0 * wp1 + wp2), 2);
	beta[3] = 1.0 / 6.0 * (beta[0] + 4.0 * beta[1] + beta[2]) + abs(beta[0] - beta[2]);
	double tau5 = 1.0 / 3.0 * (abs(beta[3] - beta[0]) + abs(beta[3] - beta[1]) + abs(beta[3] - beta[2]));

	sum_alpha = 0.0;
	for (int i = 0; i < 4; i++)
	{
		double global_div = tau5 / (beta[i] + epsilonW);
		alpha[i] = d[i] * (1.0 + global_div * global_div);
		sum_alpha += alpha[i];
	}

	for (int k = 0; k < 4; k++)
	{
		ww[k] = alpha[k] / sum_alpha;
	}
	//-- - candidate polynomial-- -
	p[0] = w0 + df[0] * (1.0 / 3.0 * wn2 - 7.0 / 6.0 * wn1 + 11.0 / 6.0 * w0 - w0);
	p[1] = w0 + df[1] * (-1.0 / 6.0 * wn1 + 5.0 / 6.0 * w0 + 1.0 / 3.0 * wp1 - w0);
	p[2] = w0 + df[2] * (1.0 / 3.0 * w0 + 5.0 / 6.0 * wp1 - 1.0 / 6.0 * wp2 - w0);
	p[3] = w0 + df[3] * ((1.0 / 60.0) * (47.0 * w0 - 13.0 * wn1 + 2.0 * wn2 + 27.0 * wp1 - 3.0 * wp2) - w0);

	px[0] = df[0] * (2.0 * w0 - 3.0 * wn1 + wn2) / h;
	px[1] = df[1] * (-w0 + wp1) / h;
	px[2] = df[2] * (-w0 + wp1) / h;
	px[3] = df[3] * (-15.0 * w0 + wn1 + 15.0 * wp1 - wp2) / (12.0 * h);

	// double b2, c2, b4, c4, d4, e4, x = 0.0;
	// b2 = wn2 - 3.0 * wn1 + 2.0 * w0;
	// c2 = 0.5 * wn2 - wn1 + 0.5 * w0;
	// p[0] = w0 + df[0] * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	// px[0] = df[0] * (b2 + 2.0 * c2 * x) / h;
	// // 3th order polynomial 2
	// b2 = wp1 - w0;
	// c2 = 0.5 * wn1 - w0 + 0.5 * wp1;
	// p[1] = w0 + df[1] * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	// px[1] = df[1] * (b2 + 2.0 * c2 * x) / h;
	// // 3th order polynomial 3
	// b2 = wp1 - w0;
	// c2 = 0.5 * w0 - wp1 + 0.5 * wp2;
	// p[2] = w0 + df[2] * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	// px[2] = df[2] * (b2 + 2.0 * c2 * x) / h;
	// // 5th order polynomial
	// // a4 = (2 * wn2 - 13 * wn1 + 47 * w0 + 27 * wp1 - 3 * wp2) / 60;
	// b4 = (wn1 - 15 * w0 + 15 * wp1 - wp2) / 12.0;
	// c4 = (-wn2 + 6 * wn1 - 8 * w0 + 2 * wp1 + wp2) / 8.0;
	// d4 = (-wn1 + 3 * w0 - 3 * wp1 + wp2) / 6.0;
	// e4 = (wn2 - 4 * wn1 + 6 * w0 - 4 * wp1 + wp2) / 24.0;
	// p[3] = w0 +
	// 	df[3] * (b4 * (x + 1.0 / 2.0) + c4 * (x * x - 1.0 / 3.0) + d4 * (x * x * x + 1.0 / 4.0) + e4 * (x * x * x * x - 1.0 / 5.0));
	// px[3] = df[3] * (b4 + 2.0 * c4 * x + 3.0 * d4 * x * x + 4.0 * e4 * x * x * x) / h;

	//-- - combination-- -
	var = 0.0;
	der1 = 0.0;
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
	}
}

void WENO7_AO_with_DF(Point1d& left, Point1d& right, Fluid1d* fluids, Block1d block)
{
	// discontiunity feedback factor
	double df[5];
	double sum_df[5];
	sum_df[0] = fluids[-2].alpha + fluids[0].alpha;
	sum_df[1] = fluids[-1].alpha + fluids[1].alpha;
	sum_df[2] = fluids[0].alpha + fluids[2].alpha;
	sum_df[3] = fluids[-2].alpha + fluids[0].alpha + fluids[2].alpha;
	sum_df[4] = fluids[-3].alpha + fluids[-1].alpha + fluids[1].alpha + fluids[3].alpha;
	for (int i = 0; i < 5; i++)
	{
		df[i] = sum_df[i] < 2.0 ? 1.0 : 2.0 / sum_df[i];
	}
	//Note: function by WENO7_AO reconstruction
	double wn3[3]; Copy_Array(wn3, fluids[-3].convar, 3);
	double wn2[3]; Copy_Array(wn2, fluids[-2].convar, 3);
	double wn1[3]; Copy_Array(wn1, fluids[-1].convar, 3);
	double w0[3];  Copy_Array(w0, fluids[0].convar, 3);
	double wp1[3]; Copy_Array(wp1, fluids[1].convar, 3);
	double wp2[3]; Copy_Array(wp2, fluids[2].convar, 3);
	double wp3[3]; Copy_Array(wp3, fluids[3].convar, 3);
	double tmp;
	//non-uniform grid was treated as uniform grid
	double  h = fluids[0].dx;

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 3; i++)
		{
			weno_7th_ao_with_df_right(right.convar[i], right.der1[i], tmp, wn3[i], wn2[i], wn1[i], w0[i], wp1[i], wp2[i], wp3[i], df, h);
			weno_7th_ao_with_df_left(left.convar[i], left.der1[i], tmp, wn3[i], wn2[i], wn1[i], w0[i], wp1[i], wp2[i], wp3[i], df, h);
		}
	}

	if (reconstruction_variable == characteristic)
	{

		double ren3[3], ren2[3], ren1[3], re0[3], rep1[3], rep2[3], rep3[3];
		double var[3], der1[3], der2[3];

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
		Convar_to_char1D(ren3, base_left, wn3);
		Convar_to_char1D(ren2, base_left, wn2);
		Convar_to_char1D(ren1, base_left, wn1);
		Convar_to_char1D(re0, base_left, w0);
		Convar_to_char1D(rep1, base_left, wp1);
		Convar_to_char1D(rep2, base_left, wp2);
		Convar_to_char1D(rep3, base_left, wp3);

		// left_reconstruction
		for (int i = 0; i < 3; i++)
		{
			weno_7th_ao_with_df_left(var[i], der1[i], der2[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], df, h);
		}
		Char_to_convar1D(left.convar, base_left, var);
		Char_to_convar1D(left.der1, base_left, der1);
	
		// right reconstruction
		Convar_to_char1D(ren3, base_right, wn3);
		Convar_to_char1D(ren2, base_right, wn2);
		Convar_to_char1D(ren1, base_right, wn1);
		Convar_to_char1D(re0, base_right, w0);
		Convar_to_char1D(rep1, base_right, wp1);
		Convar_to_char1D(rep2, base_right, wp2);
		Convar_to_char1D(rep3, base_right, wp3);

		for (int i = 0; i < 3; i++)
		{
			weno_7th_ao_with_df_right(var[i], der1[i], der2[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], df, h);
		}
		Char_to_convar1D(right.convar, base_right, var);
		Char_to_convar1D(right.der1, base_right, der1);

	}

	Check_Order_Reduce(left, right, fluids[0]);

}

void weno_7th_ao_with_df_left(double& var, double& der1, double& der2, double wn3, double wn2, double wn1, double w0, double wp1, double wp2, double wp3, double* df, double h)
{
	double dhi = 0.85;
	double dlo = 0.85;
	double davg= 0.85;
	// P(7,5,3) -- one 7th-order stencil, one 5th-order stencil, three 3rd-order stencil
	//-- - parameter of WENO-- -
	double beta[5], d[5], ww[5], alpha[5];
	double epsilonW = 1e-6;
	//-- - intermediate parameter-- -
	double p[5], px[5], pxx[5], tempvar;
	double sum_alpha;

	//three 3rd-order stencil
	d[0] = (1.0 - dhi) * (1.0 - dlo) * (1.0 - davg) / 2.0;
	d[1] = (1.0 - dhi) * (1.0 - davg) * dlo;
	d[2] = (1.0 - dhi) * (1.0 - dlo) * (1.0 - davg) / 2.0;
	//one 5th-order stencil
	d[3] = (1.0 - dhi) * davg;
	//one 7th-order stencil
	d[4] = dhi;

	beta[0] = 13.0 / 12.0 * pow((wn2 - 2.0 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4.0 * wn1 + 3.0 * w0), 2);
	beta[1] = 13.0 / 12.0 * pow((wn1 - 2.0 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
	beta[2] = 13.0 / 12.0 * pow((w0 - 2.0 * wp1 + wp2), 2) + 0.25 * pow((3.0 * w0 - 4.0 * wp1 + wp2), 2);
	beta[3] = 1.0 / 6.0 * (beta[0] + 4.0 * beta[1] + beta[2]) + abs(beta[0] - beta[2]);

	//For 7th-order stencil, the smoothness indicator
	double a1, a2, a3;
	double beta4[4];
	a1 = (-19 * wn3 + 87 * wn2 - 177 * wn1 + 109 * w0) / 60.0;
	a2 = (-wn3 + 4 * wn2 - 5 * wn1 + 2 * w0) / 2;
	a3 = (-wn3 + 3 * wn2 - 3 * wn1 + w0) / 6.0;
	beta4[0] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (11 * wn2 - 63 * wn1 + 33 * w0 + 19 * wp1) / 60.0;
	a2 = (wn1 - 2 * w0 + wp1) / 2;
	a3 = (-wn2 + 3 * wn1 - 3 * w0 + wp1) / 6.0;
	beta4[1] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (-19 * wn1 - 33 * w0 + 63 * wp1 - 11 * wp2) / 60.0;
	a2 = (wn1 - 2 * w0 + wp1) / 2;
	a3 = (-wn1 + 3 * w0 - 3 * wp1 + wp2) / 6.0;
	beta4[2] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (-109 * w0 + 177 * wp1 - 87 * wp2 + 19 * wp3) / 60.0;
	a2 = (2 * w0 - 5 * wp1 + 4 * wp2 - wp3) / 2;
	a3 = (-w0 + 3 * wp1 - 3 * wp2 + wp3) / 6.0;
	beta4[3] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	beta[4] = 1.0 / 20.0 * (beta4[0] + 9.0 * beta4[1] + 9.0 * beta4[2] + beta4[3])
			+ abs(beta4[0] - beta4[1] - beta4[2] + beta4[3]);
	// double a1, a2, a3, a4, a5, a6;
	// a1 = (-191 * wn3 + 1688* wn2 - 7843 * wn1 + 7843 * wp1 - 1688 * wp2 + 191 * wp3) / 10080.0;
	// a2 = (79 * wn3 - 1014 * wn2 + 8385 * wn1 - 14900 * w0 + 8385 * wp1 - 1014 * wp2 + 79 * wp3) / 10080.0;
	// a3 = (5 * wn3 - 38 * wn2 + 61 * wn1 - 61 * wp1 + 38 * wp2 - 5 * wp3) / 216.0;
	// a4 = (-13 * wn3 + 144 * wn2 - 459 * wn1 + 656 * w0 - 459 * wp1 + 144 * wp2 - 13 * wp3) / 1584.0;
	// a5 = (-wn3 + 4 * wn2 - 5 * wn1 + 5 * wp1 - 4 * wp2 + wp3) / 240.0;
	// a6 = (wn3 - 6 * wn2 + 15 * wn1 - 20 * w0 + 15 * wp1 - 6 * wp2 + wp3) / 720.0;
	// beta[4] = (a1 + a3 / 10.0 + a5 / 126.0) * (a1 + a3 / 10.0 + a5 / 126.0) + 13.0 / 3.0 * (a2 + 123 / 455 * a4 + 85 / 2002 * a6)* (a2 + 123 / 455 * a4 + 85 / 2002 * a6)
	// 		+ 781 / 20 * (a3 + 26045 / 49203 * a5) * (a3 + 26045 / 49203 * a5)
	// 		+ 1421461 / 2275 * (a4 + 81596225 / 93816426 * a6) * (a4 + 81596225 / 93816426 * a6)
	// 		+ 21520059541 / 1377684 * a5 * a5 + 15510384942580921 / 27582029244 * a6 * a6;
	double tau5 = 0.25 * (abs(beta[4] - beta[0]) + abs(beta[4] - beta[1]) + abs(beta[4] - beta[2]) + abs(beta[4] - beta[3]));
	sum_alpha = 0.0;
	for (int i = 0; i < 5; i++)
	{
		double global_div = tau5 / (beta[i] + epsilonW);
		alpha[i] = d[i] * (1.0 + global_div * global_div);
		sum_alpha += alpha[i];
	}

	for (int k = 0; k < 5; k++)
	{
		ww[k] = alpha[k] / sum_alpha;
	}
	//-- - candidate polynomial-- -
	p[0] = w0 + df[0] * (-1.0 / 6.0 * wn2 + 5.0 / 6.0 * wn1 + 1.0 / 3.0 * w0 - w0);
	p[1] = w0 + df[1] * (1.0 / 3.0 * wn1 + 5.0 / 6.0 * w0 - 1.0 / 6.0 * wp1 - w0);
	p[2] = w0 + df[2] * (11.0 / 6.0 * w0 - 7.0 / 6.0 * wp1 + 1.0 / 3.0 * wp2 - w0);
	p[3] = w0 + df[3] * ((1.0 / 60.0) * (47.0 * w0 + 27.0 * wn1 - 3.0 * wn2 - 13.0 * wp1 + 2.0 * wp2) - w0);
	p[4] = w0 + df[4] * (1.0 / 420.0 * (319.0 * w0 + 214.0 * wn1 - 38.0 * wn2 + 4.0 * wn3 - 101.0 * wp1 + 25.0 * wp2 - 3.0 * wp3) - w0);

	px[0] = df[0] * (w0 - wn1) / h;
	px[1] = df[1] * (w0 - wn1) / h;
	px[2] = -df[2] * ((2.0 * w0 - 3.0 * wp1 + wp2) / h);
	px[3] = df[3] * (15.0 * w0 - 15.0 * wn1 + wn2 - wp1) / (12.0 * h);
	px[4] = df[4] * (245.0 * w0 - 245.0 * wn1 + 25.0 * wn2 - 2.0 * wn3 - 25.0 * wp1 + 2.0 * wp2) / (180.0 * h);

	//-- - combination-- -
	var = 0.0;
	der1 = 0.0;
	double final_weight[5];
	final_weight[4] = ww[4] / d[4];
	for (int k = 0; k < 4; k++)
	{
		final_weight[k] = ww[k] - ww[4] / d[4] * d[k];
	}
	for (int k = 0; k < 5; k++)
	{
		var += final_weight[k] * p[k];
		der1 += final_weight[k] * px[k];
	}
}

void weno_7th_ao_with_df_right(double& var, double& der1, double& der2, double wn3, double wn2, double wn1, double w0, double wp1, double wp2, double wp3, double* df, double h)
{
	double dhi = 0.85;
	double dlo = 0.85;
	double davg= 0.85;
	// P(7,5,3) -- one 7th-order stencil, one 5th-order stencil, three 3rd-order stencil
	//-- - parameter of WENO-- -
	double beta[5], d[5], ww[5], alpha[5];
	double epsilonW = 1e-6;
	//-- - intermediate parameter-- -
	double p[5], px[5], pxx[5], tempvar;
	double sum_alpha;

	//three 3rd-order stencil
	d[0] = (1.0 - dhi) * (1.0 - dlo) * (1.0 - davg) / 2.0;
	d[1] = (1.0 - dhi) * (1.0 - davg) * dlo;
	d[2] = (1.0 - dhi) * (1.0 - dlo) * (1.0 - davg) / 2.0;
	//one 5th-order stencil
	d[3] = (1.0 - dhi) * davg;
	//one 7th-order stencil
	d[4] = dhi;

	beta[0] = 13.0 / 12.0 * pow((wn2 - 2.0 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4.0 * wn1 + 3.0 * w0), 2);
	beta[1] = 13.0 / 12.0 * pow((wn1 - 2.0 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
	beta[2] = 13.0 / 12.0 * pow((w0 - 2.0 * wp1 + wp2), 2) + 0.25 * pow((3.0 * w0 - 4.0 * wp1 + wp2), 2);
	// beta[3] = (1.0 / 5040.0) * (231153.0 * w0 * w0 + 104963.0 * wn1 * wn1 + 6908.0 * wn2 * wn2 -
	// 		38947.0 * wn2 * wp1 + 104963.0 * wp1 * wp1 +
	// 		wn1 * (-51001.0 * wn2 + 179098.0 * wp1 - 38947.0 * wp2) -
	// 		3.0 * w0 * (99692.0 * wn1 - 22641.0 * wn2 + 99692.0 * wp1 - 22641.0 * wp2) +
	// 		8209.0 * wn2 * wp2 - 51001.0 * wp1 * wp2 + 6908.0 * wp2 * wp2);
	beta[3] = 1.0 / 6.0 * (beta[0] + 4.0 * beta[1] + beta[2]) + abs(beta[0] - beta[2]);

	//For 7th-order stencil, the smoothness indicator
	double a1, a2, a3;
	double beta4[4];
	a1 = (-19 * wn3 + 87 * wn2 - 177 * wn1 + 109 * w0) / 60.0;
	a2 = (-wn3 + 4 * wn2 - 5 * wn1 + 2 * w0) / 2;
	a3 = (-wn3 + 3 * wn2 - 3 * wn1 + w0) / 6.0;
	beta4[0] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (11 * wn2 - 63 * wn1 + 33 * w0 + 19 * wp1) / 60.0;
	a2 = (wn1 - 2 * w0 + wp1) / 2;
	a3 = (-wn2 + 3 * wn1 - 3 * w0 + wp1) / 6.0;
	beta4[1] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (-19 * wn1 - 33 * w0 + 63 * wp1 - 11 * wp2) / 60.0;
	a2 = (wn1 - 2 * w0 + wp1) / 2;
	a3 = (-wn1 + 3 * w0 - 3 * wp1 + wp2) / 6.0;
	beta4[2] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (-109 * w0 + 177 * wp1 - 87 * wp2 + 19 * wp3) / 60.0;
	a2 = (2 * w0 - 5 * wp1 + 4 * wp2 - wp3) / 2;
	a3 = (-w0 + 3 * wp1 - 3 * wp2 + wp3) / 6.0;
	beta4[3] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	beta[4] = 1.0 / 20.0 * (beta4[0] + 9.0 * beta4[1] + 9.0 * beta4[2] + beta4[3])
			+ abs(beta4[0] - beta4[1] - beta4[2] + beta4[3]);
	// double a1, a2, a3, a4, a5, a6;
	// a1 = (-191 * wn3 + 1688* wn2 - 7843 * wn1 + 7843 * wp1 - 1688 * wp2 + 191 * wp3) / 10080.0;
	// a2 = (79 * wn3 - 1014 * wn2 + 8385 * wn1 - 14900 * w0 + 8385 * wp1 - 1014 * wp2 + 79 * wp3) / 10080.0;
	// a3 = (5 * wn3 - 38 * wn2 + 61 * wn1 - 61 * wp1 + 38 * wp2 - 5 * wp3) / 216.0;
	// a4 = (-13 * wn3 + 144 * wn2 - 459 * wn1 + 656 * w0 - 459 * wp1 + 144 * wp2 - 13 * wp3) / 1584.0;
	// a5 = (-wn3 + 4 * wn2 - 5 * wn1 + 5 * wp1 - 4 * wp2 + wp3) / 240.0;
	// a6 = (wn3 - 6 * wn2 + 15 * wn1 - 20 * w0 + 15 * wp1 - 6 * wp2 + wp3) / 720.0;
	// beta[4] = (a1 + a3 / 10.0 + a5 / 126.0) * (a1 + a3 / 10.0 + a5 / 126.0) + 13.0 / 3.0 * (a2 + 123 / 455 * a4 + 85 / 2002 * a6)* (a2 + 123 / 455 * a4 + 85 / 2002 * a6)
	// 		+ 781 / 20 * (a3 + 26045 / 49203 * a5) * (a3 + 26045 / 49203 * a5)
	// 		+ 1421461 / 2275 * (a4 + 81596225 / 93816426 * a6) * (a4 + 81596225 / 93816426 * a6)
	// 		+ 21520059541 / 1377684 * a5 * a5 + 15510384942580921 / 27582029244 * a6 * a6;
	double tau5 = 0.25 * (abs(beta[4] - beta[0]) + abs(beta[4] - beta[1]) + abs(beta[4] - beta[2]) + abs(beta[4] - beta[3]));
	sum_alpha = 0.0;
	for (int i = 0; i < 5; i++)
	{
		double global_div = tau5 / (beta[i] + epsilonW);
		alpha[i] = d[i] * (1.0 + global_div * global_div);
		sum_alpha += alpha[i];
	}

	for (int k = 0; k < 5; k++)
	{
		ww[k] = alpha[k] / sum_alpha;
	}
	//-- - candidate polynomial-- -
	// 3th order polynomial 1
	double b2, c2, b4, c4, d4, e4, x = 0.0;
	double b6, c6, d6, e6, f6, g6;
	p[0] = w0 + df[0] * (1.0 / 3.0 * wn2 - 7.0 / 6.0 * wn1 + 11.0 / 6.0 * w0 - w0);
	p[1] = w0 + df[1] * (-1.0 / 6.0 * wn1 + 5.0 / 6.0 * w0 + 1.0 / 3.0 * wp1 - w0);
	p[2] = w0 + df[2] * (1.0 / 3.0 * w0 + 5.0 / 6.0 * wp1 - 1.0 / 6.0 * wp2 - w0);
	p[3] = w0 + df[3] * ((1.0 / 60.0) * (47.0 * w0 - 13.0 * wn1 + 2.0 * wn2 + 27.0 * wp1 - 3.0 * wp2) - w0);
	p[4] = w0 + df[4] * (1.0 / 420.0 * (319.0 * w0 - 101.0 * wn1 + 25.0 * wn2 - 3.0 * wn3 + 214.0 * wp1 - 38.0 * wp2 + 4.0 * wp3) - w0);
	
	px[0] = df[0] * (2.0 * w0 - 3.0 * wn1 + wn2) / h;
	px[1] = df[1] * (-w0 + wp1) / h;
	px[2] = df[2] * (-w0 + wp1) / h;
	px[3] = df[3] * (-15.0 * w0 + wn1 + 15.0 * wp1 - wp2) / (12.0 * h);
	px[4] = df[4] * (-245.0 * w0 + 25.0 * wn1 - 2.0 * wn2 + 245.0 * wp1 - 25.0 * wp2 + 2.0 * wp3) / (180.0 * h);
	
	//-- - combination-- -
	var = 0.0;
	der1 = 0.0;
	double final_weight[5];
	final_weight[4] = ww[4] / d[4];
	for (int k = 0; k < 4; k++)
	{
		final_weight[k] = ww[k] - ww[4] / d[4] * d[k];
	}
	for (int k = 0; k < 5; k++)
	{
		var += final_weight[k] * p[k];
		der1 += final_weight[k] * px[k];
	}
}

void WENO9_AO_with_DF(Point1d& left, Point1d& right, Fluid1d* fluids, Block1d block)
{
	// discontiunity feedback factor
	double df[6];
	double sum_df[6];
	sum_df[0] = fluids[-2].alpha + fluids[0].alpha;
	sum_df[1] = fluids[-1].alpha + fluids[1].alpha;
	sum_df[2] = fluids[0].alpha + fluids[2].alpha;
	sum_df[3] = fluids[-2].alpha + fluids[0].alpha + fluids[2].alpha;
	sum_df[4] = fluids[-3].alpha + fluids[-1].alpha + fluids[1].alpha + fluids[3].alpha;
	sum_df[5] = fluids[-4].alpha + fluids[-2].alpha + fluids[0].alpha + fluids[2].alpha + fluids[4].alpha;
	for (int i = 0; i < 6; i++)
	{
		df[i] = sum_df[i] < 2.0 ? 1.0 : 2.0 / sum_df[i];
	}
	//Note: function by WENO9_AO reconstruction
	double wn4[3]; Copy_Array(wn4, fluids[-4].convar, 3);
	double wn3[3]; Copy_Array(wn3, fluids[-3].convar, 3);
	double wn2[3]; Copy_Array(wn2, fluids[-2].convar, 3);
	double wn1[3]; Copy_Array(wn1, fluids[-1].convar, 3);
	double w0[3];  Copy_Array(w0, fluids[0].convar, 3);
	double wp1[3]; Copy_Array(wp1, fluids[1].convar, 3);
	double wp2[3]; Copy_Array(wp2, fluids[2].convar, 3);
	double wp3[3]; Copy_Array(wp3, fluids[3].convar, 3);
	double wp4[3]; Copy_Array(wp4, fluids[4].convar, 3);
	double tmp;
	//non-uniform grid was treated as uniform grid
	double  h = fluids[0].dx;

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 3; i++)
		{
			weno_9th_ao_with_df_right(right.convar[i], right.der1[i], tmp, wn4[i], wn3[i], wn2[i], wn1[i], w0[i], wp1[i], wp2[i], wp3[i], wp4[i], df, h);
			weno_9th_ao_with_df_left(left.convar[i], left.der1[i], tmp, wn4[i], wn3[i], wn2[i], wn1[i], w0[i], wp1[i], wp2[i], wp3[i], wp4[i], df, h);
		}
	}

	if (reconstruction_variable == characteristic)
	{

		double ren4[3], ren3[3], ren2[3], ren1[3], re0[3], rep1[3], rep2[3], rep3[3], rep4[3];
		double var[3], der1[3], der2[3];

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
		Convar_to_char1D(ren4, base_left, wn4);
		Convar_to_char1D(ren3, base_left, wn3);
		Convar_to_char1D(ren2, base_left, wn2);
		Convar_to_char1D(ren1, base_left, wn1);
		Convar_to_char1D(re0, base_left, w0);
		Convar_to_char1D(rep1, base_left, wp1);
		Convar_to_char1D(rep2, base_left, wp2);
		Convar_to_char1D(rep3, base_left, wp3);
		Convar_to_char1D(rep4, base_left, wp4);

		// left_reconstruction
		for (int i = 0; i < 3; i++)
		{
			weno_9th_ao_with_df_left(var[i], der1[i], der2[i], ren4[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], rep4[i], df, h);
		}
		Char_to_convar1D(left.convar, base_left, var);
		Char_to_convar1D(left.der1, base_left, der1);
	
		// right reconstruction
		Convar_to_char1D(ren3, base_right, wn3);
		Convar_to_char1D(ren2, base_right, wn2);
		Convar_to_char1D(ren1, base_right, wn1);
		Convar_to_char1D(re0, base_right, w0);
		Convar_to_char1D(rep1, base_right, wp1);
		Convar_to_char1D(rep2, base_right, wp2);
		Convar_to_char1D(rep3, base_right, wp3);

		for (int i = 0; i < 3; i++)
		{
			weno_9th_ao_with_df_right(var[i], der1[i], der2[i], ren4[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], rep4[i], df, h);
		}
		Char_to_convar1D(right.convar, base_right, var);
		Char_to_convar1D(right.der1, base_right, der1);

	}

	Check_Order_Reduce(left, right, fluids[0]);

}

void weno_9th_ao_with_df_left(double& var, double& der1, double& der2, double wn4, double wn3, double wn2, double wn1, double w0, double wp1, double wp2, double wp3, double wp4, double* df, double h)
{
	double dhi = 0.85;
	double dlo = 0.85;
	double davg1= 0.85;
	double davg2 = 0.85;
	// P(9,5,3) -- one 9th-order stencil, one 5th-order stencil, three 3rd-order stencil
	//-- - parameter of WENO-- -
	double beta[6], d[6], ww[6], alpha[6];
	double epsilonW = 1e-6;
	//-- - intermediate parameter-- -
	double p[6], px[6], pxx[6], tempvar;
	double sum_alpha;

	//three 3rd-order stencil
	d[0] = (1.0 - dhi) * (1.0 - dlo) * (1.0 - davg1) * (1.0 - davg2) / 2.0;
	d[1] = (1.0 - dhi) * (1.0 - davg1) * (1.0 - davg2) * dlo;
	d[2] = (1.0 - dhi) * (1.0 - dlo) * (1.0 - davg1) * (1.0 - davg2) / 2.0;
	//one 5th-order stencil
	d[3] = (1.0 - dhi) * (1.0 - davg2) * davg1;
	//one 7th-order stencil
	d[4] = (1.0 - dhi) * davg2;
	//one 9th-order stencil
	d[5] = dhi;

	beta[0] = 13.0 / 12.0 * pow((wn2 - 2.0 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4.0 * wn1 + 3.0 * w0), 2);
	beta[1] = 13.0 / 12.0 * pow((wn1 - 2.0 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
	beta[2] = 13.0 / 12.0 * pow((w0 - 2.0 * wp1 + wp2), 2) + 0.25 * pow((3.0 * w0 - 4.0 * wp1 + wp2), 2);
	//beta[3] = 1.0 / 6.0 * (beta[0] + 4.0 * beta[1] + beta[2]) + abs(beta[0] - beta[2]);

	//For 7th-order stencil, the smoothness indicator
	double a1, a2, a3, a4;
	double beta4[4];
	a1 = (-19 * wn3 + 87 * wn2 - 177 * wn1 + 109 * w0) / 60.0;
	a2 = (-wn3 + 4 * wn2 - 5 * wn1 + 2 * w0) / 2;
	a3 = (-wn3 + 3 * wn2 - 3 * wn1 + w0) / 6.0;
	beta4[0] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (11 * wn2 - 63 * wn1 + 33 * w0 + 19 * wp1) / 60.0;
	a2 = (wn1 - 2 * w0 + wp1) / 2;
	a3 = (-wn2 + 3 * wn1 - 3 * w0 + wp1) / 6.0;
	beta4[1] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (-19 * wn1 - 33 * w0 + 63 * wp1 - 11 * wp2) / 60.0;
	a2 = (wn1 - 2 * w0 + wp1) / 2;
	a3 = (-wn1 + 3 * w0 - 3 * wp1 + wp2) / 6.0;
	beta4[2] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (-109 * w0 + 177 * wp1 - 87 * wp2 + 19 * wp3) / 60.0;
	a2 = (2 * w0 - 5 * wp1 + 4 * wp2 - wp3) / 2;
	a3 = (-w0 + 3 * wp1 - 3 * wp2 + wp3) / 6.0;
	beta4[3] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	beta[4] = 1.0 / 20.0 * (beta4[0] + 9.0 * beta4[1] + 9.0 * beta4[2] + beta4[3])
			+ abs(beta4[0] - beta4[1] - beta4[2] + beta4[3]);

	// For 9th-order stencil, the smoothness indicator
	double beta5[5];
	a1 = (-462 * wn1 + 336 * wn2 - 146 * wn3 + 27 * wn4 + 245 * w0) / 120.0;
	a2 = (-240 * wn1 + 262 * wn2 - 128 * wn3 + 25 * wn4 + 81 * w0) / 56.0;
	a3 = (-18 * wn1 + 24 * wn2 - 14 * wn3 + 3 * wn4 + 5 * w0) / 12.0;
	a4 = (-4 * wn1 + 6 * wn2 - 4 * wn3 + wn4 + w0) / 24.0;
	beta5[0] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13.0 / 3.0 * (a2 + 123.0 / 455.0 * a4) * (a2 + 123.0 / 455.0 * a4) 
		     + 781.0 / 20.0 * a3 * a3 + 1421461.0 / 2275.0 * a4 * a4;

	a1 = (-192 * wn1 + 66 * wn2 - 11 * wn3 + 110 * w0 + 27 * wp1) / 120;
	a2 = (10 * wn1 + 12 * wn2 - 3 * wn3 - 44 * w0 + 25 * wp1) / 56;
	a3 = (12 * wn1 - 6 * wn2 + wn3 - 10 * w0 + 3 * wp1) / 12;
	a4 = (6 * wn1 - 4 * wn2 + wn3 - 4 * w0 + wp1) / 24;
	beta5[1] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13.0 / 3.0 * (a2 + 123.0 / 455.0 * a4) * (a2 + 123.0 / 455.0 * a4) 
		     + 781.0 / 20.0 * a3 * a3 + 1421461.0 / 2275.0 * a4 * a4;

	a1 = (-82 * wn1 + 11 * wn2 + 82 * wp1 - 11 * wp2)/120;
	a2 = (40 * wn1 - 3 * wn2 - 74 * w0 + 40 * wp1 - 3 * wp2)/56;
	a3 = (2 * wn1 - wn2 - 2 * wp1 + wp2)/12;
	a4 = (-4 * wn1 + wn2 + 6 * w0 - 4 * wp1 + wp2)/24;
	beta5[2] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13.0 / 3.0 * (a2 + 123.0 / 455.0 * a4) * (a2 + 123.0 / 455.0 * a4) 
		     + 781.0 / 20.0 * a3 * a3 + 1421461.0 / 2275.0 * a4 * a4;
	
	a1 = (-27 * wn1 - 110 * w0 + 192 * wp1 - 66 * wp2 + 11 * wp3)/120;
	a2 = (25 * wn1 - 44 * w0 + 10 * wp1 + 12 * wp2 - 3 * wp3)/56;
	a3 = (-3 * wn1 + 10 * w0 - 12 * wp1 + 6 * wp2 - wp3)/12;
	a4 = (wn1 - 4 * w0 + 6 * wp1 - 4 * wp2 + wp3)/24;
	beta5[3] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13.0 / 3.0 * (a2 + 123.0 / 455.0 * a4) * (a2 + 123.0 / 455.0 * a4) 
		     + 781.0 / 20.0 * a3 * a3 + 1421461.0 / 2275.0 * a4 * a4;

	a1 = (-245 * w0 + 462 * wp1 - 336 * wp2 + 146 * wp3 - 27 * wp4)/120;
	a2 = (81 * w0 - 240 * wp1 + 262 * wp2 - 128 * wp3 + 25 * wp4)/56;
	a3 = (-5 * w0 + 18 * wp1 - 24 * wp2 + 14 * wp3 - 3 * wp4)/12;
	a4 = (w0 - 4 * wp1 + 6 * wp2 - 4 * wp3 + wp4)/24;
	beta5[4] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13.0 / 3.0 * (a2 + 123.0 / 455.0 * a4) * (a2 + 123.0 / 455.0 * a4) 
		     + 781.0 / 20.0 * a3 * a3 + 1421461.0 / 2275.0 * a4 * a4;
	
	// beta[3] represent the 5th-order smoothness indicator
	beta[3] = beta5[2];
	beta[5] = 1.0 / 70.0 * (beta5[0] + 16.0 * beta5[1] + 36.0 * beta5[2] + 16.0 * beta5[3] + beta5[4]) + abs(beta5[0] - beta5[1] - beta5[3] + beta5[4]);
			
	double tau5 = 0.2 * (abs(beta[5] - beta[0]) + abs(beta[5] - beta[1]) + abs(beta[5] - beta[2]) + abs(beta[5] - beta[3]) + abs(beta[5] - beta[4]));

	sum_alpha = 0.0;
	for (int i = 0; i < 6; i++)
	{
		double global_div = tau5 / (beta[i] + epsilonW);
		alpha[i] = d[i] * (1.0 + global_div * global_div);
		sum_alpha += alpha[i];
	}

	for (int k = 0; k < 6; k++)
	{
		ww[k] = alpha[k] / sum_alpha;
	}
	//-- - candidate polynomial-- -
	p[0] = w0 + df[0] * (-1.0 / 6.0 * wn2 + 5.0 / 6.0 * wn1 + 1.0 / 3.0 * w0 - w0);
	p[1] = w0 + df[1] * (1.0 / 3.0 * wn1 + 5.0 / 6.0 * w0 - 1.0 / 6.0 * wp1 - w0);
	p[2] = w0 + df[2] * (11.0 / 6.0 * w0 - 7.0 / 6.0 * wp1 + 1.0 / 3.0 * wp2 - w0);
	p[3] = w0 + df[3] * ((1.0 / 60.0) * (47.0 * w0 + 27.0 * wn1 - 3.0 * wn2 - 13.0 * wp1 + 2.0 * wp2) - w0);
	p[4] = w0 + df[4] * ((1.0 / 420.0) * (319.0 * w0 + 214.0 * wn1 - 38.0 * wn2 + 4.0 * wn3 - 101.0 * wp1 + 25.0 * wp2 - 3.0 * wp3) - w0);
	p[5] = w0 + df[5] * ((1.0 / 2520.0) * (1879.0 * w0 + 1375.0 * wn1 - 305.0 * wn2 + 55.0 * wn3 - 5.0 * wn4 - 641.0 * wp1 + 199.0 * wp2 - 41.0 * wp3 + 4.0 * wp4) - w0);
	
	px[0] = df[0] * (w0 - wn1) / h;
	px[1] = df[1] * (w0 - wn1) / h;
	px[2] = -df[2] * ((2.0 * w0 - 3.0 * wp1 + wp2) / h);
	px[3] = df[3] * (15.0 * w0 - 15.0 * wn1 + wn2 - wp1) / (12.0 * h);
	px[4] = df[4] * (245.0 * w0 - 245.0 * wn1 + 25.0 * wn2 - 2.0 * wn3 - 25.0 * wp1 + 2.0 * wp2) / (180.0 * h);
	px[5] = df[5] * (7175.0 * w0 - 7175.0 * wn1 + 889.0 * wn2 - 119.0 * wn3 + 9.0 * wn4 - 889.0 * wp1 + 119.0 * wp2 - 9.0 * wp3) / (5040.0 * h);
	//-- - combination-- -
	var = 0.0;
	der1 = 0.0;
	double final_weight[6];
	final_weight[5] = ww[5] / d[5];
	for (int k = 0; k < 5; k++)
	{
		final_weight[k] = ww[k] - ww[5] / d[5] * d[k];
	}

	for (int k = 0; k < 6; k++)
	{
		var += final_weight[k] * p[k];
		der1 += final_weight[k] * px[k];
	}
}

void weno_9th_ao_with_df_right(double& var, double& der1, double& der2, double wn4, double wn3, double wn2, double wn1, double w0, double wp1, double wp2, double wp3, double wp4, double* df, double h)
{
	double dhi = 0.85;
	double dlo = 0.85;
	double davg1= 0.85;
	double davg2 = 0.85;
	// P(9, 7,5,3) -- one 9th-order stencil, one 7th-order stencil, one 5th-order stencil, three 3rd-order stencil
	//-- - parameter of WENO-- -
	double beta[6], d[6], ww[6], alpha[6];
	double epsilonW = 1e-6;
	//-- - intermediate parameter-- -
	double p[6], px[6], pxx[6], tempvar;
	double sum_alpha;

	//three 3rd-order stencil
	d[0] = (1.0 - dhi) * (1.0 - dlo) * (1.0 - davg1) * (1.0 - davg2) / 2.0;
	d[1] = (1.0 - dhi) * (1.0 - davg1) * (1.0 - davg2) * dlo;
	d[2] = (1.0 - dhi) * (1.0 - dlo) * (1.0 - davg1) * (1.0 - davg2) / 2.0;
	//one 5th-order stencil
	d[3] = (1.0 - dhi) * (1.0 - davg2) * davg1;
	//one 7th-order stencil
	d[4] = (1.0 - dhi) * davg2;
	//one 9th-order stencil
	d[5] = dhi;

	beta[0] = 13.0 / 12.0 * pow((wn2 - 2.0 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4.0 * wn1 + 3.0 * w0), 2);
	beta[1] = 13.0 / 12.0 * pow((wn1 - 2.0 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
	beta[2] = 13.0 / 12.0 * pow((w0 - 2.0 * wp1 + wp2), 2) + 0.25 * pow((3.0 * w0 - 4.0 * wp1 + wp2), 2);
	//beta[3] = 1.0 / 6.0 * (beta[0] + 4.0 * beta[1] + beta[2]) + abs(beta[0] - beta[2]);

	//For 7th-order stencil, the smoothness indicator
	double a1, a2, a3, a4;
	double beta4[4];
	a1 = (-19 * wn3 + 87 * wn2 - 177 * wn1 + 109 * w0) / 60.0;
	a2 = (-wn3 + 4 * wn2 - 5 * wn1 + 2 * w0) / 2;
	a3 = (-wn3 + 3 * wn2 - 3 * wn1 + w0) / 6.0;
	beta4[0] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (11 * wn2 - 63 * wn1 + 33 * w0 + 19 * wp1) / 60.0;
	a2 = (wn1 - 2 * w0 + wp1) / 2;
	a3 = (-wn2 + 3 * wn1 - 3 * w0 + wp1) / 6.0;
	beta4[1] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (-19 * wn1 - 33 * w0 + 63 * wp1 - 11 * wp2) / 60.0;
	a2 = (wn1 - 2 * w0 + wp1) / 2;
	a3 = (-wn1 + 3 * w0 - 3 * wp1 + wp2) / 6.0;
	beta4[2] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (-109 * w0 + 177 * wp1 - 87 * wp2 + 19 * wp3) / 60.0;
	a2 = (2 * w0 - 5 * wp1 + 4 * wp2 - wp3) / 2;
	a3 = (-w0 + 3 * wp1 - 3 * wp2 + wp3) / 6.0;
	beta4[3] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	beta[4] = 1.0 / 20.0 * (beta4[0] + 9.0 * beta4[1] + 9.0 * beta4[2] + beta4[3])
			+ abs(beta4[0] - beta4[1] - beta4[2] + beta4[3]);
	
	// For 9th-order stencil, the smoothness indicator
	double beta5[5];
	a1 = (-462 * wn1 + 336 * wn2 - 146 * wn3 + 27 * wn4 + 245 * w0) / 120.0;
	a2 = (-240 * wn1 + 262 * wn2 - 128 * wn3 + 25 * wn4 + 81 * w0) / 56.0;
	a3 = (-18 * wn1 + 24 * wn2 - 14 * wn3 + 3 * wn4 + 5 * w0) / 12.0;
	a4 = (-4 * wn1 + 6 * wn2 - 4 * wn3 + wn4 + w0) / 24.0;
	beta5[0] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13.0 / 3.0 * (a2 + 123.0 / 455.0 * a4) * (a2 + 123.0 / 455.0 * a4) 
		     + 781.0 / 20.0 * a3 * a3 + 1421461.0 / 2275.0 * a4 * a4;

	a1 = (-192 * wn1 + 66 * wn2 - 11 * wn3 + 110 * w0 + 27 * wp1) / 120;
	a2 = (10 * wn1 + 12 * wn2 - 3 * wn3 - 44 * w0 + 25 * wp1) / 56;
	a3 = (12 * wn1 - 6 * wn2 + wn3 - 10 * w0 + 3 * wp1) / 12;
	a4 = (6 * wn1 - 4 * wn2 + wn3 - 4 * w0 + wp1) / 24;
	beta5[1] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13.0 / 3.0 * (a2 + 123.0 / 455.0 * a4) * (a2 + 123.0 / 455.0 * a4) 
		     + 781.0 / 20.0 * a3 * a3 + 1421461.0 / 2275.0 * a4 * a4;

	a1 = (-82 * wn1 + 11 * wn2 + 82 * wp1 - 11 * wp2)/120;
	a2 = (40 * wn1 - 3 * wn2 - 74 * w0 + 40 * wp1 - 3 * wp2)/56;
	a3 = (2 * wn1 - wn2 - 2 * wp1 + wp2)/12;
	a4 = (-4 * wn1 + wn2 + 6 * w0 - 4 * wp1 + wp2)/24;
	beta5[2] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13.0 / 3.0 * (a2 + 123.0 / 455.0 * a4) * (a2 + 123.0 / 455.0 * a4) 
		     + 781.0 / 20.0 * a3 * a3 + 1421461.0 / 2275.0 * a4 * a4;
	
	a1 = (-27 * wn1 - 110 * w0 + 192 * wp1 - 66 * wp2 + 11 * wp3)/120;
	a2 = (25 * wn1 - 44 * w0 + 10 * wp1 + 12 * wp2 - 3 * wp3)/56;
	a3 = (-3 * wn1 + 10 * w0 - 12 * wp1 + 6 * wp2 - wp3)/12;
	a4 = (wn1 - 4 * w0 + 6 * wp1 - 4 * wp2 + wp3)/24;
	beta5[3] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13.0 / 3.0 * (a2 + 123.0 / 455.0 * a4) * (a2 + 123.0 / 455.0 * a4) 
		     + 781.0 / 20.0 * a3 * a3 + 1421461.0 / 2275.0 * a4 * a4;

	a1 = (-245 * w0 + 462 * wp1 - 336 * wp2 + 146 * wp3 - 27 * wp4)/120;
	a2 = (81 * w0 - 240 * wp1 + 262 * wp2 - 128 * wp3 + 25 * wp4)/56;
	a3 = (-5 * w0 + 18 * wp1 - 24 * wp2 + 14 * wp3 - 3 * wp4)/12;
	a4 = (w0 - 4 * wp1 + 6 * wp2 - 4 * wp3 + wp4)/24;
	beta5[4] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13.0 / 3.0 * (a2 + 123.0 / 455.0 * a4) * (a2 + 123.0 / 455.0 * a4) 
		     + 781.0 / 20.0 * a3 * a3 + 1421461.0 / 2275.0 * a4 * a4;
	
	// beta[3] represent the 5th-order smoothness indicator
	beta[3] = beta5[2];
	beta[5] = 1.0 / 70.0 * (beta5[0] + 16.0 * beta5[1] + 36.0 * beta5[2] + 16.0 * beta5[3] + beta5[4]) + abs(beta5[0] - beta5[1] - beta5[3] + beta5[4]);
			
	double tau5 = 0.2 * (abs(beta[5] - beta[0]) + abs(beta[5] - beta[1]) + abs(beta[5] - beta[2]) + abs(beta[5] - beta[3]) + abs(beta[5] - beta[4]));


	sum_alpha = 0.0;
	for (int i = 0; i < 6; i++)
	{
		double global_div = tau5 / (beta[i] + epsilonW);
		alpha[i] = d[i] * (1.0 + global_div * global_div);
		sum_alpha += alpha[i];
	}

	for (int k = 0; k < 6; k++)
	{
		ww[k] = alpha[k] / sum_alpha;
	}
	//-- - candidate polynomial-- -
	double b9, c9, d9, e9, f9, g9, h9, i9, x = 0.0;
	p[0] = w0 + df[0] * (1.0 / 3.0 * wn2 - 7.0 / 6.0 * wn1 + 11.0 / 6.0 * w0 - w0);
	p[1] = w0 + df[1] * (-1.0 / 6.0 * wn1 + 5.0 / 6.0 * w0 + 1.0 / 3.0 * wp1 - w0);
	p[2] = w0 + df[2] * (1.0 / 3.0 * w0 + 5.0 / 6.0 * wp1 - 1.0 / 6.0 * wp2 - w0);
	p[3] = w0 + df[3] * ((1.0 / 60.0) * (47.0 * w0 - 13.0 * wn1 + 2.0 * wn2 + 27.0 * wp1 - 3.0 * wp2) - w0);
	p[4] = w0 + df[4] * (1.0 / 420.0 * (319.0 * w0 - 101.0 * wn1 + 25.0 * wn2 - 3.0 * wn3 + 214.0 * wp1 - 38.0 * wp2 + 4.0 * wp3) - w0);
	p[5] = w0 + df[5] * ((1.0 / 2520) * (1879.0 * w0 - 641.0 * wn1 + 199.0 * wn2 - 41.0 * wn3 + 4.0 * wn4 + 1375.0 * wp1 - 305.0 * wp2 + 55.0 * wp3 - 5.0 * wp4) - w0);
	
	px[0] = df[0] * (2.0 * w0 - 3.0 * wn1 + wn2) / h;
	px[1] = df[1] * (-w0 + wp1) / h;
	px[2] = df[2] * (-w0 + wp1) / h;
	px[3] = df[3] * (-15.0 * w0 + wn1 + 15.0 * wp1 - wp2) / (12.0 * h);
	px[4] = df[4] * (-245.0 * w0 + 25.0 * wn1 - 2.0 * wn2 + 245.0 * wp1 - 25.0 * wp2 + 2.0 * wp3) / (180.0 * h);
	px[5] = df[5] * (-7175.0 * w0 + 889.0 * wn1 - 119.0 * wn2 + 9.0 * wn3 + 7175.0 * wp1 - 889.0 * wp2 + 119.0 * wp3 - 9.0 * wp4) / (5040.0 * h);
	//-- - combination-- -
	var = 0.0;
	der1 = 0.0;
	double final_weight[6];
	final_weight[5] = ww[5] / d[5];
	for (int k = 0; k < 5; k++)
	{
		final_weight[k] = ww[k] - ww[5] / d[5] * d[k];
	}

	for (int k = 0; k < 6; k++)
	{
		var += final_weight[k] * p[k];
		der1 += final_weight[k] * px[k];
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

// two-dimensional problem
Reconstruction_within_Cell_2D_normal cellreconstruction_2D_normal = WENO5_AO_normal; //initilization
Reconstruction_within_Cell_2D_tangent cellreconstruction_2D_tangent = WENO5_AO_tangent; //initilization
Reconstruction_forG0_2D_normal g0reconstruction_2D_normal = Center_do_nothing_normal; //initilization
Reconstruction_forG0_2D_tangent g0reconstruction_2D_tangent = Center_all_collision_multi; //initilization

void Reconstruction_within_cell(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{
#pragma omp parallel  for
	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			cellreconstruction_2D_normal
			(xinterfaces[i * (block.ny + 1) + j], xinterfaces[(i + 1) * (block.ny + 1) + j],
				yinterfaces[i * (block.ny + 1) + j], yinterfaces[i * (block.ny + 1) + j + 1], &fluids[i * block.ny + j], block);
		}
	}
#pragma omp parallel  for
	for (int i = block.ghost - 1; i < block.nx - block.ghost + 1; i++)
	{
		for (int j = block.ghost - 1; j < block.ny - block.ghost + 1; j++)
		{
			cellreconstruction_2D_tangent
			(&xinterfaces[i * (block.ny + 1) + j], &xinterfaces[(i + 1) * (block.ny + 1) + j],
				&yinterfaces[i * (block.ny + 1) + j], &yinterfaces[i * (block.ny + 1) + j + 1], &fluids[i * block.ny + j], block);
		}
	}
}

void Check_Order_Reduce_by_Lambda_2D(bool& order_reduce, double* convar)
{
	order_reduce = false;
	double lambda;
	lambda = Lambda(convar[0], convar[1] / convar[0], convar[2] / convar[0], convar[3]);
	//if lambda <0, then reduce to the first order
	if (lambda <= 0.0 || (lambda == lambda) == false)
	{
		order_reduce = true;
	}
}

double Calculate_alpha_k(double* prim_left, double* prim_right)
{
	double A, Ma_left_normal, Ma_right_normal;
	double Ma_left_tangent, Ma_right_tangent;

	Ma_left_normal  = prim_left[1] / sqrt(Gamma * prim_left[3] / prim_left[0]);
	Ma_left_tangent = prim_left[2] / sqrt(Gamma * prim_left[3] / prim_left[0]);
	Ma_right_normal = prim_right[1] / sqrt(Gamma * prim_right[3] / prim_right[0]);
	Ma_right_tangent = prim_right[2] / sqrt(Gamma * prim_right[3] / prim_right[0]);

	A = abs(prim_left[3] - prim_right[3]) / prim_left[3] + abs(prim_left[3] - prim_right[3]) / prim_right[3];
	A = A + pow(Ma_left_normal - Ma_right_normal, 2) + pow(Ma_left_tangent - Ma_right_tangent, 2);
	return A;
}

void Calculate_alpha(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d& fluids)
{
	double prim_left_left[4], prim_left_right[4], prim_right_left[4], prim_right_right[4];
	double prim_down_left[4], prim_down_right[4], prim_up_left[4], prim_up_right[4];
	double alpha_left, alpha_right, alpha_down, alpha_up;
	double alpha[10][4];
	for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
	{
		// here two gaussian points are used by default
		Convar_to_primvar_2D(prim_left_left, left.gauss[num_gauss].left.convar);
		Convar_to_primvar_2D(prim_left_right, left.gauss[num_gauss].right.convar);
		Convar_to_primvar_2D(prim_right_left, right.gauss[num_gauss].left.convar);
		Convar_to_primvar_2D(prim_right_right, right.gauss[num_gauss].right.convar);
		// here
		Convar_to_primvar_2D(prim_down_left, down.gauss[num_gauss].left.convar);
		Convar_to_primvar_2D(prim_down_right, down.gauss[num_gauss].right.convar);
		Convar_to_primvar_2D(prim_up_left, up.gauss[num_gauss].left.convar);
		Convar_to_primvar_2D(prim_up_right, up.gauss[num_gauss].right.convar);
		double down_left[4], down_right[4], up_left[4], up_right[4];
		XchangetoY(down_left, prim_down_left); XchangetoY(down_right, prim_down_right); 
		XchangetoY(up_left, prim_up_left); XchangetoY(up_right, prim_up_right);

		alpha[num_gauss][0] = Calculate_alpha_k(prim_left_left, prim_left_right);
		alpha[num_gauss][1] = Calculate_alpha_k(prim_right_left, prim_right_right);
		alpha[num_gauss][2] = Calculate_alpha_k(down_left, down_right);
		alpha[num_gauss][3] = Calculate_alpha_k(up_left, up_right);
	}
	if (gausspoint == 2)
	{
		alpha_left = 0.5 * (alpha[0][0] + alpha[1][0]);
		alpha_right = 0.5 * (alpha[0][1] + alpha[1][1]);
		alpha_down = 0.5 * (alpha[0][2] + alpha[1][2]);
		alpha_up = 0.5 * (alpha[0][3] + alpha[1][3]);
	}
	if (gausspoint == 3)
	{
		alpha_left = (alpha[0][0] + alpha[1][0] + alpha[2][0]) / 3.0;
		alpha_right = (alpha[0][1] + alpha[1][1] + alpha[2][1]) / 3.0;
		alpha_down = (alpha[0][2] + alpha[1][2] + alpha[2][2]) / 3.0;
		alpha_up = (alpha[0][3] + alpha[1][3] + alpha[2][3]) / 3.0;
	}
	if (gausspoint == 4)
	{
		alpha_left = 0.25 * (alpha[0][0] + alpha[1][0] + alpha[2][0] + alpha[3][0]);
		alpha_right = 0.25 * (alpha[0][1] + alpha[1][1] + alpha[2][1] + alpha[3][1]);
		alpha_down = 0.25 * (alpha[0][2] + alpha[1][2] + alpha[2][2] + alpha[3][2]);
		alpha_up = 0.25 * (alpha[0][3] + alpha[1][3] + alpha[2][3] + alpha[3][3]);
	}
	alpha_left = alpha_left < 2.0 ? 0.0 : alpha_left;
	alpha_right = alpha_right < 2.0 ? 0.0 : alpha_right;
	alpha_down = alpha_down < 2.0 ? 0.0 : alpha_down;
	alpha_up = alpha_up < 2.0 ? 0.0 : alpha_up;

	fluids.alpha_x = alpha_left + alpha_right;
	fluids.alpha_y = alpha_down + alpha_up;
}

void Update_alpha(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodex + block.ghost; j++)
		{
			Calculate_alpha(xinterfaces[i * (block.ny + 1) + j], xinterfaces[(i + 1) * (block.ny + 1) + j],
				yinterfaces[i * (block.ny + 1) + j], yinterfaces[i * (block.ny + 1) + j + 1], fluids[i * block.ny + j]);
		}
	}
}

// cell left & right side normal reconstruction
void First_order_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block)
{
	first_order(left.line.right, right.line.left,
		left.normal, right.normal, fluids[0].convar);

	first_order(down.line.right, up.line.left,
		down.normal, up.normal, fluids[0].convar);
}

void first_order(Point2d& left, Point2d& right, double* normal_l, double* normal_r, double* w)
{
	double splus[4], sminus[4];

	double w_l[4];
	Global_to_Local(w_l, w, normal_l);
	Copy_Array(left.convar, w_l, 4);
	Array_zero(left.der1x, 4);

	double w_r[4];
	Global_to_Local(w_r, w, normal_r);
	Copy_Array(right.convar, w_r, 4);
	Array_zero(right.der1x, 4);

}

void WENO5_AO_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block)	
{
	if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
		{
			WENO5_AO(left.line.right, right.line.left, fluids[-2 * block.ny].convar, fluids[-block.ny].convar, fluids[0].convar, fluids[block.ny].convar, fluids[2 * block.ny].convar, fluids[0].dx);
		}


	if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
		{
			double wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4];
			YchangetoX(wn1tmp, fluids[-1].convar); YchangetoX(wtmp, fluids[0].convar); YchangetoX(wp1tmp, fluids[1].convar);
			YchangetoX(wn2tmp, fluids[-2].convar); YchangetoX(wp2tmp, fluids[2].convar);

			WENO5_AO(down.line.right, up.line.left, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, fluids[0].dy);
		}
}

void WENO5_AO(Point2d& left, Point2d& right, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double h)
{
	//we denote that   |left...cell-center...right|
	double ren2[4], ren1[4], re0[4], rep1[4], rep2[4];
	double var[4], der1[4], der2[4];

	double base_left[4];
	double base_right[4];
	double wn1_primvar[4], w_primvar[4], wp1_primvar[4];
	Convar_to_primvar_2D(wn1_primvar, wn1);
	Convar_to_primvar_2D(w_primvar, w);
	Convar_to_primvar_2D(wp1_primvar, wp1);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = 0.5 * (wn1_primvar[i] + w_primvar[i]);
		base_right[i] = 0.5 * (wp1_primvar[i] + w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
		}
	}
	else
	{
		Convar_to_char(ren2, base_left, wn2);
		Convar_to_char(ren1, base_left, wn1);
		Convar_to_char(re0, base_left, w);
		Convar_to_char(rep1, base_left, wp1);
		Convar_to_char(rep2, base_left, wp2);
	}

	for (int i = 0; i < 4; i++)
	{
		weno_5th_ao_left(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h);

	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			left.convar[i] = var[i];
			left.der1x[i] = der1[i];

		}
	}
	else
	{
		Char_to_convar(left.convar, base_left, var);
		Char_to_convar(left.der1x, base_left, der1);

	}

	// cell right
	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
		}
	}
	else
	{
		Convar_to_char(ren2, base_right, wn2);
		Convar_to_char(ren1, base_right, wn1);
		Convar_to_char(re0, base_right, w);
		Convar_to_char(rep1, base_right, wp1);
		Convar_to_char(rep2, base_right, wp2);
	}

	for (int i = 0; i < 4; i++)
	{
		weno_5th_ao_right(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			right.convar[i] = var[i];
			right.der1x[i] = der1[i];
		}
	}
	else
	{
		Char_to_convar(right.convar, base_right, var);
		Char_to_convar(right.der1x, base_right, der1);

	}

	
	Check_Order_Reduce_by_Lambda_2D(right.is_reduce_order, right.convar);
	Check_Order_Reduce_by_Lambda_2D(left.is_reduce_order, left.convar);

	if (left.is_reduce_order == true || right.is_reduce_order == true)
	{
		//if (is_reduce_order_warning == true)
			//cout << " WENO5-cell-splitting order reduce" << endl;
		for (int m = 0; m < 4; m++)
		{
			right.convar[m] = w[m];
			left.convar[m] = w[m];
			right.der1x[m] = 0.0;
			left.der1x[m] = 0.0;
		}
	}
	
}

void WENO5_AO_with_df_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block)
{
	if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
	{
		double alpha[5];
		alpha[0] = fluids[-2 * block.ny].alpha_x;
		alpha[1] = fluids[-block.ny].alpha_x;
		alpha[2] = fluids[0].alpha_x;
		alpha[3] = fluids[block.ny].alpha_x;
		alpha[4] = fluids[2 * block.ny].alpha_x;
		WENO5_AO_with_df(left.line.right, right.line.left, alpha, fluids[-2 * block.ny].convar, fluids[-block.ny].convar, fluids[0].convar, fluids[block.ny].convar, fluids[2 * block.ny].convar, fluids[0].dx);
	}


	if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
	{
		double wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4];
		YchangetoX(wn1tmp, fluids[-1].convar); YchangetoX(wtmp, fluids[0].convar); YchangetoX(wp1tmp, fluids[1].convar);
		YchangetoX(wn2tmp, fluids[-2].convar); YchangetoX(wp2tmp, fluids[2].convar);

		double alpha[5];
		alpha[0] = fluids[-2].alpha_y;
		alpha[1] = fluids[-1].alpha_y;
		alpha[2] = fluids[0].alpha_y;
		alpha[3] = fluids[1].alpha_y;
		alpha[4] = fluids[2].alpha_y;
		WENO5_AO_with_df(down.line.right, up.line.left, alpha, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, fluids[0].dy);
	}
}

void WENO5_AO_with_df(Point2d& left, Point2d& right, double* alpha, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double h)
{
	double df[4];
	double sum_df[4];
	sum_df[0] = alpha[0] + alpha[2];
	sum_df[1] = alpha[1] + alpha[3];
	sum_df[2] = alpha[2] + alpha[4];
	sum_df[3] = alpha[0] + alpha[2] + alpha[4];
	for (int i = 0; i < 4; i++)
	{
		df[i] = sum_df[i] < 2.0 ? 1.0 : 2.0 / sum_df[i];
	}

	//we denote that   |left...cell-center...right|
	double ren2[4], ren1[4], re0[4], rep1[4], rep2[4];
	double var[4], der1[4], der2[4];

	double base_left[4];
	double base_right[4];
	double wn1_primvar[4], w_primvar[4], wp1_primvar[4];
	Convar_to_primvar_2D(wn1_primvar, wn1);
	Convar_to_primvar_2D(w_primvar, w);
	Convar_to_primvar_2D(wp1_primvar, wp1);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = 0.5 * (wn1_primvar[i] + w_primvar[i]);
		base_right[i] = 0.5 * (wp1_primvar[i] + w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
		}
	}
	else
	{
		Convar_to_char(ren2, base_left, wn2);
		Convar_to_char(ren1, base_left, wn1);
		Convar_to_char(re0, base_left, w);
		Convar_to_char(rep1, base_left, wp1);
		Convar_to_char(rep2, base_left, wp2);
	}

	for (int i = 0; i < 4; i++)
	{
		weno_5th_ao_with_df_left(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], df, h);

	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			left.convar[i] = var[i];
			left.der1x[i] = der1[i];

		}
	}
	else
	{
		Char_to_convar(left.convar, base_left, var);
		Char_to_convar(left.der1x, base_left, der1);

	}

	// cell right
	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
		}
	}
	else
	{
		Convar_to_char(ren2, base_right, wn2);
		Convar_to_char(ren1, base_right, wn1);
		Convar_to_char(re0, base_right, w);
		Convar_to_char(rep1, base_right, wp1);
		Convar_to_char(rep2, base_right, wp2);
	}

	for (int i = 0; i < 4; i++)
	{
		weno_5th_ao_with_df_right(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], df, h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			right.convar[i] = var[i];
			right.der1x[i] = der1[i];
		}
	}
	else
	{
		Char_to_convar(right.convar, base_right, var);
		Char_to_convar(right.der1x, base_right, der1);

	}

	
	Check_Order_Reduce_by_Lambda_2D(right.is_reduce_order, right.convar);
	Check_Order_Reduce_by_Lambda_2D(left.is_reduce_order, left.convar);

	if (left.is_reduce_order == true || right.is_reduce_order == true)
	{
		for (int m = 0; m < 4; m++)
		{
			right.convar[m] = w[m];
			left.convar[m] = w[m];
			right.der1x[m] = 0.0;
			left.der1x[m] = 0.0;
		}
	}
	
}

void WENO7_AO_with_df_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block)
{
	if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
	{
		double alpha[7];
		alpha[0] = fluids[-3 * block.ny].alpha_x;
		alpha[1] = fluids[-2 * block.ny].alpha_x;
		alpha[2] = fluids[-block.ny].alpha_x;
		alpha[3] = fluids[0].alpha_x;
		alpha[4] = fluids[block.ny].alpha_x;
		alpha[5] = fluids[2 * block.ny].alpha_x;
		alpha[6] = fluids[3 * block.ny].alpha_x;
		WENO7_AO_with_df(left.line.right, right.line.left, alpha, fluids[-3 * block.ny].convar, fluids[-2 * block.ny].convar, fluids[-block.ny].convar, fluids[0].convar, fluids[block.ny].convar, fluids[2 * block.ny].convar, fluids[3 * block.ny].convar, fluids[0].dx);
	}


	if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
	{
		double wn3tmp[4], wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4], wp3tmp[4];
		YchangetoX(wn1tmp, fluids[-1].convar); YchangetoX(wtmp, fluids[0].convar); YchangetoX(wp1tmp, fluids[1].convar);
		YchangetoX(wn2tmp, fluids[-2].convar); YchangetoX(wp2tmp, fluids[2].convar);
		YchangetoX(wn3tmp, fluids[-3].convar); YchangetoX(wp3tmp, fluids[3].convar);

		double alpha[7];
		alpha[0] = fluids[-3].alpha_y;
		alpha[1] = fluids[-2].alpha_y;
		alpha[2] = fluids[-1].alpha_y;
		alpha[3] = fluids[0].alpha_y;
		alpha[4] = fluids[1].alpha_y;
		alpha[5] = fluids[2].alpha_y;
		alpha[6] = fluids[3].alpha_y;
		WENO7_AO_with_df(down.line.right, up.line.left, alpha, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, fluids[0].dy);
	}
}

void WENO7_AO_with_df(Point2d& left, Point2d& right, double* alpha, double* wn3, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double* wp3, double h)
{
	double df[5];
	double sum_df[5];
	sum_df[0] = alpha[1] + alpha[3];
	sum_df[1] = alpha[2] + alpha[4];
	sum_df[2] = alpha[3] + alpha[5];
	sum_df[3] = alpha[1] + alpha[3] + alpha[5];
	sum_df[4] = alpha[0] + alpha[2] + alpha[4] + alpha[6];
	
	for (int i = 0; i < 5; i++)
	{
		df[i] = sum_df[i] < 2.0 ? 1.0 : 2.0 / sum_df[i];
	}
	//we denote that   |left...cell-center...right|
	double ren3[4], ren2[4], ren1[4], re0[4], rep1[4], rep2[4], rep3[4];
	double var[4], der1[4], der2[4];

	double base_left[4];
	double base_right[4];
	double wn1_primvar[4], w_primvar[4], wp1_primvar[4];
	Convar_to_primvar_2D(wn1_primvar, wn1);
	Convar_to_primvar_2D(w_primvar, w);
	Convar_to_primvar_2D(wp1_primvar, wp1);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = 0.5 * (wn1_primvar[i] + w_primvar[i]);
		base_right[i] = 0.5 * (wp1_primvar[i] + w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren3[i] = wn3[i];
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
			rep3[i] = wp3[i];
		}
	}
	else
	{
		Convar_to_char(ren3, base_left, wn3);
		Convar_to_char(ren2, base_left, wn2);
		Convar_to_char(ren1, base_left, wn1);
		Convar_to_char(re0, base_left, w);
		Convar_to_char(rep1, base_left, wp1);
		Convar_to_char(rep2, base_left, wp2);
		Convar_to_char(rep3, base_left, wp3);
	}

	for (int i = 0; i < 4; i++)
	{
		weno_7th_ao_with_df_left(var[i], der1[i], der2[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], df, h);

	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			left.convar[i] = var[i];
			left.der1x[i] = der1[i];
		}
	}
	else
	{
		Char_to_convar(left.convar, base_left, var);
		Char_to_convar(left.der1x, base_left, der1);
	}

	// cell right
	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren3[i] = wn3[i];
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
			rep3[i] = wp3[i];
		}
	}
	else
	{
		Convar_to_char(ren3, base_right, wn3);
		Convar_to_char(ren2, base_right, wn2);
		Convar_to_char(ren1, base_right, wn1);
		Convar_to_char(re0, base_right, w);
		Convar_to_char(rep1, base_right, wp1);
		Convar_to_char(rep2, base_right, wp2);
		Convar_to_char(rep3, base_right, wp3);
	}

	for (int i = 0; i < 4; i++)
	{
		weno_7th_ao_with_df_right(var[i], der1[i], der2[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], df, h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			right.convar[i] = var[i];
			right.der1x[i] = der1[i];
		}
	}
	else
	{
		Char_to_convar(right.convar, base_right, var);
		Char_to_convar(right.der1x, base_right, der1);
	}
	
	Check_Order_Reduce_by_Lambda_2D(right.is_reduce_order, right.convar);
	Check_Order_Reduce_by_Lambda_2D(left.is_reduce_order, left.convar);

	if (left.is_reduce_order == true || right.is_reduce_order == true)
	{
		for (int m = 0; m < 4; m++)
		{
			right.convar[m] = w[m];
			left.convar[m] = w[m];
			right.der1x[m] = 0.0;
			left.der1x[m] = 0.0;
		}
	}
}

void WENO9_AO_with_df_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block)
{
	if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
	{
		double alpha[9];
		alpha[0] = fluids[-4 * block.ny].alpha_x;
		alpha[1] = fluids[-3 * block.ny].alpha_x;
		alpha[2] = fluids[-2 * block.ny].alpha_x;
		alpha[3] = fluids[-block.ny].alpha_x;
		alpha[4] = fluids[0].alpha_x;
		alpha[5] = fluids[block.ny].alpha_x;
		alpha[6] = fluids[2 * block.ny].alpha_x;
		alpha[7] = fluids[3 * block.ny].alpha_x;
		alpha[8] = fluids[4 * block.ny].alpha_x;
		WENO9_AO_with_df(left.line.right, right.line.left, alpha, fluids[-4 * block.ny].convar, fluids[-3 * block.ny].convar, 
		fluids[-2 * block.ny].convar, fluids[-block.ny].convar, fluids[0].convar, fluids[block.ny].convar, fluids[2 * block.ny].convar, fluids[3 * block.ny].convar,fluids[4 * block.ny].convar, fluids[0].dx);
	}

	if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
	{
		double wn4tmp[4], wn3tmp[4], wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4], wp3tmp[4], wp4tmp[4];
		YchangetoX(wn1tmp, fluids[-1].convar); YchangetoX(wtmp, fluids[0].convar); YchangetoX(wp1tmp, fluids[1].convar);
		YchangetoX(wn2tmp, fluids[-2].convar); YchangetoX(wp2tmp, fluids[2].convar);
		YchangetoX(wn3tmp, fluids[-3].convar); YchangetoX(wp3tmp, fluids[3].convar);
		YchangetoX(wn4tmp, fluids[-4].convar); YchangetoX(wp4tmp, fluids[4].convar);

		double alpha[9];
		alpha[0] = fluids[-4].alpha_y;
		alpha[1] = fluids[-3].alpha_y;
		alpha[2] = fluids[-2].alpha_y;
		alpha[3] = fluids[-1].alpha_y;
		alpha[4] = fluids[0].alpha_y;
		alpha[5] = fluids[1].alpha_y;
		alpha[6] = fluids[2].alpha_y;
		alpha[7] = fluids[3].alpha_y;
		alpha[8] = fluids[4].alpha_y;
		WENO9_AO_with_df(down.line.right, up.line.left, alpha, wn4tmp, wn3tmp, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, wp3tmp, wp4tmp, fluids[0].dy);
	}
}

void WENO9_AO_with_df(Point2d& left, Point2d& right, double* alpha, double* wn4, double* wn3, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double* wp3, double* wp4, double h)
{
	double df[6];
	double sum_df[6];
	sum_df[0] = alpha[2] + alpha[4];
	sum_df[1] = alpha[3] + alpha[5];
	sum_df[2] = alpha[4] + alpha[6];
	sum_df[3] = alpha[2] + alpha[4] + alpha[6];
	sum_df[4] = alpha[1] + alpha[3] + alpha[5] + alpha[7];
	sum_df[5] = alpha[0] + alpha[2] + alpha[4] + alpha[6] + alpha[8];
	
	for (int i = 0; i < 6; i++)
	{
		df[i] = sum_df[i] < 2.0 ? 1.0 : 2.0 / sum_df[i];
	}
	//we denote that   |left...cell-center...right|
	double ren4[4], ren3[4], ren2[4], ren1[4], re0[4], rep1[4], rep2[4], rep3[4], rep4[4];
	double var[4], der1[4], der2[4];

	double base_left[4];
	double base_right[4];
	double wn1_primvar[4], w_primvar[4], wp1_primvar[4];
	Convar_to_primvar_2D(wn1_primvar, wn1);
	Convar_to_primvar_2D(w_primvar, w);
	Convar_to_primvar_2D(wp1_primvar, wp1);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = 0.5 * (wn1_primvar[i] + w_primvar[i]);
		base_right[i] = 0.5 * (wp1_primvar[i] + w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren4[i] = wn4[i];
			ren3[i] = wn3[i];
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
			rep3[i] = wp3[i];
			rep4[i] = wp4[i];
		}
	}
	else
	{
		Convar_to_char(ren4, base_left, wn4);
		Convar_to_char(ren3, base_left, wn3);
		Convar_to_char(ren2, base_left, wn2);
		Convar_to_char(ren1, base_left, wn1);
		Convar_to_char(re0, base_left, w);
		Convar_to_char(rep1, base_left, wp1);
		Convar_to_char(rep2, base_left, wp2);
		Convar_to_char(rep3, base_left, wp3);
		Convar_to_char(rep4, base_left, wp4);
	}

	for (int i = 0; i < 4; i++)
	{
		weno_9th_ao_with_df_left(var[i], der1[i], der2[i], ren4[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], rep4[i], df, h);

	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			left.convar[i] = var[i];
			left.der1x[i] = der1[i];
		}
	}
	else
	{
		Char_to_convar(left.convar, base_left, var);
		Char_to_convar(left.der1x, base_left, der1);
	}

	// cell right
	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren4[i] = wn4[i];
			ren3[i] = wn3[i];
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
			rep3[i] = wp3[i];
			rep4[i] = wp4[i];
		}
	}
	else
	{
		Convar_to_char(ren4, base_right, wn4);
		Convar_to_char(ren3, base_right, wn3);
		Convar_to_char(ren2, base_right, wn2);
		Convar_to_char(ren1, base_right, wn1);
		Convar_to_char(re0, base_right, w);
		Convar_to_char(rep1, base_right, wp1);
		Convar_to_char(rep2, base_right, wp2);
		Convar_to_char(rep3, base_right, wp3);
		Convar_to_char(rep4, base_right, wp4);
	}

	for (int i = 0; i < 4; i++)
	{
		weno_9th_ao_with_df_right(var[i], der1[i], der2[i], ren4[i], ren3[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], rep4[i], df, h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			right.convar[i] = var[i];
			right.der1x[i] = der1[i];
		}
	}
	else
	{
		Char_to_convar(right.convar, base_right, var);
		Char_to_convar(right.der1x, base_right, der1);
	}
	
	Check_Order_Reduce_by_Lambda_2D(right.is_reduce_order, right.convar);
	Check_Order_Reduce_by_Lambda_2D(left.is_reduce_order, left.convar);

	if (left.is_reduce_order == true || right.is_reduce_order == true)
	{
		for (int m = 0; m < 4; m++)
		{
			right.convar[m] = w[m];
			left.convar[m] = w[m];
			right.der1x[m] = 0.0;
			left.der1x[m] = 0.0;
		}
	}
}

// cell left & right side tangent reconstruction
void First_order_tangent(Interface2d* left, Interface2d* right, Interface2d* down, Interface2d* up, Fluid2d* fluids, Block2d block)
{
	for (int num_gauss = 0; num_gauss < gausspoint; ++num_gauss)
	{
		first_order_tangent(left[0].gauss[num_gauss].right, left[0].line.right);
		first_order_tangent(right[0].gauss[num_gauss].left, right[0].line.left);
		first_order_tangent(down[0].gauss[num_gauss].right, down[0].line.right);
		first_order_tangent(up[0].gauss[num_gauss].left, up[0].line.left);
	}
}

void first_order_tangent(Point2d& gauss, Point2d& w0)
{
	//tangential
	Copy_Array(gauss.convar, w0.convar, 4);
	Copy_Array(gauss.der1x, w0.der1x, 4);
	Array_zero(gauss.der1y, 4);
}

void WENO5_AO_tangent(Interface2d* left, Interface2d* right, Interface2d* down, Interface2d* up, Fluid2d* fluids, Block2d block)
{
	// along x direction tangitial recontruction,
	WENO5_AO_tangential(right[0].gauss,
			right[-2].line, right[-1].line, right[0].line, right[1].line, right[2].line, right[0].length);

	//since we already do the coordinate transform, along y, no transform needed.
	WENO5_AO_tangential(up[0].gauss, up[2 * (block.ny + 1)].line, up[block.ny + 1].line,
			up[0].line, up[-(block.ny + 1)].line, up[-2 * (block.ny + 1)].line, up[0].length);
}

void WENO5_AO_tangential(Recon2d* re, Recon2d& wn2, Recon2d& wn1, Recon2d& w0, Recon2d& wp1, Recon2d& wp2, double h)
{
	//lets first reconstruction the left value

	double ren2[4], ren1[4], re0[4], rep1[4], rep2[4];
	double base_left[4];
	double base_right[4];
	double wn1_primvar[4], w_primvar[4], wp1_primvar[4];
	Convar_to_primvar_2D(w_primvar, w0.left.convar);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = (w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		double tmp[2];
		for (int i = 0; i < 4; i++)
		{
			if (gausspoint == 2)
			{
				weno_5th_ao_2gauss(re[0].left.convar[i], re[0].left.der1y[i], tmp[0],
					re[1].left.convar[i], re[1].left.der1y[i], tmp[1],
					wn2.left.convar[i], wn1.left.convar[i], w0.left.convar[i], wp1.left.convar[i], wp2.left.convar[i], h, 2);
			}
		}
	}
	else
	{
		Convar_to_char(ren2, base_left, wn2.left.convar);
		Convar_to_char(ren1, base_left, wn1.left.convar);
		Convar_to_char(re0, base_left, w0.left.convar);
		Convar_to_char(rep1, base_left, wp1.left.convar);
		Convar_to_char(rep2, base_left, wp2.left.convar);
		if (gausspoint == 2)
		{
			double var[2][4], der1[2][4], der2[2][4];
			for (int i = 0; i < 4; i++)
			{
				weno_5th_ao_2gauss(var[0][i], der1[0][i], der2[0][i],
					var[1][i], der1[1][i], der2[1][i],
					ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h, 2);
			}
			for (int igauss = 0; igauss < gausspoint; igauss++)
			{
				Char_to_convar(re[igauss].left.convar, base_left, var[igauss]);
				Char_to_convar(re[igauss].left.der1y, base_left, der1[igauss]);
			}
		}

	}

	for (int i = 0; i < 4; i++)
	{
		double tmp[4];
		if (gausspoint == 2)
		{
			weno_5th_ao_2gauss(re[0].left.der1x[i], tmp[0], tmp[1],
				re[1].left.der1x[i], tmp[2], tmp[3],
				wn2.left.der1x[i], wn1.left.der1x[i], w0.left.der1x[i], wp1.left.der1x[i], wp2.left.der1x[i], h, 1);
		}
	}

	//then let's construct the right part....
	Convar_to_primvar_2D(w_primvar, w0.right.convar);
	for (int i = 0; i < 4; i++)
	{
		base_right[i] = (w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		double tmp[2];
		for (int i = 0; i < 4; i++)
		{
			if (gausspoint == 2)
			{
				weno_5th_ao_2gauss(re[0].right.convar[i], re[0].right.der1y[i], tmp[0],
					re[1].right.convar[i], re[1].right.der1y[i], tmp[1],
					wn2.right.convar[i], wn1.right.convar[i], w0.right.convar[i], wp1.right.convar[i], wp2.right.convar[i], h, 2);
			}
		}
	}
	else
	{
		Convar_to_char(ren2, base_right, wn2.right.convar);
		Convar_to_char(ren1, base_right, wn1.right.convar);
		Convar_to_char(re0, base_right, w0.right.convar);
		Convar_to_char(rep1, base_right, wp1.right.convar);
		Convar_to_char(rep2, base_right, wp2.right.convar);
		if (gausspoint == 2)
		{
			double var[2][4], der1[2][4], der2[2][4];
			for (int i = 0; i < 4; i++)
			{
				weno_5th_ao_2gauss(var[0][i], der1[0][i], der2[0][i],
					var[1][i], der1[1][i], der2[1][i],
					ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h, 2);
			}
			for (int igauss = 0; igauss < gausspoint; igauss++)
			{
				Char_to_convar(re[igauss].right.convar, base_right, var[igauss]);
				Char_to_convar(re[igauss].right.der1y, base_right, der1[igauss]);
			}
		}
	}

	for (int i = 0; i < 4; i++)
	{
		double tmp[4];
		if (gausspoint == 2)
		{
			weno_5th_ao_2gauss(re[0].right.der1x[i], tmp[0], tmp[1],
				re[1].right.der1x[i], tmp[2], tmp[3],
				wn2.right.der1x[i], wn1.right.der1x[i], w0.right.der1x[i], wp1.right.der1x[i], wp2.right.der1x[i], h, 1);
		}
	}	
	for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
	{
		Check_Order_Reduce_by_Lambda_2D(re[num_gauss].left.is_reduce_order, re[num_gauss].left.convar);
		if (re[num_gauss].left.is_reduce_order == true)
		{
			//cout << " WENO5-cell-multi-left order reduce" << endl;
			for (int var = 0; var < 4; var++)
			{
				re[num_gauss].left.convar[var] = w0.left.convar[var];
				re[num_gauss].left.der1x[var] = w0.left.der1x[var];
				re[num_gauss].left.der1y[var] = 0.0;
			}
		}
	}

	for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
	{
		Check_Order_Reduce_by_Lambda_2D(re[num_gauss].right.is_reduce_order, re[num_gauss].right.convar);
		if (re[num_gauss].right.is_reduce_order == true)
		{
			//cout << " WENO5-cell-multi-right order reduce" << endl;
			for (int var = 0; var < 4; var++)
			{
				re[num_gauss].right.convar[var] = w0.right.convar[var];
				re[num_gauss].right.der1x[var] = w0.right.der1x[var];
				re[num_gauss].right.der1y[var] = 0.0;
			}
		}
	}
}

void weno_5th_ao_2gauss(double& g1, double& g1x, double& g1xx, double& g2, double& g2x, double& g2xx, double wn2, double wn1, double w0, double wp1, double wp2, double h, int order)
{
	//the parameter order constrols up to which order you want construct
	// order from 0, 1, 2
	if (order > 2 || order < 0)
	{
		cout << "invalid order input for the function " << __FUNCTION__ << endl;
		exit(0);
	}
	double dhi = 0.85;
	double dlo = 0.85;
	//-- - parameter of WENO-- -
	double beta[4], d[4], ww[4], alpha[4];
	double epsilonW = 1e-6;
	//-- - intermediate parameter-- -
	double p[4], px[4], pxx[4], tempvar;
	double sum_alpha;

	//three small stencil
	d[0] = (1 - dhi) * (1 - dlo) * 0.5;
	d[1] = (1 - dhi) * dlo;
	d[2] = (1 - dhi) * (1 - dlo) * 0.5;
	//one big stencil
	d[3] = dhi;

	//cout << "here" << endl;
	beta[0] = 13.0 / 12.0 * pow((wn2 - 2 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4 * wn1 + 3 * w0), 2);
	beta[1] = 13.0 / 12.0 * pow((wn1 - 2 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
	beta[2] = 13.0 / 12.0 * pow((w0 - 2 * wp1 + wp2), 2) + 0.25 * pow((3 * w0 - 4 * wp1 + wp2), 2);

	beta[3] = (1.0 / 5040.0) * (231153 * w0 * w0 + 104963 * wn1 * wn1 + 6908 * wn2 * wn2 -
		38947 * wn2 * wp1 + 104963 * wp1 * wp1 +
		wn1 * (-51001 * wn2 + 179098 * wp1 - 38947 * wp2) -
		3 * w0 * (99692 * wn1 - 22641 * wn2 + 99692 * wp1 - 22641 * wp2) +
		8209 * wn2 * wp2 - 51001 * wp1 * wp2 + 6908 * wp2 * wp2);

	double tau5 = 1.0 / 3.0 * (std::abs(beta[3] - beta[0]) + std::abs(beta[3] - beta[1]) + std::abs(beta[3] - beta[2]));


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


	//cout <<"substencil "<< h0 << " " << h1 << " " << h2 << endl;
	////-- - combination-- -

	double final_weight[4];
	final_weight[3] = ww[3] / d[3];
	for (int k = 0; k < 3; k++)
	{
		final_weight[k] = ww[k] - ww[3] / d[3] * d[k];
	}

	g1 = 0; g1x = 0; g1xx = 0;
	g2 = 0; g2x = 0; g2xx = 0;

	double sqrt3 = 1.732050807568877294;
	//-- - candidate polynomial-- for gauss 1 point
	p[0] = w0 - (sqrt3 * w0) / 4 + (4 * wn1 - wn2) / (4 * sqrt3);
	p[1] = w0 + (wn1 - wp1) / (4 * sqrt3);
	p[2] = (1.0 / 12.0) * (3 * (4 + sqrt3) * w0 + sqrt3 * (-4 * wp1 + wp2));
	p[3] = (4314 * w0 + (4 + 500 * sqrt3) * wn1 - wn2 - 70 * sqrt3 * wn2 +
		4 * wp1 - 500 * sqrt3 * wp1 - wp2 + 70 * sqrt3 * wp2) / 4320;

	for (int k = 0; k < 4; k++)
	{
		g1 += final_weight[k] * p[k];
	}

	if (order > 0)
	{
		px[0] = -((-9 + sqrt3) * w0 - 2 * (-6 + sqrt3) * wn1 + (-3 + sqrt3) * wn2) / (6 * h);
		px[1] = -(-2 * sqrt3 * w0 + (3 + sqrt3) * wn1 + (-3 + sqrt3) * wp1) / (6 * h);
		px[2] = -((9 + sqrt3) * w0 - 2 * (6 + sqrt3) * wp1 + (3 + sqrt3) * wp2) / (6 * h);
		px[3] = (48 * sqrt3 * w0 - 72 * wn1 - 26 * sqrt3 * wn1 + 9 * wn2 +
			2 * sqrt3 * wn2 + 72 * wp1 - 26 * sqrt3 * wp1 - 9 * wp2 +
			2 * sqrt3 * wp2) / (108 * h);

		for (int k = 0; k < 4; k++)
		{
			g1x += final_weight[k] * px[k];
		}
		if (order == 2)
		{
			pxx[0] = (w0 - 2 * wn1 + wn2) / h / h;
			pxx[1] = (-2 * w0 + wn1 + wp1) / h / h;
			pxx[2] = (w0 - 2 * wp1 + wp2) / h / h;
			pxx[3] = -((30 * w0 + 2 * (-8 + sqrt3) * wn1 + wn2 - sqrt3 * wn2 - 16 * wp1 -
				2 * sqrt3 * wp1 + wp2 + sqrt3 * wp2) / (12 * h * h));

			for (int k = 0; k < 4; k++)
			{
				g1xx += final_weight[k] * pxx[k];
			}
		}
	}


	//-- - candidate polynomial-- for gauss 2 point
	p[0] = (1.0 / 12.0) * (3 * (4 + sqrt3) * w0 + sqrt3 * (-4 * wn1 + wn2));
	p[1] = w0 + (-wn1 + wp1) / (4 * sqrt3);
	p[2] = w0 - (sqrt3 * w0) / 4 + (4 * wp1 - wp2) / (4 * sqrt3);
	p[3] = (4314 * w0 + (4 - 500 * sqrt3) * wn1 - wn2 + 70 * sqrt3 * wn2 +
		4 * wp1 + 500 * sqrt3 * wp1 - wp2 - 70 * sqrt3 * wp2) / 4320;
	for (int k = 0; k < 4; k++)
	{
		g2 += final_weight[k] * p[k];
	}
	if (order > 0)
	{
		px[0] = ((9 + sqrt3) * w0 - 2 * (6 + sqrt3) * wn1 + (3 + sqrt3) * wn2) /
			(6 * h);
		px[1] = (-2 * sqrt3 * w0 + (-3 + sqrt3) * wn1 + (3 + sqrt3) * wp1) / (6 * h);
		px[2] = ((-9 + sqrt3) * w0 - 2 * (-6 + sqrt3) * wp1 + (-3 + sqrt3) * wp2) / (6 * h);
		px[3] = -((48 * sqrt3 * w0 + 72 * wn1 - 26 * sqrt3 * wn1 - 9 * wn2 +
			2 * sqrt3 * wn2 - 72 * wp1 - 26 * sqrt3 * wp1 + 9 * wp2 +
			2 * sqrt3 * wp2) / (108 * h));
		for (int k = 0; k < 4; k++)
		{
			g2x += final_weight[k] * px[k];
		}
		if (order == 2)
		{
			pxx[0] = (w0 - 2 * wn1 + wn2) / h / h;
			pxx[1] = (-2 * w0 + wn1 + wp1) / h / h;
			pxx[2] = (w0 - 2 * wp1 + wp2) / h / h;
			pxx[3] = -((30 * w0 - 2 * (8 + sqrt3) * wn1 + wn2 + sqrt3 * wn2 - 16 * wp1 +
				2 * sqrt3 * wp1 + wp2 - sqrt3 * wp2) / (12 * h * h));
		}
		for (int k = 0; k < 4; k++)
		{
			g2xx += final_weight[k] * pxx[k];
		}
	}
}

void WENO5_AO_with_df_tangent(Interface2d* left, Interface2d* right, Interface2d* down, Interface2d* up, Fluid2d* fluids, Block2d block)
{
	// along x direction tangitial recontruction,
	double alpha1[5], alpha2[5];
	alpha1[0] = fluids[-2].alpha_y; 
	alpha1[1] = fluids[-1].alpha_y; alpha1[2] = fluids[0].alpha_y;
	alpha1[3] = fluids[1].alpha_y;  alpha1[4] = fluids[2].alpha_y;
	alpha2[0] = fluids[block.ny - 2].alpha_y; 
	alpha2[1] = fluids[block.ny - 1].alpha_y; alpha2[2] = fluids[block.ny].alpha_y;
	alpha2[3] = fluids[block.ny + 1].alpha_y;  alpha2[4] = fluids[block.ny + 2].alpha_y;
	weno_5th_ao_with_df_tangential(right[0].gauss,
		right[-2].line, right[-1].line, right[0].line, right[1].line, right[2].line, alpha1, alpha2, right[0].length);

	//since we already do the coordinate transform, along y, no transform needed.
	alpha1[0] = fluids[2 * block.ny].alpha_x; 
	alpha1[1] = fluids[block.ny].alpha_x; alpha1[2] = fluids[0].alpha_x;
	alpha1[3] = fluids[-block.ny].alpha_x;  alpha1[4] = fluids[-2 * block.ny].alpha_x;

	alpha2[0] = fluids[2 * block.ny + 1].alpha_x; 
	alpha2[1] = fluids[block.ny + 1].alpha_x; alpha2[2] = fluids[1].alpha_x;
	alpha2[3] = fluids[-block.ny + 1].alpha_x;  alpha2[4] = fluids[-2 * block.ny + 1].alpha_x;
	weno_5th_ao_with_df_tangential(up[0].gauss, up[2 * (block.ny + 1)].line, up[block.ny + 1].line,
		up[0].line, up[-(block.ny + 1)].line, up[-2 * (block.ny + 1)].line, alpha1, alpha2, up[0].length);
}

void weno_5th_ao_with_df_tangential(Recon2d* re, Recon2d& wn2, Recon2d& wn1, Recon2d& w0, Recon2d& wp1, Recon2d& wp2, double* alpha1, double* alpha2, double h)
{
	double df1[4], df2[4];
	double sum_df1[4], sum_df2[4];
	sum_df1[0] = alpha1[0] + alpha1[2];
	sum_df1[1] = alpha1[1] + alpha1[3];
	sum_df1[2] = alpha1[2] + alpha1[4];
	sum_df1[3] = alpha1[0] + alpha1[2] + alpha1[4];

	sum_df2[0] = alpha2[0] + alpha2[2];
	sum_df2[1] = alpha2[1] + alpha2[3];
	sum_df2[2] = alpha2[2] + alpha2[4];
	sum_df2[3] = alpha2[0] + alpha2[2] + alpha2[4];
	for (int i = 0; i < 4; i++)
	{
		df1[i] = sum_df1[i] < 2.0 ? 1.0 : 2.0 / sum_df1[i];
		df2[i] = sum_df2[i] < 2.0 ? 1.0 : 2.0 / sum_df2[i];
	}
	//lets first reconstruction the left value
	double ren2[4], ren1[4], re0[4], rep1[4], rep2[4];
	double base_left[4];
	double base_right[4];
	double wn1_primvar[4], w_primvar[4], wp1_primvar[4];
	Convar_to_primvar_2D(w_primvar, w0.left.convar);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = (w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		double tmp[2];
		for (int i = 0; i < 4; i++)
		{
			if (gausspoint == 2)
			{
				weno_5th_ao_with_df_2gauss(re[0].left.convar[i], re[0].left.der1y[i], tmp[0],
					re[1].left.convar[i], re[1].left.der1y[i], tmp[1], df1,
					wn2.left.convar[i], wn1.left.convar[i], w0.left.convar[i], wp1.left.convar[i], wp2.left.convar[i], h, 2);
			}
		}
	}
	else
	{
		Convar_to_char(ren2, base_left, wn2.left.convar);
		Convar_to_char(ren1, base_left, wn1.left.convar);
		Convar_to_char(re0, base_left, w0.left.convar);
		Convar_to_char(rep1, base_left, wp1.left.convar);
		Convar_to_char(rep2, base_left, wp2.left.convar);
		if (gausspoint == 2)
		{
			double var[2][4], der1[2][4], der2[2][4];
			for (int i = 0; i < 4; i++)
			{
				weno_5th_ao_with_df_2gauss(var[0][i], der1[0][i], der2[0][i],
					var[1][i], der1[1][i], der2[1][i], df1,
					ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h, 2);
			}
			for (int igauss = 0; igauss < gausspoint; igauss++)
			{
				Char_to_convar(re[igauss].left.convar, base_left, var[igauss]);
				Char_to_convar(re[igauss].left.der1y, base_left, der1[igauss]);
			}
		}

	}

	for (int i = 0; i < 4; i++)
	{
		double tmp[4];
		if (gausspoint == 2)
		{
			weno_5th_ao_with_df_2gauss(re[0].left.der1x[i], tmp[0], tmp[1],
				re[1].left.der1x[i], tmp[2], tmp[3], df1,
				wn2.left.der1x[i], wn1.left.der1x[i], w0.left.der1x[i], wp1.left.der1x[i], wp2.left.der1x[i], h, 1);
		}
	}

	//then let's construct the right part....
	Convar_to_primvar_2D(w_primvar, w0.right.convar);
	for (int i = 0; i < 4; i++)
	{
		base_right[i] = (w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		double tmp[2];
		for (int i = 0; i < 4; i++)
		{
			if (gausspoint == 2)
			{
				weno_5th_ao_with_df_2gauss(re[0].right.convar[i], re[0].right.der1y[i], tmp[0],
					re[1].right.convar[i], re[1].right.der1y[i], tmp[1], df2,
					wn2.right.convar[i], wn1.right.convar[i], w0.right.convar[i], wp1.right.convar[i], wp2.right.convar[i], h, 2);
			}
		}
	}
	else
	{
		Convar_to_char(ren2, base_right, wn2.right.convar);
		Convar_to_char(ren1, base_right, wn1.right.convar);
		Convar_to_char(re0, base_right, w0.right.convar);
		Convar_to_char(rep1, base_right, wp1.right.convar);
		Convar_to_char(rep2, base_right, wp2.right.convar);
		if (gausspoint == 2)
		{
			double var[2][4], der1[2][4], der2[2][4];
			for (int i = 0; i < 4; i++)
			{
				weno_5th_ao_with_df_2gauss(var[0][i], der1[0][i], der2[0][i],
					var[1][i], der1[1][i], der2[1][i], df2,
					ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h, 2);
			}
			for (int igauss = 0; igauss < gausspoint; igauss++)
			{
				Char_to_convar(re[igauss].right.convar, base_right, var[igauss]);
				Char_to_convar(re[igauss].right.der1y, base_right, der1[igauss]);
			}
		}
	}

	for (int i = 0; i < 4; i++)
	{
		double tmp[4];
		if (gausspoint == 2)
		{
			weno_5th_ao_with_df_2gauss(re[0].right.der1x[i], tmp[0], tmp[1],
				re[1].right.der1x[i], tmp[2], tmp[3], df2,
				wn2.right.der1x[i], wn1.right.der1x[i], w0.right.der1x[i], wp1.right.der1x[i], wp2.right.der1x[i], h, 1);
		}
	}


	for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
	{
		Check_Order_Reduce_by_Lambda_2D(re[num_gauss].left.is_reduce_order, re[num_gauss].left.convar);
		if (re[num_gauss].left.is_reduce_order == true)
		{
			for (int var = 0; var < 4; var++)
			{
				re[num_gauss].left.convar[var] = w0.left.convar[var];
				re[num_gauss].left.der1x[var] = w0.left.der1x[var];
				re[num_gauss].left.der1y[var] = 0.0;
			}
		}
	}

	for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
	{
		Check_Order_Reduce_by_Lambda_2D(re[num_gauss].right.is_reduce_order, re[num_gauss].right.convar);
		if (re[num_gauss].right.is_reduce_order == true)
		{
			for (int var = 0; var < 4; var++)
			{
				re[num_gauss].right.convar[var] = w0.right.convar[var];
				re[num_gauss].right.der1x[var] = w0.right.der1x[var];
				re[num_gauss].right.der1y[var] = 0.0;
			}
		}
	}
}

void weno_5th_ao_with_df_2gauss(double& g1, double& g1x, double& g1xx, double& g2, double& g2x, double& g2xx, double* df, double wn2, double wn1, double w0, double wp1, double wp2, double h, int order)
{
	//the parameter order constrols up to which order you want construct
	// order from 0, 1, 2
	if (order > 2 || order < 0)
	{
		cout << "invalid order input for the function " << __FUNCTION__ << endl;
		exit(0);
	}
	double dhi = 0.85;
	double dlo = 0.85;
	//-- - parameter of WENO-- -
	double beta[4], d[4], ww[4], alpha[4];
	double epsilonW = 1e-6;
	//-- - intermediate parameter-- -
	double p[4], px[4], pxx[4], tempvar;
	double sum_alpha;

	//three small stencil
	d[0] = (1.0 - dhi) * (1.0 - dlo) / 2.0;
	d[1] = (1.0 - dhi) * dlo;
	d[2] = (1.0 - dhi) * (1.0 - dlo) / 2.0;
	//one big stencil
	d[3] = dhi;

	beta[0] = 13.0 / 12.0 * pow((wn2 - 2.0 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4.0 * wn1 + 3.0 * w0), 2);
	beta[1] = 13.0 / 12.0 * pow((wn1 - 2.0 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
	beta[2] = 13.0 / 12.0 * pow((w0 - 2.0 * wp1 + wp2), 2) + 0.25 * pow((3.0 * w0 - 4.0 * wp1 + wp2), 2);
	beta[3] = 1.0 / 6.0 * (beta[0] + 4.0 * beta[1] + beta[2]) + abs(beta[0] - beta[2]);
	double tau5 = 1.0 / 3.0 * (abs(beta[3] - beta[0]) + abs(beta[3] - beta[1]) + abs(beta[3] - beta[2]));

	sum_alpha = 0.0;
	for (int i = 0; i < 4; i++)
	{
		double global_div = tau5 / (beta[i] + epsilonW);
		alpha[i] = d[i] * (1.0 + global_div * global_div);
		sum_alpha += alpha[i];
	}

	for (int k = 0; k < 4; k++)
	{
		ww[k] = alpha[k] / sum_alpha;
	}
	//-- - candidate polynomial-- -
	double sqrt3 = 1.732050807568877294;
	// gauss point 1
	p[0] = w0 + df[0] * (w0 - (sqrt3 * w0) / 4 + (4 * wn1 - wn2) / (4 * sqrt3) - w0);
	p[1] = w0 + df[1] * (w0 + (wn1 - wp1) / (4 * sqrt3) - w0);
	p[2] = w0 + df[2] * ((1.0 / 12.0) * (3 * (4 + sqrt3) * w0 + sqrt3 * (-4 * wp1 + wp2)) - w0);
	p[3] = w0 + df[3] * ((4314 * w0 + (4 + 500 * sqrt3) * wn1 - wn2 - 70 * sqrt3 * wn2 +
		4 * wp1 - 500 * sqrt3 * wp1 - wp2 + 70 * sqrt3 * wp2) / 4320 - w0);
	
	px[0] = -df[0] * ((-9 + sqrt3) * w0 - 2 * (-6 + sqrt3) * wn1 + (-3 + sqrt3) * wn2) / (6 * h);
	px[1] = -df[1] * (-2 * sqrt3 * w0 + (3 + sqrt3) * wn1 + (-3 + sqrt3) * wp1) / (6 * h);
	px[2] = -df[2] * ((9 + sqrt3) * w0 - 2 * (6 + sqrt3) * wp1 + (3 + sqrt3) * wp2) / (6 * h);
	px[3] = df[3] * (48 * sqrt3 * w0 - 72 * wn1 - 26 * sqrt3 * wn1 + 9 * wn2 +
		2 * sqrt3 * wn2 + 72 * wp1 - 26 * sqrt3 * wp1 - 9 * wp2 +
		2 * sqrt3 * wp2) / (108 * h);
	//-- - combination-- -
	g1 = 0.0;
	g1x = 0.0;
	double final_weight[4];
	final_weight[3] = ww[3] / d[3];
	for (int k = 0; k < 3; k++)
	{
		final_weight[k] = ww[k] - ww[3] / d[3] * d[k];
	}

	for (int k = 0; k < 4; k++)
	{
		g1 += final_weight[k] * p[k];
		g1x += final_weight[k] * px[k];
	}
	// gauss point 2
	p[0] = w0 + df[0] * ((1.0 / 12.0) * (3 * (4 + sqrt3) * w0 + sqrt3 * (-4 * wn1 + wn2)) - w0);
	p[1] = w0 + df[1] * (w0 + (-wn1 + wp1) / (4 * sqrt3) - w0);
	p[2] = w0 + df[2] * (w0 - (sqrt3 * w0) / 4 + (4 * wp1 - wp2) / (4 * sqrt3) - w0);
	p[3] = w0 + df[3] * ((4314 * w0 + (4 - 500 * sqrt3) * wn1 - wn2 + 70 * sqrt3 * wn2 +
		4 * wp1 + 500 * sqrt3 * wp1 - wp2 - 70 * sqrt3 * wp2) / 4320 - w0);

	px[0] = df[0] * ((9 + sqrt3) * w0 - 2 * (6 + sqrt3) * wn1 + (3 + sqrt3) * wn2) / (6 * h);
	px[1] = df[1] * (-2 * sqrt3 * w0 + (-3 + sqrt3) * wn1 + (3 + sqrt3) * wp1) / (6 * h);
	px[2] = df[2] * ((-9 + sqrt3) * w0 - 2 * (-6 + sqrt3) * wp1 + (-3 + sqrt3) * wp2) / (6 * h);
	px[3] = -df[3] * ((48 * sqrt3 * w0 + 72 * wn1 - 26 * sqrt3 * wn1 - 9 * wn2 +
		2 * sqrt3 * wn2 - 72 * wp1 - 26 * sqrt3 * wp1 + 9 * wp2 +
		2 * sqrt3 * wp2) / (108 * h));
	
	//-- - combination-- -
	g2 = 0.0;
	g2x = 0.0;

	for (int k = 0; k < 4; k++)
	{
		g2 += final_weight[k] * p[k];
		g2x += final_weight[k] * px[k];
	}
}

void WENO7_AO_with_df_tangent(Interface2d* left, Interface2d* right, Interface2d* down, Interface2d* up, Fluid2d* fluids, Block2d block)
{
	// along x direction tangitial recontruction,
	double alpha1[7], alpha2[7];
	alpha1[0] = fluids[-3].alpha_y;
	alpha1[1] = fluids[-2].alpha_y; alpha1[2] = fluids[-1].alpha_y;
	alpha1[3] = fluids[0].alpha_y;  alpha1[4] = fluids[1].alpha_y;
	alpha1[5] = fluids[2].alpha_y; alpha1[6] = fluids[3].alpha_y;
	alpha2[0] = fluids[block.ny - 3].alpha_y; 
	alpha2[1] = fluids[block.ny - 2].alpha_y; alpha2[2] = fluids[block.ny - 1].alpha_y;
	alpha2[3] = fluids[block.ny].alpha_y;  alpha2[4] = fluids[block.ny + 1].alpha_y;
	alpha2[5] = fluids[block.ny + 2].alpha_y; alpha2[6] = fluids[block.ny + 3].alpha_y;
	weno_7th_ao_with_df_tangential(right[0].gauss, right[-3].line, 
		right[-2].line, right[-1].line, right[0].line, right[1].line, right[2].line, right[3].line, alpha1, alpha2, right[0].length);

	//since we already do the coordinate transform, along y, no transform needed.
	alpha1[0] = fluids[3 * block.ny].alpha_x; 
	alpha1[1] = fluids[2 * block.ny].alpha_x; alpha1[2] = fluids[block.ny].alpha_x;
	alpha1[3] = fluids[0].alpha_x;  alpha1[4] = fluids[-block.ny].alpha_x;
	alpha1[5] = fluids[-2 * block.ny].alpha_x; alpha1[6] = fluids[-3 * block.ny].alpha_x;
	alpha2[0] = fluids[3 * block.ny + 1].alpha_x; 
	alpha2[1] = fluids[2 * block.ny + 1].alpha_x; alpha2[2] = fluids[block.ny + 1].alpha_x;
	alpha2[3] = fluids[1].alpha_x;  alpha2[4] = fluids[-block.ny + 1].alpha_x;
	alpha2[5] = fluids[-2 * block.ny + 1].alpha_x; alpha2[6] = fluids[-3 * block.ny + 1].alpha_x;
	weno_7th_ao_with_df_tangential(up[0].gauss, up[3 * (block.ny + 1)].line, up[2 * (block.ny + 1)].line, up[block.ny + 1].line,
		up[0].line, up[-(block.ny + 1)].line, up[-2 * (block.ny + 1)].line, up[-3 * (block.ny + 1)].line, alpha1, alpha2, up[0].length);
}

void weno_7th_ao_with_df_tangential(Recon2d* re, Recon2d& wn3, Recon2d& wn2, Recon2d& wn1, Recon2d& w0, Recon2d& wp1, Recon2d& wp2, Recon2d& wp3, double* alpha1, double* alpha2, double h)
{
	double df1[5], df2[5];
	double sum_df1[5], sum_df2[5];
	sum_df1[0] = alpha1[1] + alpha1[3];
	sum_df1[1] = alpha1[2] + alpha1[4];
	sum_df1[2] = alpha1[3] + alpha1[5];
	sum_df1[3] = alpha1[1] + alpha1[3] + alpha1[5];
	sum_df1[4] = alpha1[0] + alpha1[2] + alpha1[4] + alpha1[6];

	sum_df2[0] = alpha2[1] + alpha2[3];
	sum_df2[1] = alpha2[2] + alpha2[4];
	sum_df2[2] = alpha2[3] + alpha2[5];
	sum_df2[3] = alpha2[1] + alpha2[3] + alpha2[5];
	sum_df2[4] = alpha2[0] + alpha2[2] + alpha2[4] + alpha2[6];
	for (int i = 0; i < 5; i++)
	{
		df1[i] = sum_df1[i] < 2.0 ? 1.0 : 2.0 / sum_df1[i];
		df2[i] = sum_df2[i] < 2.0 ? 1.0 : 2.0 / sum_df2[i];
	}
	//lets first reconstruction the left value
	double ren3[4], ren2[4], ren1[4], re0[4], rep1[4], rep2[4], rep3[4];
	double base_left[4];
	double base_right[4];
	double wn1_primvar[4], w_primvar[4], wp1_primvar[4];
	Convar_to_primvar_2D(w_primvar, w0.left.convar);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = (w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		double tmp[2];
		for (int i = 0; i < 4; i++)
		{
			// 7th-order we use 4 gauss points
			weno_7th_ao_with_df_2gauss(re[0].left.convar[i], re[0].left.der1y[i],
					re[1].left.convar[i], re[1].left.der1y[i], re[2].left.convar[i], re[2].left.der1y[i], df1, wn3.left.convar[i],
					wn2.left.convar[i], wn1.left.convar[i], w0.left.convar[i], wp1.left.convar[i], wp2.left.convar[i], wp3.left.convar[i], h, 2);
		}
	}
	else
	{
		Convar_to_char(ren3, base_left, wn3.left.convar);
		Convar_to_char(ren2, base_left, wn2.left.convar);
		Convar_to_char(ren1, base_left, wn1.left.convar);
		Convar_to_char(re0, base_left, w0.left.convar);
		Convar_to_char(rep1, base_left, wp1.left.convar);
		Convar_to_char(rep2, base_left, wp2.left.convar);
		Convar_to_char(rep3, base_left, wp3.left.convar);
		double var[4][4], der1[4][4], der2[4][4];
		for (int i = 0; i < 4; i++)
		{
			weno_7th_ao_with_df_2gauss(var[0][i], der1[0][i],
				var[1][i], der1[1][i], var[2][i], der1[2][i], df1, ren3[i],
				ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], h, 2);
		}
		for (int igauss = 0; igauss < gausspoint; igauss++)
		{
			Char_to_convar(re[igauss].left.convar, base_left, var[igauss]);
			Char_to_convar(re[igauss].left.der1y, base_left, der1[igauss]);
		}

	}

	for (int i = 0; i < 4; i++)
	{
		double tmp[4];
		weno_7th_ao_with_df_2gauss(re[0].left.der1x[i], tmp[0],
				re[1].left.der1x[i], tmp[1], re[2].left.der1x[i], tmp[2], df1, wn3.left.der1x[i], 
				wn2.left.der1x[i], wn1.left.der1x[i], w0.left.der1x[i], wp1.left.der1x[i], wp2.left.der1x[i], wp3.left.der1x[i], h, 1);
	}

	//then let's construct the right part....
	Convar_to_primvar_2D(w_primvar, w0.right.convar);
	for (int i = 0; i < 4; i++)
	{
		base_right[i] = (w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		double tmp[2];
		for (int i = 0; i < 4; i++)
		{
			weno_7th_ao_with_df_2gauss(re[0].right.convar[i], re[0].right.der1y[i],
					re[1].right.convar[i], re[1].right.der1y[i], re[2].right.convar[i], re[2].right.der1y[i], df2, wn3.right.convar[i], 
					wn2.right.convar[i], wn1.right.convar[i], w0.right.convar[i], wp1.right.convar[i], wp2.right.convar[i], wp3.right.convar[i], h, 2);
		}
	}
	else
	{
		Convar_to_char(ren3, base_right, wn3.right.convar);
		Convar_to_char(ren2, base_right, wn2.right.convar);
		Convar_to_char(ren1, base_right, wn1.right.convar);
		Convar_to_char(re0, base_right, w0.right.convar);
		Convar_to_char(rep1, base_right, wp1.right.convar);
		Convar_to_char(rep2, base_right, wp2.right.convar);
		Convar_to_char(rep3, base_right, wp3.right.convar);
		double var[4][4], der1[4][4], der2[4][4];
		for (int i = 0; i < 4; i++)
		{
			weno_7th_ao_with_df_2gauss(var[0][i], der1[0][i], 
				var[1][i], der1[1][i], var[2][i], der1[2][i], df2, ren3[i],
				ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], h, 2);
		}
		for (int igauss = 0; igauss < gausspoint; igauss++)
		{
			Char_to_convar(re[igauss].right.convar, base_right, var[igauss]);
			Char_to_convar(re[igauss].right.der1y, base_right, der1[igauss]);
		}
	}

	for (int i = 0; i < 4; i++)
	{
		double tmp[4];
		weno_7th_ao_with_df_2gauss(re[0].right.der1x[i], tmp[0],
				re[1].right.der1x[i], tmp[1], re[2].right.der1x[i], tmp[2], df2, wn3.right.der1x[i], 
				wn2.right.der1x[i], wn1.right.der1x[i], w0.right.der1x[i], wp1.right.der1x[i], wp2.right.der1x[i], wp3.right.der1x[i], h, 1);
	}


	for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
	{
		Check_Order_Reduce_by_Lambda_2D(re[num_gauss].left.is_reduce_order, re[num_gauss].left.convar);
		if (re[num_gauss].left.is_reduce_order == true)
		{
			for (int var = 0; var < 4; var++)
			{
				re[num_gauss].left.convar[var] = w0.left.convar[var];
				re[num_gauss].left.der1x[var] = w0.left.der1x[var];
				re[num_gauss].left.der1y[var] = 0.0;
			}
		}
	}

	for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
	{
		Check_Order_Reduce_by_Lambda_2D(re[num_gauss].right.is_reduce_order, re[num_gauss].right.convar);
		if (re[num_gauss].right.is_reduce_order == true)
		{
			for (int var = 0; var < 4; var++)
			{
				re[num_gauss].right.convar[var] = w0.right.convar[var];
				re[num_gauss].right.der1x[var] = w0.right.der1x[var];
				re[num_gauss].right.der1y[var] = 0.0;
			}
		}
	}
}

void weno_7th_ao_with_df_2gauss(double& g1, double& g1x, double& g2, double& g2x, double& g3, double& g3x, double* df, double wn3, double wn2, double wn1, double w0, double wp1, double wp2, double wp3, double h, int order)
{
	//the parameter order constrols up to which order you want construct
	// order from 0, 1, 2
	if (order > 2 || order < 0)
	{
		cout << "invalid order input for the function " << __FUNCTION__ << endl;
		exit(0);
	}
	double dhi = 0.85;
	double dlo = 0.85;
	double davg= 0.85;
	// P(7,5,3) -- one 7th-order stencil, one 5th-order stencil, three 3rd-order stencil
	//-- - parameter of WENO-- -
	double beta[5], d[5], ww[5], alpha[5];
	double epsilonW = 1e-6;
	//-- - intermediate parameter-- -
	double p[5], px[5], pxx[5], tempvar;
	double sum_alpha;

	//three 3rd-order stencil
	d[0] = (1.0 - dhi) * (1.0 - dlo) * (1.0 - davg) / 2.0;
	d[1] = (1.0 - dhi) * (1.0 - davg) * dlo;
	d[2] = (1.0 - dhi) * (1.0 - dlo) * (1.0 - davg) / 2.0;
	//one 5th-order stencil
	d[3] = (1.0 - dhi) * davg;
	//one 7th-order stencil
	d[4] = dhi;

	beta[0] = 13.0 / 12.0 * pow((wn2 - 2.0 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4.0 * wn1 + 3.0 * w0), 2);
	beta[1] = 13.0 / 12.0 * pow((wn1 - 2.0 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
	beta[2] = 13.0 / 12.0 * pow((w0 - 2.0 * wp1 + wp2), 2) + 0.25 * pow((3.0 * w0 - 4.0 * wp1 + wp2), 2);

	beta[3] = 1.0 / 6.0 * (beta[0] + 4.0 * beta[1] + beta[2]) + abs(beta[0] - beta[2]);
	// For 7th-order stencil, the smoothness indicator
	double a1, a2, a3;
	double beta4[4];
	a1 = (-19 * wn3 + 87 * wn2 - 177 * wn1 + 109 * w0) / 60.0;
	a2 = (-wn3 + 4 * wn2 - 5 * wn1 + 2 * w0) / 2;
	a3 = (-wn3 + 3 * wn2 - 3 * wn1 + w0) / 6.0;
	beta4[0] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (11 * wn2 - 63 * wn1 + 33 * w0 + 19 * wp1) / 60.0;
	a2 = (wn1 - 2 * w0 + wp1) / 2;
	a3 = (-wn2 + 3 * wn1 - 3 * w0 + wp1) / 6.0;
	beta4[1] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (-19 * wn1 - 33 * w0 + 63 * wp1 - 11 * wp2) / 60.0;
	a2 = (wn1 - 2 * w0 + wp1) / 2;
	a3 = (-wn1 + 3 * w0 - 3 * wp1 + wp2) / 6.0;
	beta4[2] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (-109 * w0 + 177 * wp1 - 87 * wp2 + 19 * wp3) / 60.0;
	a2 = (2 * w0 - 5 * wp1 + 4 * wp2 - wp3) / 2;
	a3 = (-w0 + 3 * wp1 - 3 * wp2 + wp3) / 6.0;
	beta4[3] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	beta[4] = 1.0 / 20.0 * (beta4[0] + 9.0 * beta4[1] + 9.0 * beta4[2] + beta4[3])
			+ abs(beta4[0] - beta4[1] - beta4[2] + beta4[3]);
	
	double tau5 = 0.25 * (abs(beta[4] - beta[0]) + abs(beta[4] - beta[1]) + abs(beta[4] - beta[2]) + abs(beta[4] - beta[3]));

	sum_alpha = 0.0;
	for (int i = 0; i < 5; i++)
	{
		double global_div = tau5 / (beta[i] + epsilonW);
		alpha[i] = d[i] * (1.0 + global_div * global_div);
		sum_alpha += alpha[i];
	}

	for (int k = 0; k < 5; k++)
	{
		ww[k] = alpha[k] / sum_alpha;
	}
	double sqrt15 = 3.872983346207417;
	// gauss point 1
	p[0] = w0 + df[0] * ((23.0 * w0 + 2.0 * wn1 - wn2) / 24.0 - w0);
	p[1] = w0 + df[1] * ((26.0 * w0 - wn1 - wp1) / 24.0 - w0);
	p[2] = w0 + df[2] * ((23.0 * w0 + 2.0 * wp1 - wp2) / 24.0 - w0);
	p[3] = w0 + df[3] * ((2134.0 * w0 - 116.0 * wn1 + 9.0 * wn2 - 116.0 * wp1 + 9.0 * wp2) / 1920.0 - w0);
	p[4] = w0 + df[4] * ((121004 * w0 - 7621 * wn1 + 954 * wn2 - 75 * wn3 - 7621 * wp1 + 954 * wp2 - 75 * wp3) / 107520.0 - w0);
	
	px[0] = df[0] * (3.0 * w0 - 4.0 * wn1 + wn2) / (2.0 * h);
	px[1] = df[1] * (-wn1 + wp1) / (2.0 * h);
	px[2] = df[2] * (-3.0 * w0 + 4.0 * wp1 - wp2) / (2.0 * h);
	px[3] = df[3] * (-34.0 * wn1 + 5.0 * wn2 + 34.0 * wp1 - 5.0 * wp2) / (48.0 * h);
	px[4] = df[4] * (-9455 * wn1 + 2236 * wn2 - 259 * wn3 + 9455 * wp1 - 2236 * wp2 + 259 * wp3) / (11520.0 * h);
	//-- - combination-- -
	g1 = 0.0;
	g1x = 0.0;
	double final_weight[5];
	final_weight[4] = ww[4] / d[4];
	for (int k = 0; k < 4; k++)
	{
		final_weight[k] = ww[k] - ww[4] / d[4] * d[k];
	}

	for (int k = 0; k < 5; k++)
	{
		g1 += final_weight[k] * p[k];
		g1x += final_weight[k] * px[k];
	}
	// gauss point 2
	p[0] = w0 + df[0] * (((62 - 9 * sqrt15) * w0 + 4 * (-1 + 3 * sqrt15) * wn1 + (2 - 3 * sqrt15) * wn2) / 60.0 - w0);
	p[1] = w0 + df[1] * ((56 * w0 + (2 + 3 * sqrt15) * wn1 + (2 - 3 * sqrt15) * wp1) / 60.0 - w0);
	p[2] = w0 + df[2] * (((62 + 9 * sqrt15) * w0 - 4 * (1 + 3 * sqrt15) * wp1 + (2 + 3 * sqrt15) * wp2) / 60.0 - w0);
	p[3] = w0 + df[3] * ((2186 * w0 + 4 * (29 + 41 * sqrt15) * wn1 - 9 * wn2 - 22 * sqrt15 * wn2 + 116 * wp1 - 164 * sqrt15 * wp1 - 9 * wp2 + 22 * sqrt15 * wp2) / 2400.0 - w0);
	p[4] = w0 + df[4] * ((4534440 * w0 + 5 * (57144 + 78421 * sqrt15) * wn1 - 35748 * wn2 - 84364 * sqrt15 * wn2 + 2808 * wn3 + 9541 * sqrt15 * wn3 + 285720 * wp1 - 
 		   392105 * sqrt15 * wp1 - 35748 * wp2 + 84364 * sqrt15 * wp2 + 2808 * wp3 - 9541 * sqrt15 * wp3) / 5040000.0 - w0);

	px[0] = df[0] * (-((-15 + sqrt15) * w0) + 2 * (-10 + sqrt15) * wn1 - (-5 + sqrt15) * wn2) / (10.0 * h);
	px[1] = df[1] * (2 * sqrt15 * w0 - (5 + sqrt15) * wn1 - (-5 + sqrt15) * wp1) / (h * 10.0);
	px[2] = df[2] * (-((15 + sqrt15) * w0) + 2 * (10 + sqrt15) * wp1 - (5 + sqrt15) * wp2) / (10.0 * h);
	px[3] = df[3] * (78 * sqrt15 * w0 - 2 * (95 + 21 * sqrt15) * wn1 + 20 * wn2 + 3 * sqrt15 * wn2 + 190 * wp1 - 42 * sqrt15 * wp1 - 20 * wp2 + 3 * sqrt15 * wp2) / (300.0 * h);
	px[4] = df[4] * (103860 * sqrt15 * w0 - 5 * (49925 + 11619 * sqrt15) * wn1 + 41300 * wn2 + 6678 * sqrt15 * wn2 - 4325 * wn3 - 513 * sqrt15 * wn3 + 249625 * wp1 - 
 		    58095 * sqrt15 * wp1 - 41300 * wp2 + 6678 * sqrt15 * wp2 + 4325 * wp3 - 513 * sqrt15 * wp3) / (360000.0 * h);
	//-- - combination-- -
	g2 = 0.0;
	g2x = 0.0;

	for (int k = 0; k < 5; k++)
	{
		g2 += final_weight[k] * p[k];
		g2x += final_weight[k] * px[k];
	}

	// gauss point 3
	p[0] = w0 + df[0] * (((62 + 9 * sqrt15) * w0 - 4 * (1 + 3 * sqrt15) * wn1 + (2 + 3 * sqrt15) * wn2) / 60.0 - w0);
	p[1] = w0 + df[1] * ((56 * w0 + (2 - 3 * sqrt15) * wn1 + (2 + 3 * sqrt15) * wp1) / 60.0 - w0);
	p[2] = w0 + df[2] * (((62 - 9 * sqrt15) * w0 + 4 * (-1 + 3 * sqrt15) * wp1 + (2 - 3 * sqrt15) * wp2) / 60.0 - w0);
	p[3] = w0 + df[3] * ((2186 * w0 - 4 * (-29 + 41 * sqrt15) * wn1 - 9 * wn2 + 22 * sqrt15 * wn2 + 116 * wp1 + 164 * sqrt15 * wp1 - 9 * wp2 - 22 * sqrt15 * wp2) / 2400.0 - w0);
	p[4] = w0 + df[4] * ((4534440 * w0 - 5 * (-57144 + 78421 * sqrt15) * wn1 - 35748 * wn2 + 84364 * sqrt15 * wn2 + 2808 * wn3 - 9541 * sqrt15 * wn3 + 285720 * wp1 + 
 		   392105 * sqrt15 * wp1 - 35748 * wp2 - 84364 * sqrt15 * wp2 + 2808 * wp3 + 9541 * sqrt15 * wp3) / 5040000.0 - w0);

	px[0] = df[0] * ((15 + sqrt15) * w0 - 2 * (10 + sqrt15) * wn1 + (5 + sqrt15) * wn2) / (10.0 * h);
	px[1] = df[1] * (-2 * sqrt15 * w0 - (5 - sqrt15) * wn1 + (5 + sqrt15) * wp1) / (h * 10.0);
	px[2] = df[2] * ((-15 + sqrt15) * w0 - 2 * (-10 + sqrt15) * wp1 + (-5 + sqrt15) * wp2) / (10.0 * h);
	px[3] = df[3] * (-78 * sqrt15 * w0 + 2 * (-95 + 21 * sqrt15) * wn1 + 20 * wn2 - 3 * sqrt15 * wn2 + 190 * wp1 + 42 * sqrt15 * wp1 - 20 * wp2 - 3 * sqrt15 * wp2) / (300.0 * h);
	px[4] = df[4] * (-103860 * sqrt15 * w0 + 5 * (-49925 + 11619 * sqrt15) * wn1 + 41300 * wn2 - 6678 * sqrt15 * wn2 - 4325 * wn3 + 513 * sqrt15 * wn3 + 249625 * wp1 + 
 			58095 * sqrt15 * wp1 - 41300 * wp2 - 6678 * sqrt15 * wp2 + 4325 * wp3 + 513 * sqrt15 * wp3) / (360000.0 * h);
	//-- - combination-- -
	g3 = 0.0;
	g3x = 0.0;

	for (int k = 0; k < 5; k++)
	{
		g3 += final_weight[k] * p[k];
		g3x += final_weight[k] * px[k];
	}
}

void WENO9_AO_with_df_tangent(Interface2d* left, Interface2d* right, Interface2d* down, Interface2d* up, Fluid2d* fluids, Block2d block)
{
	// along x direction tangitial recontruction,
	double alpha1[9], alpha2[9];
	alpha1[0] = fluids[-4].alpha_y; alpha1[1] = fluids[-3].alpha_y;
	alpha1[2] = fluids[-2].alpha_y; alpha1[3] = fluids[-1].alpha_y;
	alpha1[4] = fluids[0].alpha_y;  alpha1[5] = fluids[1].alpha_y;
	alpha1[6] = fluids[2].alpha_y; alpha1[7] = fluids[3].alpha_y;
	alpha1[8] = fluids[4].alpha_y;
	alpha2[0] = fluids[block.ny - 4].alpha_y; alpha2[1] = fluids[block.ny - 3].alpha_y; 
	alpha2[2] = fluids[block.ny - 2].alpha_y; alpha2[3] = fluids[block.ny - 1].alpha_y;
	alpha2[4] = fluids[block.ny].alpha_y;  alpha2[5] = fluids[block.ny + 1].alpha_y;
	alpha2[6] = fluids[block.ny + 2].alpha_y; alpha2[7] = fluids[block.ny + 3].alpha_y;
	alpha2[8] = fluids[block.ny + 4].alpha_y;
	weno_9th_ao_with_df_tangential(right[0].gauss, right[-4].line, right[-3].line, 
		right[-2].line, right[-1].line, right[0].line, right[1].line, right[2].line, right[3].line, right[4].line, alpha1, alpha2, right[0].length);

	//since we already do the coordinate transform, along y, no transform needed.
	alpha1[0] = fluids[4 * block.ny].alpha_x; alpha1[1] = fluids[3 * block.ny].alpha_x; 
	alpha1[2] = fluids[2 * block.ny].alpha_x; alpha1[3] = fluids[block.ny].alpha_x;
	alpha1[4] = fluids[0].alpha_x;  alpha1[5] = fluids[-block.ny].alpha_x;
	alpha1[6] = fluids[-2 * block.ny].alpha_x; alpha1[7] = fluids[-3 * block.ny].alpha_x;
	alpha1[8] = fluids[-4 * block.ny].alpha_x;
	alpha2[0] = fluids[4 * block.ny + 1].alpha_x; alpha2[1] = fluids[3 * block.ny + 1].alpha_x;
	alpha2[2] = fluids[2 * block.ny + 1].alpha_x; alpha2[3] = fluids[block.ny + 1].alpha_x;
	alpha2[4] = fluids[1].alpha_x;  alpha2[5] = fluids[-block.ny + 1].alpha_x;
	alpha2[6] = fluids[-2 * block.ny + 1].alpha_x; alpha2[7] = fluids[-3 * block.ny + 1].alpha_x;
	alpha2[8] = fluids[-4 * block.ny + 1].alpha_x;
	weno_9th_ao_with_df_tangential(up[0].gauss, up[4 * (block.ny + 1)].line, up[3 * (block.ny + 1)].line, up[2 * (block.ny + 1)].line, up[block.ny + 1].line,
		up[0].line, up[-(block.ny + 1)].line, up[-2 * (block.ny + 1)].line, up[-3 * (block.ny + 1)].line, up[-4 * (block.ny + 1)].line, alpha1, alpha2, up[0].length);
}

void weno_9th_ao_with_df_tangential(Recon2d* re, Recon2d& wn4, Recon2d& wn3, Recon2d& wn2, Recon2d& wn1, Recon2d& w0, Recon2d& wp1, Recon2d& wp2, Recon2d& wp3, Recon2d& wp4, double* alpha1, double* alpha2, double h)
{
	double df1[6], df2[6];
	double sum_df1[6], sum_df2[6];
	sum_df1[0] = alpha1[2] + alpha1[4];
	sum_df1[1] = alpha1[3] + alpha1[5];
	sum_df1[2] = alpha1[4] + alpha1[6];
	sum_df1[3] = alpha1[2] + alpha1[4] + alpha1[6];
	sum_df1[4] = alpha1[1] + alpha1[3] + alpha1[5] + alpha1[7];
	sum_df1[5] = alpha1[0] + alpha1[2] + alpha1[4] + alpha1[6] + alpha1[8];

	sum_df2[0] = alpha2[2] + alpha2[4];
	sum_df2[1] = alpha2[3] + alpha2[5];
	sum_df2[2] = alpha2[4] + alpha2[6];
	sum_df2[3] = alpha2[2] + alpha2[4] + alpha2[6];
	sum_df2[4] = alpha2[1] + alpha2[3] + alpha2[5] + alpha2[7];
	sum_df2[5] = alpha2[0] + alpha2[2] + alpha2[4] + alpha2[6] + alpha2[8];
	for (int i = 0; i < 6; i++)
	{
		df1[i] = sum_df1[i] < 2.0 ? 1.0 : 2.0 / sum_df1[i];
		df2[i] = sum_df2[i] < 2.0 ? 1.0 : 2.0 / sum_df2[i];
	}
	//lets first reconstruction the left value
	double ren4[4], ren3[4], ren2[4], ren1[4], re0[4], rep1[4], rep2[4], rep3[4], rep4[4];
	double base_left[4];
	double base_right[4];
	double wn1_primvar[4], w_primvar[4], wp1_primvar[4];
	Convar_to_primvar_2D(w_primvar, w0.left.convar);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = (w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		double tmp[2];
		for (int i = 0; i < 4; i++)
		{
			// 9th-order we use 4 gauss points
			weno_9th_ao_with_df_2gauss(re[0].left.convar[i], re[0].left.der1y[i],
					re[1].left.convar[i], re[1].left.der1y[i], re[2].left.convar[i], re[2].left.der1y[i], re[3].left.convar[i], re[3].left.der1y[i], df1, wn4.left.convar[i], wn3.left.convar[i],
					wn2.left.convar[i], wn1.left.convar[i], w0.left.convar[i], wp1.left.convar[i], wp2.left.convar[i], wp3.left.convar[i], wp4.left.convar[i], h, 2);
		}
	}
	else
	{
		Convar_to_char(ren4, base_left, wn4.left.convar);
		Convar_to_char(ren3, base_left, wn3.left.convar);
		Convar_to_char(ren2, base_left, wn2.left.convar);
		Convar_to_char(ren1, base_left, wn1.left.convar);
		Convar_to_char(re0, base_left, w0.left.convar);
		Convar_to_char(rep1, base_left, wp1.left.convar);
		Convar_to_char(rep2, base_left, wp2.left.convar);
		Convar_to_char(rep3, base_left, wp3.left.convar);
		Convar_to_char(rep4, base_left, wp4.left.convar);
		double var[4][4], der1[4][4], der2[4][4];
		for (int i = 0; i < 4; i++)
		{
			weno_9th_ao_with_df_2gauss(var[0][i], der1[0][i],
				var[1][i], der1[1][i], var[2][i], der1[2][i], var[3][i], der1[3][i], df1, ren4[i], ren3[i],
				ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], rep4[i], h, 2);
		}
		for (int igauss = 0; igauss < gausspoint; igauss++)
		{
			Char_to_convar(re[igauss].left.convar, base_left, var[igauss]);
			Char_to_convar(re[igauss].left.der1y, base_left, der1[igauss]);
		}

	}

	for (int i = 0; i < 4; i++)
	{
		double tmp[4];
		weno_9th_ao_with_df_2gauss(re[0].left.der1x[i], tmp[0],
				re[1].left.der1x[i], tmp[1], re[2].left.der1x[i], tmp[2], re[3].left.der1x[i], tmp[3], df1, wn4.left.der1x[i], wn3.left.der1x[i], 
				wn2.left.der1x[i], wn1.left.der1x[i], w0.left.der1x[i], wp1.left.der1x[i], wp2.left.der1x[i], wp3.left.der1x[i], wp4.left.der1x[i], h, 1);
	}

	//then let's construct the right part....
	Convar_to_primvar_2D(w_primvar, w0.right.convar);
	for (int i = 0; i < 4; i++)
	{
		base_right[i] = (w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		double tmp[2];
		for (int i = 0; i < 4; i++)
		{
			weno_9th_ao_with_df_2gauss(re[0].right.convar[i], re[0].right.der1y[i],
					re[1].right.convar[i], re[1].right.der1y[i], re[2].right.convar[i], re[2].right.der1y[i], re[3].right.convar[i], re[3].right.der1y[i], df2, wn4.right.convar[i], wn3.right.convar[i], 
					wn2.right.convar[i], wn1.right.convar[i], w0.right.convar[i], wp1.right.convar[i], wp2.right.convar[i], wp3.right.convar[i], wp4.right.convar[i], h, 2);
		}
	}
	else
	{
		Convar_to_char(ren3, base_right, wn3.right.convar);
		Convar_to_char(ren2, base_right, wn2.right.convar);
		Convar_to_char(ren1, base_right, wn1.right.convar);
		Convar_to_char(re0, base_right, w0.right.convar);
		Convar_to_char(rep1, base_right, wp1.right.convar);
		Convar_to_char(rep2, base_right, wp2.right.convar);
		Convar_to_char(rep3, base_right, wp3.right.convar);
		double var[4][4], der1[4][4], der2[4][4];
		for (int i = 0; i < 4; i++)
		{
			weno_9th_ao_with_df_2gauss(var[0][i], der1[0][i], 
				var[1][i], der1[1][i], var[2][i], der1[2][i], var[3][i], der1[3][i], df2, ren4[i], ren3[i],
				ren2[i], ren1[i], re0[i], rep1[i], rep2[i], rep3[i], rep4[i], h, 2);
		}
		for (int igauss = 0; igauss < gausspoint; igauss++)
		{
			Char_to_convar(re[igauss].right.convar, base_right, var[igauss]);
			Char_to_convar(re[igauss].right.der1y, base_right, der1[igauss]);
		}
	}

	for (int i = 0; i < 4; i++)
	{
		double tmp[4];
		weno_9th_ao_with_df_2gauss(re[0].right.der1x[i], tmp[0],
				re[1].right.der1x[i], tmp[1], re[2].right.der1x[i], tmp[2], re[3].right.der1x[i], tmp[3], df2, wn4.right.der1x[i], wn3.right.der1x[i], 
				wn2.right.der1x[i], wn1.right.der1x[i], w0.right.der1x[i], wp1.right.der1x[i], wp2.right.der1x[i], wp3.right.der1x[i], wp4.right.der1x[i], h, 1);
	}


	for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
	{
		Check_Order_Reduce_by_Lambda_2D(re[num_gauss].left.is_reduce_order, re[num_gauss].left.convar);
		if (re[num_gauss].left.is_reduce_order == true)
		{
			for (int var = 0; var < 4; var++)
			{
				re[num_gauss].left.convar[var] = w0.left.convar[var];
				re[num_gauss].left.der1x[var] = w0.left.der1x[var];
				re[num_gauss].left.der1y[var] = 0.0;
			}
		}
	}

	for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
	{
		Check_Order_Reduce_by_Lambda_2D(re[num_gauss].right.is_reduce_order, re[num_gauss].right.convar);
		if (re[num_gauss].right.is_reduce_order == true)
		{
			for (int var = 0; var < 4; var++)
			{
				re[num_gauss].right.convar[var] = w0.right.convar[var];
				re[num_gauss].right.der1x[var] = w0.right.der1x[var];
				re[num_gauss].right.der1y[var] = 0.0;
			}
		}
	}
}

void weno_9th_ao_with_df_2gauss(double& g1, double& g1x, double& g2, double& g2x, double& g3, double& g3x, double& g4, double& g4x, double* df, double wn4, double wn3, double wn2, double wn1, double w0, double wp1, double wp2, double wp3, double wp4, double h, int order)
{
	//the parameter order constrols up to which order you want construct
	// order from 0, 1, 2
	if (order > 2 || order < 0)
	{
		cout << "invalid order input for the function " << __FUNCTION__ << endl;
		exit(0);
	}
	double dhi = 0.85;
	double dlo = 0.85;
	double davg1= 0.85;
	double davg2 = 0.85;
	// P(9, 7,5,3) -- one 9th-order stencil, one 7th-order stencil, one 5th-order stencil, three 3rd-order stencil
	//-- - parameter of WENO-- -
	double beta[6], d[6], ww[6], alpha[6];
	double epsilonW = 1e-6;
	//-- - intermediate parameter-- -
	double p[6], px[6], pxx[6], tempvar;
	double sum_alpha;

	//three 3rd-order stencil
	d[0] = (1.0 - dhi) * (1.0 - dlo) * (1.0 - davg1) * (1.0 - davg2) / 2.0;
	d[1] = (1.0 - dhi) * (1.0 - davg1) * (1.0 - davg2) * dlo;
	d[2] = (1.0 - dhi) * (1.0 - dlo) * (1.0 - davg1) * (1.0 - davg2) / 2.0;
	//one 5th-order stencil
	d[3] = (1.0 - dhi) * (1.0 - davg2) * davg1;
	//one 7th-order stencil
	d[4] = (1.0 - dhi) * davg2;
	//one 9th-order stencil
	d[5] = dhi;

	beta[0] = 13.0 / 12.0 * pow((wn2 - 2.0 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4.0 * wn1 + 3.0 * w0), 2);
	beta[1] = 13.0 / 12.0 * pow((wn1 - 2.0 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
	beta[2] = 13.0 / 12.0 * pow((w0 - 2.0 * wp1 + wp2), 2) + 0.25 * pow((3.0 * w0 - 4.0 * wp1 + wp2), 2);
	//beta[3] = 1.0 / 6.0 * (beta[0] + 4.0 * beta[1] + beta[2]) + abs(beta[0] - beta[2]);

	//For 7th-order stencil, the smoothness indicator
	double a1, a2, a3, a4;
	double beta4[4];
	a1 = (-19 * wn3 + 87 * wn2 - 177 * wn1 + 109 * w0) / 60.0;
	a2 = (-wn3 + 4 * wn2 - 5 * wn1 + 2 * w0) / 2;
	a3 = (-wn3 + 3 * wn2 - 3 * wn1 + w0) / 6.0;
	beta4[0] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (11 * wn2 - 63 * wn1 + 33 * w0 + 19 * wp1) / 60.0;
	a2 = (wn1 - 2 * w0 + wp1) / 2;
	a3 = (-wn2 + 3 * wn1 - 3 * w0 + wp1) / 6.0;
	beta4[1] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (-19 * wn1 - 33 * w0 + 63 * wp1 - 11 * wp2) / 60.0;
	a2 = (wn1 - 2 * w0 + wp1) / 2;
	a3 = (-wn1 + 3 * w0 - 3 * wp1 + wp2) / 6.0;
	beta4[2] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	a1 = (-109 * w0 + 177 * wp1 - 87 * wp2 + 19 * wp3) / 60.0;
	a2 = (2 * w0 - 5 * wp1 + 4 * wp2 - wp3) / 2;
	a3 = (-w0 + 3 * wp1 - 3 * wp2 + wp3) / 6.0;
	beta4[3] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13 / 3 * a2 * a2
			 + 781 / 20 * a3 * a3;

	beta[4] = 1.0 / 20.0 * (beta4[0] + 9.0 * beta4[1] + 9.0 * beta4[2] + beta4[3])
			+ abs(beta4[0] - beta4[1] - beta4[2] + beta4[3]);

	// For 9th-order stencil, the smoothness indicator
	double beta5[5];
	a1 = (-462 * wn1 + 336 * wn2 - 146 * wn3 + 27 * wn4 + 245 * w0) / 120.0;
	a2 = (-240 * wn1 + 262 * wn2 - 128 * wn3 + 25 * wn4 + 81 * w0) / 56.0;
	a3 = (-18 * wn1 + 24 * wn2 - 14 * wn3 + 3 * wn4 + 5 * w0) / 12.0;
	a4 = (-4 * wn1 + 6 * wn2 - 4 * wn3 + wn4 + w0) / 24.0;
	beta5[0] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13.0 / 3.0 * (a2 + 123.0 / 455.0 * a4) * (a2 + 123.0 / 455.0 * a4) 
		     + 781.0 / 20.0 * a3 * a3 + 1421461.0 / 2275.0 * a4 * a4;

	a1 = (-192 * wn1 + 66 * wn2 - 11 * wn3 + 110 * w0 + 27 * wp1) / 120;
	a2 = (10 * wn1 + 12 * wn2 - 3 * wn3 - 44 * w0 + 25 * wp1) / 56;
	a3 = (12 * wn1 - 6 * wn2 + wn3 - 10 * w0 + 3 * wp1) / 12;
	a4 = (6 * wn1 - 4 * wn2 + wn3 - 4 * w0 + wp1) / 24;
	beta5[1] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13.0 / 3.0 * (a2 + 123.0 / 455.0 * a4) * (a2 + 123.0 / 455.0 * a4) 
		     + 781.0 / 20.0 * a3 * a3 + 1421461.0 / 2275.0 * a4 * a4;

	a1 = (-82 * wn1 + 11 * wn2 + 82 * wp1 - 11 * wp2)/120;
	a2 = (40 * wn1 - 3 * wn2 - 74 * w0 + 40 * wp1 - 3 * wp2)/56;
	a3 = (2 * wn1 - wn2 - 2 * wp1 + wp2)/12;
	a4 = (-4 * wn1 + wn2 + 6 * w0 - 4 * wp1 + wp2)/24;
	beta5[2] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13.0 / 3.0 * (a2 + 123.0 / 455.0 * a4) * (a2 + 123.0 / 455.0 * a4) 
		     + 781.0 / 20.0 * a3 * a3 + 1421461.0 / 2275.0 * a4 * a4;
	
	a1 = (-27 * wn1 - 110 * w0 + 192 * wp1 - 66 * wp2 + 11 * wp3)/120;
	a2 = (25 * wn1 - 44 * w0 + 10 * wp1 + 12 * wp2 - 3 * wp3)/56;
	a3 = (-3 * wn1 + 10 * w0 - 12 * wp1 + 6 * wp2 - wp3)/12;
	a4 = (wn1 - 4 * w0 + 6 * wp1 - 4 * wp2 + wp3)/24;
	beta5[3] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13.0 / 3.0 * (a2 + 123.0 / 455.0 * a4) * (a2 + 123.0 / 455.0 * a4) 
		     + 781.0 / 20.0 * a3 * a3 + 1421461.0 / 2275.0 * a4 * a4;

	a1 = (-245 * w0 + 462 * wp1 - 336 * wp2 + 146 * wp3 - 27 * wp4)/120;
	a2 = (81 * w0 - 240 * wp1 + 262 * wp2 - 128 * wp3 + 25 * wp4)/56;
	a3 = (-5 * w0 + 18 * wp1 - 24 * wp2 + 14 * wp3 - 3 * wp4)/12;
	a4 = (w0 - 4 * wp1 + 6 * wp2 - 4 * wp3 + wp4)/24;
	beta5[4] = (a1 + a3 / 10.0) * (a1 + a3 / 10.0) + 13.0 / 3.0 * (a2 + 123.0 / 455.0 * a4) * (a2 + 123.0 / 455.0 * a4) 
		     + 781.0 / 20.0 * a3 * a3 + 1421461.0 / 2275.0 * a4 * a4;
	
	// beta[3] represent the 5th-order smoothness indicator
	beta[3] = beta5[2];
	beta[5] = 1.0 / 70.0 * (beta5[0] + 16.0 * beta5[1] + 36.0 * beta5[2] + 16.0 * beta5[3] + beta5[4]) + abs(beta5[0] - beta5[1] - beta5[3] + beta5[4]);
	
	double tau5 = 0.2 * (abs(beta[5] - beta[0]) + abs(beta[5] - beta[1]) + abs(beta[5] - beta[2]) + abs(beta[5] - beta[3]) + abs(beta[5] - beta[4]));

	sum_alpha = 0.0;
	for (int i = 0; i < 6; i++)
	{
		double global_div = tau5 / (beta[i] + epsilonW);
		alpha[i] = d[i] * (1.0 + global_div * global_div);
		sum_alpha += alpha[i];
	}

	for (int k = 0; k < 6; k++)
	{
		ww[k] = alpha[k] / sum_alpha;
	}
	double final_weight[6];
	final_weight[5] = ww[5] / d[5];
	for (int k = 0; k < 5; k++)
	{
		final_weight[k] = ww[k] - ww[5] / d[5] * d[k];
	}

	double sqrt5 = 2.236067977499790;
	// gauss point 1
	double x = -1.0;
	Polynomial_9th(p, px, df, wn4, wn3, wn2, wn1, w0, wp1, wp2, wp3, wp4, x, h);
	
	//-- - combination-- -
	g1 = 0.0;
	g1x = 0.0;

	for (int k = 0; k < 6; k++)
	{
		g1 += final_weight[k] * p[k];
		g1x += final_weight[k] * px[k];
	}
	
	// gauss point 2
	x = 0.0;
	Polynomial_9th(p, px, df, wn4, wn3, wn2, wn1, w0, wp1, wp2, wp3, wp4, x, h);
	//-- - combination-- -
	g2 = 0.0;
	g2x = 0.0;

	for (int k = 0; k < 6; k++)
	{
		g2 += final_weight[k] * p[k];
		g2x += final_weight[k] * px[k];
	}

	// gauss point 3
	x = -0.5 - sqrt5 / 10.0;
	Polynomial_9th(p, px, df, wn4, wn3, wn2, wn1, w0, wp1, wp2, wp3, wp4, x, h);
	
	//-- - combination-- -
	g3 = 0.0;
	g3x = 0.0;

	for (int k = 0; k < 6; k++)
	{
		g3 += final_weight[k] * p[k];
		g3x += final_weight[k] * px[k];
	}

	// gauss point 4
	x = -0.5 + sqrt5 / 10.0;
	Polynomial_9th(p, px, df, wn4, wn3, wn2, wn1, w0, wp1, wp2, wp3, wp4, x, h);
	
	//-- - combination-- -
	g4 = 0.0;
	g4x = 0.0;

	for (int k = 0; k < 6; k++)
	{
		g4 += final_weight[k] * p[k];
		g4x += final_weight[k] * px[k];
	}
}

void Polynomial_9th(double* p, double* px, double* df, double wn4, double wn3, double wn2, double wn1, double w0, double wp1, double wp2, double wp3, double wp4, double x, double h)
{
	// 3rd-order polynomial
	double a3, b3, c3;
	a3 = (11 * w0 - 7 * wn1 + 2 * wn2) / 6.0;
	b3 = (2 * w0 - 3 * wn1 + wn2);
	c3 = (w0 - 2 * wn1 + wn2) / 2.0;
	p[0] = w0 + df[0] * (a3 + b3 * x + c3 * x * x - w0);
	px[0] = df[0] * (b3 + 2.0 * c3 * x) / h;

	a3 = (5 * w0 - wn1 + 2 * wp1) / 6.0;
	b3 = (-w0 + wp1);
	c3 = (-2 * w0 + wn1 + wp1) / 2.0;
	p[1] = w0 + df[1] * (a3 + b3 * x + c3 * x * x - w0);
	px[1] = df[1] * (b3 + 2.0 * c3 * x) / h;

	a3 = (2 * w0 + 5 * wp1 - wp2) / 6.0;
	b3 = (-w0 + wp1);
	c3 = (w0 - 2 * wp1 + wp2) / 2.0;
	p[2] = w0 + df[2] * (a3 + b3 * x + c3 * x * x - w0);
	px[2] = df[2] * (b3 + 2.0 * c3 * x) / h;
	// 5th-order polynomial
	double a5, b5, c5, d5, e5;
	a5 = (47 * w0 - 13 * wn1 + 2 * wn2 + 27 * wp1 - 3 * wp2) / 60.0;
	b5 = (-15 * w0 + wn1 + 15 * wp1 - wp2) / 12.0;
	c5 = (-8 * w0 + 6 * wn1 - wn2 + 2 * wp1 + wp2) / 8.0;
	d5 = (3 * w0 - wn1 - 3 * wp1 + wp2) / 6.0;
	e5 = (6 * w0 - 4 * wn1 + wn2 - 4 * wp1 + wp2) / 24.0;
	p[3] = w0 + df[3] * (a5 + b5 * x + c5 * x * x + d5 * x * x * x + e5 * x * x * x * x - w0);
	px[3] = df[3] * (b5 + 2.0 * c5 * x + 3.0 * d5 * x * x + 4.0 * e5 * x * x * x) / h;
	// 7th-order polynomial
	double a7, b7, c7, d7, e7, f7, g7;
	a7 = (319 * w0 - 101 * wn1 + 25 * wn2 - 3 * wn3 + 214 * wp1 - 38 * wp2 + 4 * wp3) / 420.0;
	b7 = (-245 * w0 + 25 * wn1 - 2 * wn2 + 245 * wp1 - 25 * wp2 + 2 * wp3) / 180.0;
	c7 = (-230 * w0 + 210 * wn1 - 57 * wn2 + 7 * wn3 + 15 * wp1 + 63 * wp2 - 8 * wp3) / 240.0;
	d7 = (28 * w0 - 11 * wn1 + wn2 - 28 * wp1 + 11 * wp2 - wp3) / 36.0;
	e7 = (46 * w0 - 39 * wn1 + 15 * wn2 - 2 * wn3 - 24 * wp1 + 3 * wp2 + wp3) / 144.0;
	f7 = (-10 * w0 + 5 * wn1 - wn2 + 10 * wp1 - 5 * wp2 + wp3) / 120.0;
	g7 = (-20 * w0 + 15 * wn1 - 6 * wn2 + wn3 + 15 * wp1 - 6 * wp2 + wp3) / 720.0;
	p[4] = w0 + df[4] * (a7 + b7 * x + c7 * x * x + d7 * x * x * x + e7 * x * x * x * x + f7 * x * x * x * x * x
		 + g7 * x * x * x * x * x * x - w0);
	px[4] = df[4] * (b7 + 2.0 * c7 * x + 3.0 * d7 * x * x + 4.0 * e7 * x * x * x + 5.0 * f7 * x * x * x * x
		 + 6.0 * g7 * x * x * x * x * x) / h;
	// 9th-order polynomial
	double a9, b9, c9, d9, e9, f9, g9, h9, i9;
	a9 = (1879 * w0 - 641 * wn1 + 199 * wn2 - 41 * wn3 + 4 * wn4 + 1375 * wp1 - 305 * wp2 + 55 * wp3 - 5 * wp4) / 2520.0;
	b9 = (-7175 * w0 + 889 * wn1 - 119 * wn2 + 9 * wn3 + 7175 * wp1 - 889 * wp2 + 119 * wp3 - 9 * wp4) / 5040.0;
	c9 = (-27895 * w0 + 28679 * wn1 - 9835 * wn2 + 2081 * wn3 - 205 * wn4 - 2065 * wp1 + 11459 * wp2 - 2455 * wp3 + 236 * wp4) / 30240.0;
	d9 = (1365 * w0 - 587 * wn1 + 89 * wn2 - 7 * wn3 - 1365 * wp1 + 587 * wp2 - 89 * wp3 + 7 * wp4) / 1440.0;
	e9 = (1174 * w0 - 1160 * wn1 + 556 * wn2 - 128 * wn3 + 13 * wn4 - 464 * wp1 - 68 * wp2 + 88 * wp3 - 11 * wp4) / 3456.0;
	f9 = (-75 * w0 + 41 * wn1 - 11 * wn2 + wn3 + 75 * wp1 - 41 * wp2 + 11 * wp3 - wp4) / 480.0;
	g9 = (-380 * w0 + 334 * wn1 - 170 * wn2 + 46 * wn3 - 5 * wn4 + 250 * wp1 - 86 * wp2 + 10 * wp3 + wp4) / 8640.0;
	h9 = (35 * w0 - 21 * wn1 + 7 * wn2 - wn3 - 35 * wp1 + 21 * wp2 - 7 * wp3 + wp4) / 5040.0;
	i9 = (70 * w0 - 56 * wn1 + 28 * wn2 - 8 * wn3 + wn4 - 56 * wp1 + 28 * wp2 - 8 * wp3 + wp4) / 40320.0;
	p[5] = w0 + df[5] * (a9 + b9 * x + c9 * x * x + d9 * x * x * x + e9 * x * x * x * x + f9 * x * x * x * x * x
		 + g9 * x * x * x * x * x * x + h9 * x * x * x * x * x * x * x + i9 * x * x * x * x * x * x * x * x - w0);
	px[5] = df[5] * (b9 + 2.0 * c9 * x + 3.0 * d9 * x * x + 4.0 * e9 * x * x * x + 5.0 * f9 * x * x * x * x
		 + 6.0 * g9 * x * x * x * x * x + 7.0 * h9 * x * x * x * x * x * x + 8.0 * i9 * x * x * x * x * x * x * x) / h;
}

// cell center rreconstruction
void Reconstruction_forg0(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{

#pragma omp parallel  for
	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			g0reconstruction_2D_normal(&xinterfaces[i * (block.ny + 1) + j], &yinterfaces[i * (block.ny + 1) + j], &fluids[i * (block.ny) + j], block);
		}
	}

	// then get the guass point value. That is, so called multi-dimensional property
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nx - block.ghost + 1; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost + 1; j++)
		{
			g0reconstruction_2D_tangent(&xinterfaces[i * (block.ny + 1) + j], &yinterfaces[i * (block.ny + 1) + j], &fluids[i * (block.ny) + j], block);
		}
	}
}

void Center_do_nothing_normal(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{
	// do nothing;
}

void Center_all_collision_multi(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{
	if (flux_function_2d == GKS2D)
	{
		for (int m = 0; m < gausspoint; m++)
		{
			Center_all_collision_2d_multi(xinterfaces[0].gauss[m]);
			Center_all_collision_2d_multi(yinterfaces[0].gauss[m]);
		}
	}
}

void Center_all_collision_2d_multi(Recon2d& gauss)
{
	double prim_left[4], prim_right[4];
	Convar_to_ULambda_2d(prim_left, gauss.left.convar);
	Convar_to_ULambda_2d(prim_right, gauss.right.convar);

	MMDF2d ml(prim_left[1], prim_left[2], prim_left[3]);
	MMDF2d mr(prim_right[1], prim_right[2], prim_right[3]);
	Collision(gauss.center.convar, prim_left[0], prim_right[0], ml, mr);

	double a0[4] = { 1.0, 0.0, 0.0, 0.0 };
	double alx[4], arx[4];
	//w_x
	A(alx, gauss.left.der1x, prim_left);
	A(arx, gauss.right.der1x, prim_right);

	double al0x[4];
	double ar0x[4];
	GL_address(0, 0, 0, al0x, alx, ml);
	GR_address(0, 0, 0, ar0x, arx, mr);

	double aly[4], ary[4];
	//w_y
	A(aly, gauss.left.der1y, prim_left);
	A(ary, gauss.right.der1y, prim_right);
	double al0y[4];
	double ar0y[4];
	GL_address(0, 0, 0, al0y, aly, ml);
	GR_address(0, 0, 0, ar0y, ary, mr);

	for (int var = 0; var < 4; var++)
	{
		gauss.center.der1x[var] = prim_left[0] * al0x[var] + prim_right[0] * ar0x[var];
		gauss.center.der1y[var] = prim_left[0] * al0y[var] + prim_right[0] * ar0y[var];
	}
}

// three-dimensional problem
G0_construct_type g0type = all_collisionn;
Reconstruction_within_Cell_3D_of_face_value cellreconstruction_3D_of_face_value = WENO5_AO_splitting_3d;
Reconstruction_forG0_3D_of_face_value g0reconstruction_3D_of_face_value = Do_nothing_splitting_3d;
Reconstruction_within_Cell_3D_of_line_value cellreconstruction_3D_of_line_value= WENO5_AO_of_line_value;
Reconstruction_within_Cell_3D_of_point_value cellreconstruction_3D_of_point_value= WENO5_AO_of_point_value;
Reconstruction_forG0_3D_of_line_value g0reconstruction_3D_of_line_value= Do_nothing_reconstruction_of_line_value;
Reconstruction_forG0_3D_of_point_value g0reconstruction_3D_of_point_value= Do_nothing_reconstruction_of_point_value;

void Reconstruction_within_cell(Interface3d *xinterfaces, Interface3d *yinterfaces, Interface3d *zinterfaces, Fluid3d *fluids, Block3d block)
{
#pragma omp parallel  for
	for (int i = block.ghost - 2; i < block.nx - block.ghost + 2; i++)
	{
		for (int j = block.ghost - 2; j < block.ny - block.ghost + 2; j++)
		{
			for (int k = block.ghost - 2; k < block.nz - block.ghost + 2; k++)
			{
				int index = i*(block.ny)*(block.nz) + j*(block.nz) + k;
				int index_face = i*(block.ny + 1)*(block.nz + 1) + j*(block.nz + 1) + k;
				int facex = index_face + (block.ny + 1)*(block.nz + 1);
				int facey = index_face + (block.nz + 1);
				int facez = index_face + 1;
				cellreconstruction_3D_of_face_value
				(xinterfaces[index_face], xinterfaces[facex],
					yinterfaces[index_face], yinterfaces[facey],
					zinterfaces[index_face], zinterfaces[facez], &fluids[index], block);
			}
		}
	}
#pragma omp parallel  for
	for (int i = block.ghost - 2; i < block.nx - block.ghost + 2; i++)
	{
		for (int j = block.ghost - 2; j < block.ny - block.ghost + 2; j++)
		{
			for (int k = block.ghost - 2; k < block.nz - block.ghost + 2; k++)
			{
				int index = i * (block.ny)*(block.nz) + j * (block.nz) + k;
				int index_face = i * (block.ny + 1)*(block.nz + 1) + j * (block.nz + 1) + k;
				int facex = index_face + (block.ny + 1)*(block.nz + 1);
				int facey = index_face + (block.nz + 1);
				int facez = index_face + 1;
				cellreconstruction_3D_of_line_value
				(&xinterfaces[facex], &yinterfaces[facey], &zinterfaces[facez], block);
			}
		}
	}
#pragma omp parallel  for
	for (int i = block.ghost - 1; i < block.nx - block.ghost + 1; i++)
	{
		for (int j = block.ghost - 1; j < block.ny - block.ghost + 1; j++)
		{
			for (int k = block.ghost - 1; k < block.nz - block.ghost + 1; k++)
			{
				int index = i * (block.ny)*(block.nz) + j * (block.nz) + k;
				int index_face = i * (block.ny + 1)*(block.nz + 1) + j * (block.nz + 1) + k;
				int facex = index_face + (block.ny + 1)*(block.nz + 1);
				int facey = index_face + (block.nz + 1);
				int facez = index_face + 1;
				cellreconstruction_3D_of_point_value
				(&xinterfaces[facex], &yinterfaces[facey], &zinterfaces[facez], block);
			}
		}
	}
}

void Check_Order_Reduce_by_Lambda_3D(bool &order_reduce, double *convar)
	{

		order_reduce = false;
		double lambda;
		lambda = Lambda(convar[0],
			convar[1] / convar[0],
			convar[2] / convar[0],
			convar[3] / convar[0],
			convar[4]);

		//if lambda <0, then reduce to the first order
		if (lambda <= 0.0 ||isnan(lambda)==true)
		{
			order_reduce = true;
		}
	}

// face reconstruction
void WENO5_AO_splitting_3d(Interface3d &left, Interface3d &right, Interface3d &down, Interface3d &up, Interface3d &back, Interface3d &front, Fluid3d *fluids, Block3d block)
{
	//DoNothing(left.line.right, fluids[0].convar);
	//DoNothing(right.line.left, fluids[0].convar);
	int index = block.ny* block.nz;
	if ((fluids[0].index[0] > block.ghost - 2) && (fluids[0].index[0] < block.nx - block.ghost + 1))
	{
		WENO5_AO(left.face.right, right.face.left,
			fluids[-2 * index].convar, fluids[-index].convar, fluids[0].convar, fluids[index].convar, fluids[2 * index].convar
			, block.dx);
	}


	index = block.nz;
	double wn2tmp[5], wn1tmp[5], wtmp[5], wp1tmp[5], wp2tmp[5];
	if ((fluids[0].index[1] > block.ghost - 2) && (fluids[0].index[1] < block.ny - block.ghost + 1))
	{
		Ydirection(wn1tmp, fluids[-index].convar); Ydirection(wtmp, fluids[0].convar); Ydirection(wp1tmp, fluids[index].convar);
		Ydirection(wn2tmp, fluids[-2 * index].convar); Ydirection(wp2tmp, fluids[2 * index].convar);
		WENO5_AO(down.face.right, up.face.left, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, block.dy);
	}
	if ((fluids[0].index[2] > block.ghost - 2) && (fluids[0].index[2] < block.nz - block.ghost + 1))
	{
		Zdirection(wn1tmp, fluids[-1].convar); Zdirection(wtmp, fluids[0].convar); Zdirection(wp1tmp, fluids[1].convar);
		Zdirection(wn2tmp, fluids[-2].convar); Zdirection(wp2tmp, fluids[2].convar);
		WENO5_AO(back.face.right, front.face.left, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, block.dz);
	}
}

void WENO5_AO(Point3d &left, Point3d &right, double * wn2, double * wn1, double * w, double * wp1, double * wp2, double h)
	{
		//we denote that   |left...cell-center...right|
		double ren2[5], ren1[5], re0[5], rep1[5], rep2[5];
		double var[5], der1[5],der2[5];

		double base_left[5];
		double base_right[5];
		double wn1_primvar[5], w_primvar[5], wp1_primvar[5];
		Convar_to_primvar_3D(wn1_primvar, wn1);
		Convar_to_primvar_3D(w_primvar, w);
		Convar_to_primvar_3D(wp1_primvar, wp1);

		for (int i = 0; i < 5; i++)
		{
			base_left[i] = 0.5*(wn1_primvar[i] + w_primvar[i]);
			base_right[i] = 0.5*(wp1_primvar[i] + w_primvar[i]);
		}

		if (reconstruction_variable == conservative)
		{
			for (int i = 0; i < 5; i++)
			{
				ren2[i] = wn2[i];
				ren1[i] = wn1[i];
				re0[i] = w[i];
				rep1[i] = wp1[i];
				rep2[i] = wp2[i];
			}
		}
		else
		{
			double s[25];
			Char_base_3D(s, base_left);
			Convar_to_char_3D(ren2, s, wn2);
			Convar_to_char_3D(ren1, s, wn1);
			Convar_to_char_3D(re0, s, w);
			Convar_to_char_3D(rep1, s, wp1);
			Convar_to_char_3D(rep2, s, wp2);
		}
		for (int i = 0; i < 5; i++)
		{
			weno_5th_ao_left(var[i],der1[i],der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h);
		}

		if (reconstruction_variable == conservative)
		{
			for (int i = 0; i < 5; i++)
			{
				left.convar[i] = var[i];
				left.der1x[i] = der1[i];
			}
		}
		else
		{
			Char_to_convar_3D(left.convar, base_left, var);
			Char_to_convar_3D(left.der1x, base_left, der1);
		}

		// cell right
		if (reconstruction_variable == conservative)
		{
			for (int i = 0; i < 5; i++)
			{
				ren2[i] = wn2[i];
				ren1[i] = wn1[i];
				re0[i] = w[i];
				rep1[i] = wp1[i];
				rep2[i] = wp2[i];
			}
		}
		else
		{
			double s[25];
			Char_base_3D(s, base_right);
			Convar_to_char_3D(ren2, s, wn2);
			Convar_to_char_3D(ren1, s, wn1);
			Convar_to_char_3D(re0, s, w);
			Convar_to_char_3D(rep1, s, wp1);
			Convar_to_char_3D(rep2, s, wp2);
		}

		for (int i = 0; i < 5; i++)
		{
			weno_5th_ao_right(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h);
		}

		if (reconstruction_variable == conservative)
		{
			for (int i = 0; i < 5; i++)
			{
				right.convar[i] = var[i];
				right.der1x[i] = der1[i];
			}
		}
		else
		{
			Char_to_convar_3D(right.convar, base_right, var);
			Char_to_convar_3D(right.der1x, base_right, der1);
		}

		bool reduce_order_left, reduce_order_right;
		Check_Order_Reduce_by_Lambda_3D(reduce_order_left, right.convar);
		Check_Order_Reduce_by_Lambda_3D(reduce_order_right, left.convar);
		if (reduce_order_left == true || reduce_order_right == true)
		{
			if (is_reduce_order_warning == true)
			{
				cout << "WENO-AO-cell-splitting order reduce" << endl;
			}	
			
			for (int m = 0; m < 5; m++)
			{
				left.convar[m] = w[m];
				right.convar[m] = w[m];
				left.der1x[m] = 0.0;
				right.der1x[m] = 0.0;
			}
		}
	}

// line reconstruction
void WENO5_AO_of_line_value(Interface3d *right, Interface3d *up, Interface3d *front, Block3d block)
{
	// for y-z plane
	int index = block.nz + 1;
	if ((right[0].index[1] > block.ghost - 1) && (right[0].index[1] < block.ny - block.ghost)&&
		(right[0].index[0] > block.ghost - 1) && (right[0].index[0] < block.nx - block.ghost + 1))
	{
		WENO5_AO_for_line_value(right[0].line, right[-2 * index].face, right[-index].face,
			right[0].face, right[index].face, right[2 * index].face, block.dy);
	}

	if ((up[0].index[0] > block.ghost - 1) && (up[0].index[0] < block.nx - block.ghost)&&
		(up[0].index[1] > block.ghost - 1) && (up[0].index[1] < block.ny - block.ghost + 1))
	{
		// x-z plane
		index = -(block.nz + 1)*(block.ny + 1);
		WENO5_AO_for_line_value(up[0].line, up[-2 * index].face, up[-index].face,
			up[0].face, up[index].face, up[2 * index].face, block.dx);
	}

	if ((front[0].index[1] > block.ghost - 1) && (front[0].index[1] < block.ny - block.ghost)&&
		(front[0].index[2] > block.ghost - 1) && (front[0].index[2] < block.nz - block.ghost + 1))
	{
		// x-y plane
		index = block.nz + 1;
		WENO5_AO_for_line_value(front[0].line, front[-2 * index].face, front[-index].face,
			front[0].face, front[index].face, front[2 * index].face, block.dy);
	}
}

void WENO5_AO_for_line_value(Recon3d *line, Recon3d& wn2, Recon3d& wn1, Recon3d& w0, Recon3d& wp1, Recon3d& wp2, double h)
{
	int sqrtgausspoint = sqrt(gausspoint);
	//we first reconstruction the left value;
	double tmp;
	double ren2[5], ren1[5], re0[5], rep1[5], rep2[5];
	double var[5], der1[5];
	double left[5], right[5];
	double base_left[5];
	double base_right[5];
	double wn1_primvar[5], w_primvar[5], wp1_primvar[5];
	Convar_to_primvar_3D(w_primvar, w0.left.convar);

	for (int i = 0; i < 5; i++)
	{
		base_left[i] = (w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{

		for (int i = 0; i < 5; i++)
		{
			if (sqrtgausspoint == 2)
			{
				weno_5th_ao_2gauss(line[0].left.convar[i], line[0].left.der1y[i], tmp,
					line[1].left.convar[i], line[1].left.der1y[i], tmp,
					wn2.left.convar[i], wn1.left.convar[i], w0.left.convar[i], wp1.left.convar[i], wp2.left.convar[i], h,1);
				weno_5th_ao_2gauss(line[0].left.der1x[i], tmp, tmp,
					line[1].left.der1x[i], tmp, tmp,
					wn2.left.der1x[i], wn1.left.der1x[i], w0.left.der1x[i], wp1.left.der1x[i], wp2.left.der1x[i], h,0);
			}
		}
	}
	else
	{
		double s[25];
		Char_base_3D(s, base_left);
		Convar_to_char_3D(ren2, s, wn2.left.convar);
		Convar_to_char_3D(ren1, s, wn1.left.convar);
		Convar_to_char_3D(re0, s, w0.left.convar);
		Convar_to_char_3D(rep1, s, wp1.left.convar);
		Convar_to_char_3D(rep2, s, wp2.left.convar);

		if (sqrtgausspoint == 2)
		{
			double var[2][5], der1[2][5];
			for (int i = 0; i < 5; i++)
			{
				weno_5th_ao_2gauss(var[0][i], der1[0][i], tmp,
					var[1][i], der1[1][i], tmp,
					ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h,1);
			}
			for (int igauss = 0; igauss < sqrtgausspoint; igauss++)
			{
				Char_to_convar_3D(line[igauss].left.convar, base_left, var[igauss]);
				Char_to_convar_3D(line[igauss].left.der1y, base_left, der1[igauss]);
			}
		}

		Convar_to_char_3D(ren2, s, wn2.left.der1x);
		Convar_to_char_3D(ren1, s, wn1.left.der1x);
		Convar_to_char_3D(re0, s, w0.left.der1x);
		Convar_to_char_3D(rep1, s, wp1.left.der1x);
		Convar_to_char_3D(rep2, s, wp2.left.der1x);

		if (sqrtgausspoint == 2)
		{
			double var[2][5], der1[2][5];
			for (int i = 0; i < 5; i++)
			{
				weno_5th_ao_2gauss(var[0][i], tmp, tmp,
					var[1][i], tmp, tmp,
					ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h,0);
			}
			for (int igauss = 0; igauss < sqrtgausspoint; igauss++)
			{
				Char_to_convar_3D(line[igauss].left.der1x, base_left, var[igauss]);
			}
		}


	}

	Convar_to_primvar_3D(w_primvar, w0.right.convar);

	for (int i = 0; i < 5; i++)
	{
		base_right[i] = (w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{

		for (int i = 0; i < 5; i++)
		{
			if (sqrtgausspoint == 2)
			{
				weno_5th_ao_2gauss(line[0].right.convar[i], line[0].right.der1y[i], tmp,
					line[1].right.convar[i], line[1].right.der1y[i], tmp,
					wn2.right.convar[i], wn1.right.convar[i], w0.right.convar[i], wp1.right.convar[i], wp2.right.convar[i], h,1);
				weno_5th_ao_2gauss(line[0].right.der1x[i], tmp, tmp,
					line[1].right.der1x[i], tmp, tmp,
					wn2.right.der1x[i], wn1.right.der1x[i], w0.right.der1x[i], wp1.right.der1x[i], wp2.right.der1x[i], h,0);
			}
		}
	}
	else
	{
		double s[25];
		Char_base_3D(s, base_right);
		Convar_to_char_3D(ren2, s, wn2.right.convar);
		Convar_to_char_3D(ren1, s, wn1.right.convar);
		Convar_to_char_3D(re0, s, w0.right.convar);
		Convar_to_char_3D(rep1, s, wp1.right.convar);
		Convar_to_char_3D(rep2, s, wp2.right.convar);

		if (sqrtgausspoint == 2)
		{
			double var[2][5], der1[2][5];
			for (int i = 0; i < 5; i++)
			{
				weno_5th_ao_2gauss(var[0][i], der1[0][i], tmp,
					var[1][i], der1[1][i], tmp,
					ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h,1);
			}
			for (int igauss = 0; igauss < sqrtgausspoint; igauss++)
			{
				Char_to_convar_3D(line[igauss].right.convar, base_right, var[igauss]);
				Char_to_convar_3D(line[igauss].right.der1y, base_right, der1[igauss]);
			}
		}

		Convar_to_char_3D(ren2, s, wn2.right.der1x);
		Convar_to_char_3D(ren1, s, wn1.right.der1x);
		Convar_to_char_3D(re0, s, w0.right.der1x);
		Convar_to_char_3D(rep1, s, wp1.right.der1x);
		Convar_to_char_3D(rep2, s, wp2.right.der1x);

		if (sqrtgausspoint == 2)
		{
			double var[2][5], der1[2][5];
			for (int i = 0; i < 5; i++)
			{
				weno_5th_ao_2gauss(var[0][i], tmp, tmp,
					var[1][i], tmp, tmp,
					ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h,0);
			}
			for (int igauss = 0; igauss < sqrtgausspoint; igauss++)
			{
				Char_to_convar_3D(line[igauss].right.der1x, base_right, var[igauss]);
			}
		}
	}

	bool reduce_order_left[3];
	bool reduce_order_right[3];
	bool reduce_order=false;

	for (int igauss = 0; igauss < sqrtgausspoint; igauss++)
	{
		Check_Order_Reduce_by_Lambda_3D(reduce_order_left[igauss], line[igauss].left.convar);
		Check_Order_Reduce_by_Lambda_3D(reduce_order_right[igauss], line[igauss].right.convar);
		if (reduce_order_left[igauss] == true|| reduce_order_right[igauss] == true)
		{
			reduce_order = true;
		}
	}

	if (reduce_order==true||line[0].order_reduce==true)
	{
		//if (is_reduce_order_warning == true)
		//{
		//	cout << "WENO-AO-face-spiltting order reduce" << endl;
		//}
		for (int igauss = 0; igauss < sqrtgausspoint; igauss++)
		{
			line[igauss].order_reduce = true;
			for (int m = 0; m < 5; m++)
			{
				line[igauss].left.convar[m] = w0.left.convar[m];
				line[igauss].right.convar[m] = w0.right.convar[m];
				line[igauss].left.der1x[m] = w0.left.der1x[m];
				line[igauss].right.der1x[m] = w0.right.der1x[m];
				line[igauss].left.der1y[m] = 0.0;
				line[igauss].right.der1y[m] = 0.0;
			}
		}
	}
}

// point reconstruction
void WENO5_AO_of_point_value(Interface3d *right, Interface3d *up, Interface3d *front, Block3d block)
{
	int sqrtgausspoint = sqrt(gausspoint);
	// for y-z plane
	for (int i = 0; i < sqrtgausspoint; i++)
	{
		WENO5_AO_for_point_value(&right[0].gauss[i*sqrtgausspoint], right[-2].line[i], right[-1].line[i],
			right[0].line[i], right[1].line[i], right[2].line[i], block.dz);
		WENO5_AO_for_point_value(&up[0].gauss[i*sqrtgausspoint], up[-2].line[i], up[-1].line[i],
			up[0].line[i], up[1].line[i], up[2].line[i], block.dz);
		int index = -(block.nz + 1)*(block.ny + 1);
		WENO5_AO_for_point_value(&front[0].gauss[i*sqrtgausspoint], front[-2 * index].line[i], front[-index].line[i],
			front[0].line[i], front[index].line[i], front[2 * index].line[i], block.dx);
	}
}

void WENO5_AO_for_point_value(Recon3d *point, Recon3d& wn2, Recon3d& wn1, Recon3d& w0, Recon3d& wp1, Recon3d& wp2, double h)
{
	int sqrtgausspoint = sqrt(gausspoint);
	//if (sqrtgausspoint != 2) { cout << "problem here" << endl; }
	//we first reconstruction the left value;
	double tmp;
	double ren2[5], ren1[5], re0[5], rep1[5], rep2[5];
	double var[5], der1[5];
	double left[5], right[5];
	double base_left[5];
	double base_right[5];
	double wn1_primvar[5], w_primvar[5], wp1_primvar[5];
	Convar_to_primvar_3D(w_primvar, w0.left.convar);
	for (int i = 0; i < 5; i++)
	{
		base_left[i] = (w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{

		for (int i = 0; i < 5; i++)
		{
			if (sqrtgausspoint == 2)
			{
				weno_5th_ao_2gauss(point[0].left.convar[i], point[0].left.der1z[i], tmp,
					point[1].left.convar[i], point[1].left.der1z[i], tmp,
					wn2.left.convar[i], wn1.left.convar[i], w0.left.convar[i], wp1.left.convar[i], wp2.left.convar[i], h,1);
				weno_5th_ao_2gauss(point[0].left.der1x[i], tmp, tmp,
					point[1].left.der1x[i], tmp, tmp,
					wn2.left.der1x[i], wn1.left.der1x[i], w0.left.der1x[i], wp1.left.der1x[i], wp2.left.der1x[i], h,0);
				weno_5th_ao_2gauss(point[0].left.der1y[i], tmp, tmp,
					point[1].left.der1y[i], tmp, tmp,
					wn2.left.der1y[i], wn1.left.der1y[i], w0.left.der1y[i], wp1.left.der1y[i], wp2.left.der1y[i], h,0);
			}
		}
	}
	else
	{
		double s[25];
		Char_base_3D(s, base_left);
		Convar_to_char_3D(ren2, s, wn2.left.convar);
		Convar_to_char_3D(ren1, s, wn1.left.convar);
		Convar_to_char_3D(re0, s, w0.left.convar);
		Convar_to_char_3D(rep1, s, wp1.left.convar);
		Convar_to_char_3D(rep2, s, wp2.left.convar);

		if (sqrtgausspoint == 2)
		{
			double var[2][5], der1[2][5];
			for (int i = 0; i < 5; i++)
			{
				weno_5th_ao_2gauss(var[0][i], der1[0][i], tmp,
					var[1][i], der1[1][i], tmp,
					ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h,1);
			}
			for (int igauss = 0; igauss < sqrtgausspoint; igauss++)
			{
				Char_to_convar_3D(point[igauss].left.convar, base_left, var[igauss]);
				Char_to_convar_3D(point[igauss].left.der1z, base_left, der1[igauss]);
			}
		}

		Convar_to_char_3D(ren2, s, wn2.left.der1x);
		Convar_to_char_3D(ren1, s, wn1.left.der1x);
		Convar_to_char_3D(re0, s, w0.left.der1x);
		Convar_to_char_3D(rep1, s, wp1.left.der1x);
		Convar_to_char_3D(rep2, s, wp2.left.der1x);

		if (sqrtgausspoint == 2)
		{
			double var[2][5], der1[2][5];
			for (int i = 0; i < 5; i++)
			{
				weno_5th_ao_2gauss(var[0][i], tmp, tmp,
					var[1][i], tmp, tmp,
					ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h,0);
			}
			for (int igauss = 0; igauss < sqrtgausspoint; igauss++)
			{
				Char_to_convar_3D(point[igauss].left.der1x, base_left, var[igauss]);
			}
		}

		Convar_to_char_3D(ren2, s, wn2.left.der1y);
		Convar_to_char_3D(ren1, s, wn1.left.der1y);
		Convar_to_char_3D(re0, s, w0.left.der1y);
		Convar_to_char_3D(rep1, s, wp1.left.der1y);
		Convar_to_char_3D(rep2, s, wp2.left.der1y);

		if (sqrtgausspoint == 2)
		{
			double var[2][5], der1[2][5];
			for (int i = 0; i < 5; i++)
			{
				weno_5th_ao_2gauss(var[0][i], tmp, tmp,
					var[1][i], tmp, tmp,
					ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h,0);
			}
			for (int igauss = 0; igauss < sqrtgausspoint; igauss++)
			{
				Char_to_convar_3D(point[igauss].left.der1y, base_left, var[igauss]);
			}
		}

	}

	Convar_to_primvar_3D(w_primvar, w0.right.convar);

	for (int i = 0; i < 5; i++)
	{
		base_right[i] = (w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{

		for (int i = 0; i < 5; i++)
		{
			if (sqrtgausspoint == 2)
			{
				weno_5th_ao_2gauss(point[0].right.convar[i], point[0].right.der1z[i], tmp,
					point[1].right.convar[i], point[1].right.der1z[i], tmp,
					wn2.right.convar[i], wn1.right.convar[i], w0.right.convar[i], wp1.right.convar[i], wp2.right.convar[i], h,1);
				weno_5th_ao_2gauss(point[0].right.der1x[i], tmp, tmp,
					point[1].right.der1x[i], tmp, tmp,
					wn2.right.der1x[i], wn1.right.der1x[i], w0.right.der1x[i], wp1.right.der1x[i], wp2.right.der1x[i], h,0);
				weno_5th_ao_2gauss(point[0].right.der1y[i], tmp, tmp,
					point[1].right.der1y[i], tmp, tmp,
					wn2.right.der1y[i], wn1.right.der1y[i], w0.right.der1y[i], wp1.right.der1y[i], wp2.right.der1y[i], h,0);
			}
		}
	}
	else
	{
		double s[25];
		Char_base_3D(s, base_right);
		Convar_to_char_3D(ren2, s, wn2.right.convar);
		Convar_to_char_3D(ren1, s, wn1.right.convar);
		Convar_to_char_3D(re0, s, w0.right.convar);
		Convar_to_char_3D(rep1, s, wp1.right.convar);
		Convar_to_char_3D(rep2, s, wp2.right.convar);

		if (sqrtgausspoint == 2)
		{
			double var[2][5], der1[2][5];
			for (int i = 0; i < 5; i++)
			{
				weno_5th_ao_2gauss(var[0][i], der1[0][i], tmp,
					var[1][i], der1[1][i], tmp,
					ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h,1);
			}
			for (int igauss = 0; igauss < sqrtgausspoint; igauss++)
			{
				Char_to_convar_3D(point[igauss].right.convar, base_right, var[igauss]);
				Char_to_convar_3D(point[igauss].right.der1z, base_right, der1[igauss]);
			}
		}

		Convar_to_char_3D(ren2, s, wn2.right.der1x);
		Convar_to_char_3D(ren1, s, wn1.right.der1x);
		Convar_to_char_3D(re0, s, w0.right.der1x);
		Convar_to_char_3D(rep1, s, wp1.right.der1x);
		Convar_to_char_3D(rep2, s, wp2.right.der1x);

		if (sqrtgausspoint == 2)
		{
			double var[2][5], der1[2][5];
			for (int i = 0; i < 5; i++)
			{
				weno_5th_ao_2gauss(var[0][i], tmp, tmp,
					var[1][i], tmp, tmp,
					ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h,0);
			}
			for (int igauss = 0; igauss < sqrtgausspoint; igauss++)
			{
				Char_to_convar_3D(point[igauss].right.der1x, base_right, var[igauss]);
			}
		}

		Convar_to_char_3D(ren2, s, wn2.right.der1y);
		Convar_to_char_3D(ren1, s, wn1.right.der1y);
		Convar_to_char_3D(re0, s, w0.right.der1y);
		Convar_to_char_3D(rep1, s, wp1.right.der1y);
		Convar_to_char_3D(rep2, s, wp2.right.der1y);

		if (sqrtgausspoint == 2)
		{
			double var[2][5], der1[2][5];
			for (int i = 0; i < 5; i++)
			{
				weno_5th_ao_2gauss(var[0][i], tmp, tmp,
					var[1][i], tmp, tmp,
					ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h,0);
			}
			for (int igauss = 0; igauss < sqrtgausspoint; igauss++)
			{
				Char_to_convar_3D(point[igauss].right.der1y, base_right, var[igauss]);
			}
		}

	}

	bool reduce_order_left[3];
	bool reduce_order_right[3];
	bool reduce_order = false;

	for (int igauss = 0; igauss < sqrtgausspoint; igauss++)
	{
		Check_Order_Reduce_by_Lambda_3D(reduce_order_left[igauss], point[igauss].left.convar);
		Check_Order_Reduce_by_Lambda_3D(reduce_order_right[igauss], point[igauss].right.convar);
		if (reduce_order_left[igauss] == true || reduce_order_right[igauss] == true)
		{
			reduce_order = true;
		}
	}

	if (reduce_order == true || point[0].order_reduce == true)
	{		
		for (int igauss = 0; igauss < sqrtgausspoint; igauss++)
		{
			point[igauss].order_reduce = true;
			for (int m = 0; m < 5; m++)
			{
				point[igauss].left.convar[m] = w0.left.convar[m];
				point[igauss].right.convar[m] = w0.right.convar[m];
				point[igauss].left.der1x[m] = w0.left.der1x[m];
				point[igauss].right.der1x[m] = w0.right.der1x[m];
				point[igauss].left.der1y[m] = w0.left.der1y[m];
				point[igauss].right.der1y[m] = w0.right.der1y[m];
				point[igauss].left.der1z[m] = 0.0;
				point[igauss].right.der1z[m] = 0.0;
			}
		}
	}
}

// center reconstruction
void Reconstruction_forg0(Interface3d *xinterfaces, Interface3d *yinterfaces, Interface3d *zinterfaces, Fluid3d *fluids, Block3d block)
{
	//go-through to obtain the avg value first
#pragma omp parallel  for
	for (int i = block.ghost - 2; i < block.nx - block.ghost + 2; i++)
	{
		for (int j = block.ghost - 2; j < block.ny - block.ghost + 2; j++)
		{
				for (int k = block.ghost - 2; k < block.nz - block.ghost + 2; k++)
				{
					int index_cell = i*(block.ny)*(block.nz) + j*(block.nz) + k;
					int index_face= i*(block.ny+1)*(block.nz + 1) + j*(block.nz + 1) + k;
					
					g0reconstruction_3D_of_face_value
					(&xinterfaces[index_face], &yinterfaces[index_face], &zinterfaces[index_face], &fluids[index_cell], block);
				}
			}
		}

#pragma omp parallel  for
	for (int i = block.ghost-2; i < block.nx - block.ghost + 2; i++)
	{
		for (int j = block.ghost-2; j < block.ny - block.ghost + 2; j++)
		{
			for (int k = block.ghost-2; k < block.nz - block.ghost + 2; k++)
			{
				int index_cell = i * (block.ny)*(block.nz) + j * (block.nz) + k;
				int index_face = i * (block.ny + 1)*(block.nz + 1) + j * (block.nz + 1) + k;

				g0reconstruction_3D_of_line_value
				(&xinterfaces[index_face], &yinterfaces[index_face], &zinterfaces[index_face], block);
			}
		}
	}
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nx - block.ghost + 1; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost + 1; j++)
		{
			for (int k = block.ghost; k < block.nz - block.ghost + 1; k++)
			{
				int index_cell = i * (block.ny)*(block.nz) + j * (block.nz) + k;
				int index_face = i * (block.ny + 1)*(block.nz + 1) + j * (block.nz + 1) + k;

				g0reconstruction_3D_of_point_value
				(&xinterfaces[index_face], &yinterfaces[index_face], &zinterfaces[index_face], block);
			}
		}
	}
}

void Do_nothing_splitting_3d(Interface3d *xinterfaces, Interface3d *yinterfaces, Interface3d *zinterfaces, 
	Fluid3d *fluids, Block3d block)
{
    // nothing need to do
}

void Do_nothing_reconstruction_of_line_value
	(Interface3d *xinterfaces, Interface3d *yinterfaces, Interface3d *zinterfaces, Block3d block)
{
	// nothing need to do
}

void Do_nothing_reconstruction_of_point_value
	(Interface3d *xinterfaces, Interface3d *yinterfaces, Interface3d *zinterfaces, Block3d block) 
{
	// nothing need to do
}





