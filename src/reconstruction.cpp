#include "reconstruction.h"

// one-dimensional problem
Reconstruction_within_Cell cellreconstruction = Vanleer; //initialization
Reconstruction_forG0 g0reconstruction = Center_collision; //initialization
Reconstruction_variable reconstruction_variable = conservative; //initialization
WENOtype wenotype = wenojs; //initialization
bool is_reduce_order_warning = false; //initialization
bool is_using_df_factor = false;

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

		if (alpha_left + alpha_right < 1.0)
		{
			fluids[i].alpha = 1.0;
		}
		else
		{
			fluids[i].alpha = 2.0 / (1.0 + (alpha_left + alpha_right) * (alpha_left + alpha_right));
		}
	}
}

void WENO5_AO_with_DF(Point1d& left, Point1d& right, Fluid1d* fluids, Block1d block)
{
	// discontiunity feedback factor
	double df[5];
	df[0] = fluids[-2].alpha;
	df[1] = fluids[-1].alpha;
	df[2] = fluids[0].alpha;
	df[3] = fluids[1].alpha;
	df[4] = fluids[2].alpha;
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

	beta[0] = 13.0 / 12.0 * pow((wn2 - 2.0 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4.0 * wn1 + 3.0 * w0), 2);
	beta[1] = 13.0 / 12.0 * pow((wn1 - 2.0 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
	beta[2] = 13.0 / 12.0 * pow((w0 - 2.0 * wp1 + wp2), 2) + 0.25 * pow((3.0 * w0 - 4.0 * wp1 + wp2), 2);

	// beta[0] = (wn2 - w0) * (wn2 - w0) + (wn1 - w0) * (wn1 - w0);
	// beta[1] = (wn1 - w0) * (wn1 - w0) + (wp1 - w0) * (wp1 - w0);
	// beta[2] = (wp2 - w0) * (wp2 - w0) + (wp1 - w0) * (wp1 - w0);
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
	double b2, c2, b4, c4, d4, e4, x = -1.0;
	b2 = wn2 - 3.0 * wn1 + 2.0 * w0;
	c2 = 0.5 * wn2 - wn1 + 0.5 * w0;
	p[0] = w0 + df[2] * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	px[0] = df[2] * (b2 + 2.0 * c2 * x) / h;
	// 3th order polynomial 2
	b2 = wp1 - w0;
	c2 = 0.5 * wn1 - w0 + 0.5 * wp1;
	p[1] = w0 + df[2] * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	px[1] = df[2] * (b2 + 2.0 * c2 * x) / h;
	// 3th order polynomial 3
	b2 = wp1 - w0;
	c2 = 0.5 * w0 - wp1 + 0.5 * wp2;
	p[2] = w0 + df[2] * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	px[2] = df[2] * (b2 + 2.0 * c2 * x) / h;
	// 5th order polynomial
	// a4 = (2 * wn2 - 13 * wn1 + 47 * w0 + 27 * wp1 - 3 * wp2) / 60;
	b4 = (wn1 - 15 * w0 + 15 * wp1 - wp2) / 12.0;
	c4 = (-wn2 + 6 * wn1 - 8 * w0 + 2 * wp1 + wp2) / 8.0;
	d4 = (-wn1 + 3 * w0 - 3 * wp1 + wp2) / 6.0;
	e4 = (wn2 - 4 * wn1 + 6 * w0 - 4 * wp1 + wp2) / 24.0;
	p[3] = w0 +
		df[2] * (b4 * (x + 1.0 / 2.0) + c4 * (x * x - 1.0 / 3.0) + d4 * (x * x * x + 1.0 / 4.0) + e4 * (x * x * x * x - 1.0 / 5.0));
	px[3] = df[2] * (b4 + 2.0 * c4 * x + 3.0 * d4 * x * x + 4.0 * e4 * x * x * x) / h;

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

	beta[0] = 13.0 / 12.0 * pow((wn2 - 2.0 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4.0 * wn1 + 3.0 * w0), 2);
	beta[1] = 13.0 / 12.0 * pow((wn1 - 2.0 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
	beta[2] = 13.0 / 12.0 * pow((w0 - 2.0 * wp1 + wp2), 2) + 0.25 * pow((3.0 * w0 - 4.0 * wp1 + wp2), 2);

	// beta[0] = (wn2 - w0) * (wn2 - w0) + (wn1 - w0) * (wn1 - w0);
	// beta[1] = (wn1 - w0) * (wn1 - w0) + (wp1 - w0) * (wp1 - w0);
	// beta[2] = (wp2 - w0) * (wp2 - w0) + (wp1 - w0) * (wp1 - w0);

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
	double b2, c2, b4, c4, d4, e4, x = 0.0;
	b2 = wn2 - 3.0 * wn1 + 2.0 * w0;
	c2 = 0.5 * wn2 - wn1 + 0.5 * w0;
	p[0] = w0 + df[2] * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	px[0] = df[2] * (b2 + 2.0 * c2 * x) / h;
	// 3th order polynomial 2
	b2 = wp1 - w0;
	c2 = 0.5 * wn1 - w0 + 0.5 * wp1;
	p[1] = w0 + df[2] * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	px[1] = df[2] * (b2 + 2.0 * c2 * x) / h;
	// 3th order polynomial 3
	b2 = wp1 - w0;
	c2 = 0.5 * w0 - wp1 + 0.5 * wp2;
	p[2] = w0 + df[2] * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	px[2] = df[2] * (b2 + 2.0 * c2 * x) / h;
	// 5th order polynomial
	// a4 = (2 * wn2 - 13 * wn1 + 47 * w0 + 27 * wp1 - 3 * wp2) / 60;
	b4 = (wn1 - 15 * w0 + 15 * wp1 - wp2) / 12.0;
	c4 = (-wn2 + 6 * wn1 - 8 * w0 + 2 * wp1 + wp2) / 8.0;
	d4 = (-wn1 + 3 * w0 - 3 * wp1 + wp2) / 6.0;
	e4 = (wn2 - 4 * wn1 + 6 * w0 - 4 * wp1 + wp2) / 24.0;
	p[3] = w0 +
		df[2] * (b4 * (x + 1.0 / 2.0) + c4 * (x * x - 1.0 / 3.0) + d4 * (x * x * x + 1.0 / 4.0) + e4 * (x * x * x * x - 1.0 / 5.0));
	px[3] = df[2] * (b4 + 2.0 * c4 * x + 3.0 * d4 * x * x + 4.0 * e4 * x * x * x) / h;

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

void WENO5_AO_with_single_weight(Point1d& left, Point1d& right, Fluid1d* fluids, Block1d block)
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
	double beta[4];
	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 3; i++)
		{
			beta[0] = 13.0 / 12.0 * pow((wn2[2] - 2.0 * wn1[2] + w0[2]), 2) + 0.25 * pow((wn2[2] - 4.0 * wn1[2] + 3.0 * w0[2]), 2);
			beta[1] = 13.0 / 12.0 * pow((wn1[2] - 2.0 * w0[2] + wp1[2]), 2) + 0.25 * pow((wn1[2] - wp1[2]), 2);
			beta[2] = 13.0 / 12.0 * pow((w0[2] - 2.0 * wp1[2] + wp2[2]), 2) + 0.25 * pow((3.0 * w0[2] - 4.0 * wp1[2] + wp2[2]), 2);
			beta[3] = (1.0 / 5040.0) * (231153.0 * w0[2] * w0[2] + 104963.0 * wn1[2] * wn1[2] + 6908.0 * wn2[2] * wn2[2] -
			38947.0 * wn2[2] * wp1[2] + 104963.0 * wp1[2] * wp1[2] +
			wn1[2] * (-51001.0 * wn2[2] + 179098.0 * wp1[2] - 38947.0 * wp2[2]) -
			3.0 * w0[2] * (99692.0 * wn1[2] - 22641.0 * wn2[2] + 99692.0 * wp1[2] - 22641.0 * wp2[2]) +
			8209.0 * wn2[2] * wp2[2] - 51001.0 * wp1[2] * wp2[2] + 6908.0 * wp2[2] * wp2[2]);
			weno_5th_ao_with_single_weight_right(right.convar[i], right.der1[i], tmp, wn2[i], wn1[i], w0[i], wp1[i], wp2[i], beta, h);
			weno_5th_ao_with_single_weight_left(left.convar[i], left.der1[i], tmp, wn2[i], wn1[i], w0[i], wp1[i], wp2[i], beta, h);
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
		beta[0] = 13.0 / 12.0 * pow((ren2[2] - 2.0 * ren1[2] + re0[2]), 2) + 0.25 * pow((ren2[2] - 4.0 * ren1[2] + 3.0 * re0[2]), 2);
		beta[1] = 13.0 / 12.0 * pow((ren1[2] - 2.0 * re0[2] + rep1[2]), 2) + 0.25 * pow((ren1[2] - rep1[2]), 2);
		beta[2] = 13.0 / 12.0 * pow((re0[2] - 2.0 * rep1[2] + rep2[2]), 2) + 0.25 * pow((3.0 * re0[2] - 4.0 * rep1[2] + rep2[2]), 2);
		beta[3] = (1.0 / 5040.0) * (231153.0 * re0[2] * re0[2] + 104963.0 * ren1[2] * ren1[2] + 6908.0 * ren2[2] * ren2[2] -
			38947.0 * ren2[2] * rep1[2] + 104963.0 * rep1[2] * rep1[2] +
			ren1[2] * (-51001.0 * ren2[2] + 179098.0 * rep1[2] - 38947.0 * rep2[2]) -
			3.0 * re0[2] * (99692.0 * ren1[2] - 22641.0 * ren2[2] + 99692.0 * rep1[2] - 22641.0 * rep2[2]) +
			8209.0 * ren2[2] * rep2[2] - 51001.0 * rep1[2] * rep2[2] + 6908.0 * rep2[2] * rep2[2]);
		for (int i = 0; i < 3; i++)
		{
			weno_5th_ao_with_single_weight_left(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], beta, h);
		}
		Char_to_convar1D(left.convar, base_left, var);
		Char_to_convar1D(left.der1, base_left, der1);
	
		// right reconstruction

		Convar_to_char1D(ren2, base_right, wn2);
		Convar_to_char1D(ren1, base_right, wn1);
		Convar_to_char1D(re0, base_right, w0);
		Convar_to_char1D(rep1, base_right, wp1);
		Convar_to_char1D(rep2, base_right, wp2);

		beta[0] = 13.0 / 12.0 * pow((ren2[2] - 2.0 * ren1[2] + re0[2]), 2) + 0.25 * pow((ren2[2] - 4.0 * ren1[2] + 3.0 * re0[2]), 2);
		beta[1] = 13.0 / 12.0 * pow((ren1[2] - 2.0 * re0[2] + rep1[2]), 2) + 0.25 * pow((ren1[2] - rep1[2]), 2);
		beta[2] = 13.0 / 12.0 * pow((re0[2] - 2.0 * rep1[2] + rep2[2]), 2) + 0.25 * pow((3.0 * re0[2] - 4.0 * rep1[2] + rep2[2]), 2);
		beta[3] = (1.0 / 5040.0) * (231153.0 * re0[2] * re0[2] + 104963.0 * ren1[2] * ren1[2] + 6908.0 * ren2[2] * ren2[2] -
			38947.0 * ren2[2] * rep1[2] + 104963.0 * rep1[2] * rep1[2] +
			ren1[2] * (-51001.0 * ren2[2] + 179098.0 * rep1[2] - 38947.0 * rep2[2]) -
			3.0 * re0[2] * (99692.0 * ren1[2] - 22641.0 * ren2[2] + 99692.0 * rep1[2] - 22641.0 * rep2[2]) +
			8209.0 * ren2[2] * rep2[2] - 51001.0 * rep1[2] * rep2[2] + 6908.0 * rep2[2] * rep2[2]);
		for (int i = 0; i < 3; i++)
		{
			weno_5th_ao_with_single_weight_right(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], beta, h);
		}
		Char_to_convar1D(right.convar, base_right, var);
		Char_to_convar1D(right.der1, base_right, der1);

	}

	Check_Order_Reduce(left, right, fluids[0]);

}

void weno_5th_ao_with_single_weight_left(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double* beta, double h)
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

void weno_5th_ao_with_single_weight_right(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double* beta, double h)
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
	//A = A + abs(prim_left[0] - prim_right[0]) / prim_left[0] + abs(prim_left[0] - prim_right[0]) / prim_right[0];
	A = A + pow(Ma_left_normal - Ma_right_normal, 2) + pow(Ma_left_tangent - Ma_right_tangent, 2);
	return A;
}

void Calculate_alpha(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d& fluids)
{
	double prim_left_left[4], prim_left_right[4], prim_right_left[4], prim_right_right[4];
	double prim_down_left[4], prim_down_right[4], prim_up_left[4], prim_up_right[4];
	double alpha_left, alpha_right, alpha_down, alpha_up;
	double alpha[2][4];
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
	alpha_left = 0.5 * (alpha[0][0] + alpha[1][0]);
	alpha_right = 0.5 * (alpha[0][1] + alpha[1][1]);
	alpha_down = 0.5 * (alpha[0][2] + alpha[1][2]);
	alpha_up = 0.5 * (alpha[0][3] + alpha[1][3]);
	double sum_alpha_x, sum_alpha_y;
	sum_alpha_x = alpha_left + alpha_right;
	sum_alpha_y = alpha_down + alpha_up;
	if (sum_alpha_x < 1.0)
	{
		fluids.alpha_x = 1.0;
	}
	else
	{
		fluids.alpha_x = 2.0 / (1.0 + sum_alpha_x * sum_alpha_x);
	}
	if (sum_alpha_y < 1.0)
	{
		fluids.alpha_y = 1.0;
	}
	else
	{
		fluids.alpha_y = 2.0 / (1.0 + sum_alpha_y * sum_alpha_y);
	}
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
		weno_5th_ao_with_df_left(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], alpha, h);

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
		weno_5th_ao_with_df_right(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], alpha, h);
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

	//double sqrt3 = sqrt(3);
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
	alpha1[0] = fluids[-2].alpha_x; alpha1[1] = fluids[-1].alpha_x;
	alpha1[2] = fluids[0].alpha_x;  alpha1[3] = fluids[1].alpha_x;
	alpha1[4] = fluids[2].alpha_x;
	alpha2[0] = fluids[block.ny - 2].alpha_x; alpha2[1] = fluids[block.ny - 1].alpha_x;
	alpha2[2] = fluids[block.ny].alpha_x;  alpha2[3] = fluids[block.ny + 1].alpha_x;
	alpha2[4] = fluids[block.ny + 2].alpha_x;
	weno_5th_ao_with_df_tangential(right[0].gauss,
		right[-2].line, right[-1].line, right[0].line, right[1].line, right[2].line, alpha1, alpha2, right[0].length);

	//since we already do the coordinate transform, along y, no transform needed.
	alpha1[0] = fluids[2 * block.ny].alpha_y; alpha1[1] = fluids[block.ny].alpha_y;
	alpha1[2] = fluids[0].alpha_y;  alpha1[3] = fluids[-block.ny].alpha_y;
	alpha1[4] = fluids[-2 * block.ny].alpha_y;
	alpha2[0] = fluids[2 * block.ny + 1].alpha_y; alpha2[1] = fluids[block.ny + 1].alpha_y;
	alpha2[2] = fluids[1].alpha_y;  alpha2[3] = fluids[-block.ny + 1].alpha_y;
	alpha2[4] = fluids[-2 * block.ny + 1].alpha_y;
	weno_5th_ao_with_df_tangential(up[0].gauss, up[2 * (block.ny + 1)].line, up[block.ny + 1].line,
		up[0].line, up[-(block.ny + 1)].line, up[-2 * (block.ny + 1)].line, alpha1, alpha2, up[0].length);
}

void weno_5th_ao_with_df_tangential(Recon2d* re, Recon2d& wn2, Recon2d& wn1, Recon2d& w0, Recon2d& wp1, Recon2d& wp2, double* alpha1, double* alpha2, double h)
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
				weno_5th_ao_with_df_2gauss(re[0].left.convar[i], re[0].left.der1y[i], tmp[0],
					re[1].left.convar[i], re[1].left.der1y[i], tmp[1], alpha1,
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
					var[1][i], der1[1][i], der2[1][i], alpha1,
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
				re[1].left.der1x[i], tmp[2], tmp[3], alpha1,
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
					re[1].right.convar[i], re[1].right.der1y[i], tmp[1], alpha2,
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
					var[1][i], der1[1][i], der2[1][i], alpha2,
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
				re[1].right.der1x[i], tmp[2], tmp[3], alpha2,
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

void weno_5th_ao_with_df_2gauss(double& g1, double& g1x, double& g1xx, double& g2, double& g2x, double& g2xx, double* alpham, double wn2, double wn1, double w0, double wp1, double wp2, double h, int order)
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

	beta[0] = 13.0 / 12.0 * pow((wn2 - 2.0 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4.0 * wn1 + 3.0 * w0), 2);
	beta[1] = 13.0 / 12.0 * pow((wn1 - 2.0 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
	beta[2] = 13.0 / 12.0 * pow((w0 - 2.0 * wp1 + wp2), 2) + 0.25 * pow((3.0 * w0 - 4.0 * wp1 + wp2), 2);
	//beta[3] = 0.25 * ((wn2 - w0) * (wn2 - w0) + (wn1 - w0) * (wn1 - w0) + (wp1 - w0) * (wp1 - w0) + (wp2 - w0) * (wp2 - w0));
	beta[3] = 1.0 / 6.0 * (beta[0] + 4.0 * beta[1] + beta[2]) + abs(beta[0] - beta[2]);

	double tau5 = (abs(beta[3] - beta[0]) + abs(beta[3] - beta[1]) + abs(beta[3] - beta[2])) / 3.0;

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
	double b2, c2, b4, c4, d4, e4, x;
	double DF31, DF32, DF33, DF5;
	DF31 = alpham[2];
	DF32 = alpham[2];
	DF33 = alpham[2];
	DF5 = alpham[2];
	// gauss point 1
	x = -0.5 - sqrt(3) / 6.0;
	// 3th order polynomial 1
	b2 = wn2 - 3.0 * wn1 + 2.0 * w0;
	c2 = 0.5 * wn2 - wn1 + 0.5 * w0;
	p[0] = w0 + DF31 * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	px[0] = DF31 * (b2 + 2.0 * c2 * x) / h;
	// 3th order polynomial 2
	b2 = wp1 - w0;
	c2 = 0.5 * wn1 - w0 + 0.5 * wp1;
	p[1] = w0 + DF32 * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	px[1] = DF32 * (b2 + 2.0 * c2 * x) / h;
	// 3th order polynomial 3
	b2 = wp1 - w0;
	c2 = 0.5 * w0 - wp1 + 0.5 * wp2;
	p[2] = w0 + DF33 * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	px[2] = DF33 * (b2 + 2.0 * c2 * x) / h;
	// 5th order polynomial
	// a4 = (2 * wn2 - 13 * wn1 + 47 * w0 + 27 * wp1 - 3 * wp2) / 60;
	b4 = (wn1 - 15 * w0 + 15 * wp1 - wp2) / 12.0;
	c4 = (-wn2 + 6 * wn1 - 8 * w0 + 2 * wp1 + wp2) / 8.0;
	d4 = (-wn1 + 3 * w0 - 3 * wp1 + wp2) / 6.0;
	e4 = (wn2 - 4 * wn1 + 6 * w0 - 4 * wp1 + wp2) / 24.0;
	p[3] = w0 +
		DF5 * (b4 * (x + 1.0 / 2.0) + c4 * (x * x - 1.0 / 3.0) + d4 * (x * x * x + 1.0 / 4.0) + e4 * (x * x * x * x - 1.0 / 5.0));
	px[3] = DF5 * (b4 + 2.0 * c4 * x + 3.0 * d4 * x * x + 4.0 * e4 * x * x * x) / h;
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
	x = -0.5 + sqrt(3) / 6.0;
	// 3th order polynomial 1
	b2 = wn2 - 3.0 * wn1 + 2.0 * w0;
	c2 = 0.5 * wn2 - wn1 + 0.5 * w0;
	p[0] = w0 + DF31 * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	px[0] = DF31 * (b2 + 2.0 * c2 * x) / h;
	// 3th order polynomial 2
	b2 = wp1 - w0;
	c2 = 0.5 * wn1 - w0 + 0.5 * wp1;
	p[1] = w0 + DF32 * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	px[1] = DF32 * (b2 + 2.0 * c2 * x) / h;
	// 3th order polynomial 3
	b2 = wp1 - w0;
	c2 = 0.5 * w0 - wp1 + 0.5 * wp2;
	p[2] = w0 + DF33 * (b2 * (x + 0.5) + c2 * (x * x - 1.0 / 3.0));
	px[2] = DF33 * (b2 + 2.0 * c2 * x) / h;
	// 5th order polynomial
	// a4 = (2 * wn2 - 13 * wn1 + 47 * w0 + 27 * wp1 - 3 * wp2) / 60;
	b4 = (wn1 - 15 * w0 + 15 * wp1 - wp2) / 12.0;
	c4 = (-wn2 + 6 * wn1 - 8 * w0 + 2 * wp1 + wp2) / 8.0;
	d4 = (-wn1 + 3 * w0 - 3 * wp1 + wp2) / 6.0;
	e4 = (wn2 - 4 * wn1 + 6 * w0 - 4 * wp1 + wp2) / 24.0;
	p[3] = w0 +
		DF5 * (b4 * (x + 1.0 / 2.0) + c4 * (x * x - 1.0 / 3.0) + d4 * (x * x * x + 1.0 / 4.0) + e4 * (x * x * x * x - 1.0 / 5.0));
	px[3] = DF5 * (b4 + 2.0 * c4 * x + 3.0 * d4 * x * x + 4.0 * e4 * x * x * x) / h;
	//-- - combination-- -
	g2 = 0.0;
	g2x = 0.0;

	for (int k = 0; k < 4; k++)
	{
		g2 += final_weight[k] * p[k];
		g2x += final_weight[k] * px[k];
	}
}

// cell center rreconstructioon
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





