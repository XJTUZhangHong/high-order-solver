#pragma once
#include "fluid_mesh.h"
#include <assert.h>

extern int K;
extern double Gamma;
extern double Mu;
extern double Nu;
extern double T_inf; //the real temperature
extern double c1_euler;
extern double c2_euler;
extern bool is_Prandtl_fix;
extern double Pr;
extern double R_gas; //the gas constant for real ideal air
enum TAU_TYPE { Euler, NS};
extern TAU_TYPE tau_type;

//basic gks function
// to store the moment
class MMDF1d
{
private:
	double u;
	double lambda;

public:
	double uwhole[10];
	double uplus[10];
	double uminus[10];
	double upxi[10][4];
	double unxi[10][4];
	double uxi[10][4];
	double xi2;
	double xi4;
	double xi6;
	MMDF1d();
	MMDF1d(double u_in, double lambda_in);
	void calcualte_MMDF1d();
};

class MMDF2d
{
private:
	double u;
	double v;
	double lambda;

public:
	double uwhole[7];
	double uplus[7];
	double uminus[7];
	double vwhole[7];
	double upvxi[7][7][3];
	double unvxi[7][7][3];
	double uvxi[7][7][3];
	double xi2;
	double xi4;
	MMDF2d();
	MMDF2d(double u_in, double v_in, double lambda_in);
	void calcualte_MMDF();
};

class MMDF3d_left
{
private:
	double u;
	double v;
	double w;
	double lambda;

public:
	double uwhole[7];
	double uplus[7];
	double vwhole[6];
	double wwhole[6];
	double upvwxi[7][6][6][3];
	double uvwxi[7][6][6][3];
	double xi[3];
	double xi2;
	double xi4;
	MMDF3d_left();
	MMDF3d_left(double u_in, double v_in, double w_in, double lambda_in);
	void calcualte_MMDF();
};

class MMDF3d_right
{
private:
	double u;
	double v;
	double w;
	double lambda;

public:
	double uwhole[7];
	double uminus[7];
	double vwhole[6];
	double wwhole[6];
	double unvwxi[7][6][6][3];
	double uvwxi[7][6][6][3];
	double xi[3];
	double xi2;
	double xi4;
	MMDF3d_right();
	MMDF3d_right(double u_in, double v_in, double w_in, double lambda_in);
	void calcualte_MMDF();
};

class MMDF3d
{
private:
	double u;
	double v;
	double w;
	double lambda;

public:
	double uwhole[7];
	double uplus[7];
	double uminus[7];
	double vwhole[6];
	double wwhole[6];
	double upvwxi[7][6][6][3];
	double unvwxi[7][6][6][3];
	double uvwxi[7][6][6][3];

	double xi2;
	double xi4;
	MMDF3d();
	MMDF3d(double u_in, double v_in, double w_in, double lambda_in);
	void calcualte_MMDF();
};

class MMDF3d_left_speed
{
private:
	double u;
	double v;
	double w;
	double lambda;

public:
	double uwhole[7];
	double uplus[7];
	double vwhole[6];
	double wwhole[6];
	double xi[3];

	MMDF3d_left_speed();
	MMDF3d_left_speed(double u_in, double v_in, double w_in, double lambda_in);
	void calcualte_MMDF();
};

class MMDF3d_right_speed
{
private:
	double u;
	double v;
	double w;
	double lambda;

public:
	double uwhole[7];
	double uminus[7];
	double vwhole[6];
	double wwhole[6];
	double xi[3];

	MMDF3d_right_speed();
	MMDF3d_right_speed(double u_in, double v_in, double w_in, double lambda_in);
	void calcualte_MMDF();
};

class MMDF3d_speed
{
private:
	double u;
	double v;
	double w;
	double lambda;

public:
	double uwhole[7];
	double uplus[7] ;
	double uminus[7];
	double vwhole[6];
	double wwhole[6];
	double xi[3];

	MMDF3d_speed();
	MMDF3d_speed(double u_in, double v_in, double w_in, double lambda_in);
	void calcualte_MMDF();
};
// to calculate the microsolpe moment
void G(int no_u, int no_xi, double* psi, double a[3], MMDF1d m);
void GL(int no_u, int no_xi, double* psi, double a[3], MMDF1d m);
void GR(int no_u, int no_xi, double* psi, double a[3], MMDF1d m);

void Microslope(double* a, double der[3], double prim[3]);

void Collision(double* w0, double left, double right, MMDF2d& m2, MMDF2d& m3);

template <class T1, class T2>
void Collision(double *w0, double left, double right, T1& m2, T2& m3)
{
	w0[0] = left * m2.uplus[0] + right * m3.uminus[0];
	w0[1] = left * m2.uplus[1] + right * m3.uminus[1];
	w0[2] = left * m2.vwhole[1] * m2.uplus[0] + right * m3.vwhole[1] * m3.uminus[0];
	w0[3] = left * m2.wwhole[1] * m2.uplus[0] + right * m3.wwhole[1] * m3.uminus[0];
	w0[4] = 0.5*left * (m2.uplus[2] + m2.uplus[0] * m2.vwhole[2] + m2.uplus[0] * m2.wwhole[2] + m2.uplus[0] * m2.xi[1]) +
		0.5*right * (m3.uminus[2] + m3.uminus[0] * m3.vwhole[2] + m3.uminus[0] * m3.wwhole[2] + m3.uminus[0] * m3.xi[1]);
}

void A(double* a, double der[4], double prim[4]);

void A(double *a, double der[5], double density, double u, double v,double w, double lambda);

void GL_address(int no_u, int no_v, int no_xi, double* psi, double a[4], MMDF2d& m);

void GR_address(int no_u, int no_v, int no_xi, double* psi, double a[4], MMDF2d& m);

template <class T>
void GL_address(int no_u, int no_v, int no_w, int no_xi, double *psi, double* a, T& m)
{
	psi[0] = a[0] * m.upvwxi[no_u][no_v][no_w][no_xi]
		+ a[1] * m.upvwxi[no_u + 1][no_v][no_w][no_xi]
		+ a[2] * m.upvwxi[no_u][no_v + 1][no_w][no_xi]
		+ a[3] * m.upvwxi[no_u][no_v][no_w + 1][no_xi]
		+ a[4] * 0.5*(m.upvwxi[no_u + 2][no_v][no_w][no_xi] + m.upvwxi[no_u][no_v + 2][no_w][no_xi] + m.upvwxi[no_u][no_v][no_w + 2][no_xi] + m.upvwxi[no_u][no_v][no_w][no_xi + 1]);

	psi[1] = a[0] * m.upvwxi[no_u + 1][no_v][no_w][no_xi]
		+ a[1] * m.upvwxi[no_u + 2][no_v][no_w][no_xi]
		+ a[2] * m.upvwxi[no_u + 1][no_v + 1][no_w][no_xi]
		+ a[3] * m.upvwxi[no_u + 1][no_v][no_w + 1][no_xi]
		+ a[4] * 0.5*(m.upvwxi[no_u + 3][no_v][no_w][no_xi] + m.upvwxi[no_u + 1][no_v + 2][no_w][no_xi] + m.upvwxi[no_u + 1][no_v][no_w + 2][no_xi] + m.upvwxi[no_u + 1][no_v][no_w][no_xi + 1]);

	psi[2] = a[0] * m.upvwxi[no_u][no_v + 1][no_w][no_xi]
		+ a[1] * m.upvwxi[no_u + 1][no_v + 1][no_w][no_xi]
		+ a[2] * m.upvwxi[no_u][no_v + 2][no_w][no_xi]
		+ a[3] * m.upvwxi[no_u][no_v + 1][no_w + 1][no_xi]
		+ a[4] * 0.5*(m.upvwxi[no_u + 2][no_v + 1][no_w][no_xi] + m.upvwxi[no_u][no_v + 3][no_w][no_xi] + m.upvwxi[no_u][no_v + 1][no_w + 2][no_xi] + m.upvwxi[no_u][no_v + 1][no_w][no_xi + 1]);

	psi[3] = a[0] * m.upvwxi[no_u][no_v][no_w + 1][no_xi]
		+ a[1] * m.upvwxi[no_u + 1][no_v][no_w + 1][no_xi]
		+ a[2] * m.upvwxi[no_u][no_v + 1][no_w + 1][no_xi]
		+ a[3] * m.upvwxi[no_u][no_v][no_w + 2][no_xi]
		+ a[4] * 0.5*(m.upvwxi[no_u + 2][no_v][no_w + 1][no_xi] + m.upvwxi[no_u][no_v + 2][no_w + 1][no_xi] + m.upvwxi[no_u][no_v][no_w + 3][no_xi] + m.upvwxi[no_u][no_v][no_w + 1][no_xi + 1]);

	psi[4] = 0.5*(
		a[0] * (m.upvwxi[no_u + 2][no_v][no_w][no_xi] + m.upvwxi[no_u][no_v + 2][no_w][no_xi] + m.upvwxi[no_u][no_v][no_w + 2][no_xi] + m.upvwxi[no_u][no_v][no_w][no_xi + 1]) +
		a[1] * (m.upvwxi[no_u + 3][no_v][no_w][no_xi] + m.upvwxi[no_u + 1][no_v + 2][no_w][no_xi] + m.upvwxi[no_u + 1][no_v][no_w + 2][no_xi] + m.upvwxi[no_u + 1][no_v][no_w][no_xi + 1]) +
		a[2] * (m.upvwxi[no_u + 2][no_v + 1][no_w][no_xi] + m.upvwxi[no_u][no_v + 3][no_w][no_xi] + m.upvwxi[no_u][no_v + 1][no_w + 2][no_xi] + m.upvwxi[no_u][no_v + 1][no_w][no_xi + 1]) +
		a[3] * (m.upvwxi[no_u + 2][no_v][no_w + 1][no_xi] + m.upvwxi[no_u][no_v + 2][no_w + 1][no_xi] + m.upvwxi[no_u][no_v][no_w + 3][no_xi] + m.upvwxi[no_u][no_v][no_w + 1][no_xi + 1]) +
		a[4] * 0.5*(m.upvwxi[no_u + 4][no_v][no_w][no_xi] +
			m.upvwxi[no_u][no_v + 4][no_w][no_xi] +
			m.upvwxi[no_u][no_v][no_w + 4][no_xi] +
			m.upvwxi[no_u][no_v][no_w][no_xi + 2] +
			2 * m.upvwxi[no_u + 2][no_v + 2][no_w][no_xi] +
			2 * m.upvwxi[no_u + 2][no_v][no_w + 2][no_xi] +
			2 * m.upvwxi[no_u + 2][no_v][no_w][no_xi + 1] +
			2 * m.upvwxi[no_u][no_v + 2][no_w + 2][no_xi] +
			2 * m.upvwxi[no_u][no_v + 2][no_w][no_xi + 1] +
			2 * m.upvwxi[no_u][no_v][no_w + 2][no_xi + 1])
		);
}

template <class T>
void GR_address(int no_u, int no_v, int no_w, int no_xi, double *psi, double *a, T &m)
{
	psi[0] = a[0] * m.unvwxi[no_u][no_v][no_w][no_xi]
		+ a[1] * m.unvwxi[no_u + 1][no_v][no_w][no_xi]
		+ a[2] * m.unvwxi[no_u][no_v + 1][no_w][no_xi]
		+ a[3] * m.unvwxi[no_u][no_v][no_w + 1][no_xi]
		+ a[4] * 0.5*(m.unvwxi[no_u + 2][no_v][no_w][no_xi] + m.unvwxi[no_u][no_v + 2][no_w][no_xi] + m.unvwxi[no_u][no_v][no_w + 2][no_xi] + m.unvwxi[no_u][no_v][no_w][no_xi + 1]);

	psi[1] = a[0] * m.unvwxi[no_u + 1][no_v][no_w][no_xi]
		+ a[1] * m.unvwxi[no_u + 2][no_v][no_w][no_xi]
		+ a[2] * m.unvwxi[no_u + 1][no_v + 1][no_w][no_xi]
		+ a[3] * m.unvwxi[no_u + 1][no_v][no_w + 1][no_xi]
		+ a[4] * 0.5*(m.unvwxi[no_u + 3][no_v][no_w][no_xi] + m.unvwxi[no_u + 1][no_v + 2][no_w][no_xi] + m.unvwxi[no_u + 1][no_v][no_w + 2][no_xi] + m.unvwxi[no_u + 1][no_v][no_w][no_xi + 1]);

	psi[2] = a[0] * m.unvwxi[no_u][no_v + 1][no_w][no_xi]
		+ a[1] * m.unvwxi[no_u + 1][no_v + 1][no_w][no_xi]
		+ a[2] * m.unvwxi[no_u][no_v + 2][no_w][no_xi]
		+ a[3] * m.unvwxi[no_u][no_v + 1][no_w + 1][no_xi]
		+ a[4] * 0.5*(m.unvwxi[no_u + 2][no_v + 1][no_w][no_xi] + m.unvwxi[no_u][no_v + 3][no_w][no_xi] + m.unvwxi[no_u][no_v + 1][no_w + 2][no_xi] + m.unvwxi[no_u][no_v + 1][no_w][no_xi + 1]);

	psi[3] = a[0] * m.unvwxi[no_u][no_v][no_w + 1][no_xi]
		+ a[1] * m.unvwxi[no_u + 1][no_v][no_w + 1][no_xi]
		+ a[2] * m.unvwxi[no_u][no_v + 1][no_w + 1][no_xi]
		+ a[3] * m.unvwxi[no_u][no_v][no_w + 2][no_xi]
		+ a[4] * 0.5*(m.unvwxi[no_u + 2][no_v][no_w + 1][no_xi] + m.unvwxi[no_u][no_v + 2][no_w + 1][no_xi] + m.unvwxi[no_u][no_v][no_w + 3][no_xi] + m.unvwxi[no_u][no_v][no_w + 1][no_xi + 1]);

	psi[4] = 0.5*(
		a[0] * (m.unvwxi[no_u + 2][no_v][no_w][no_xi] + m.unvwxi[no_u][no_v + 2][no_w][no_xi] + m.unvwxi[no_u][no_v][no_w + 2][no_xi] + m.unvwxi[no_u][no_v][no_w][no_xi + 1]) +
		a[1] * (m.unvwxi[no_u + 3][no_v][no_w][no_xi] + m.unvwxi[no_u + 1][no_v + 2][no_w][no_xi] + m.unvwxi[no_u + 1][no_v][no_w + 2][no_xi] + m.unvwxi[no_u + 1][no_v][no_w][no_xi + 1]) +
		a[2] * (m.unvwxi[no_u + 2][no_v + 1][no_w][no_xi] + m.unvwxi[no_u][no_v + 3][no_w][no_xi] + m.unvwxi[no_u][no_v + 1][no_w + 2][no_xi] + m.unvwxi[no_u][no_v + 1][no_w][no_xi + 1]) +
		a[3] * (m.unvwxi[no_u + 2][no_v][no_w + 1][no_xi] + m.unvwxi[no_u][no_v + 2][no_w + 1][no_xi] + m.unvwxi[no_u][no_v][no_w + 3][no_xi] + m.unvwxi[no_u][no_v][no_w + 1][no_xi + 1]) +
		a[4] * 0.5*(m.unvwxi[no_u + 4][no_v][no_w][no_xi] +
			m.unvwxi[no_u][no_v + 4][no_w][no_xi] +
			m.unvwxi[no_u][no_v][no_w + 4][no_xi] +
			m.unvwxi[no_u][no_v][no_w][no_xi + 2] +
			2 * m.unvwxi[no_u + 2][no_v + 2][no_w][no_xi] +
			2 * m.unvwxi[no_u + 2][no_v][no_w + 2][no_xi] +
			2 * m.unvwxi[no_u + 2][no_v][no_w][no_xi + 1] +
			2 * m.unvwxi[no_u][no_v + 2][no_w + 2][no_xi] +
			2 * m.unvwxi[no_u][no_v + 2][no_w][no_xi + 1] +
			2 * m.unvwxi[no_u][no_v][no_w + 2][no_xi + 1])
		);
}

void G_address(int no_u, int no_v, int no_xi, double* psi, double a[4], MMDF2d& m);

template <class T>
void G_address(int no_u, int no_v, int no_w,int no_xi, double *psi, double a[5], T &m)
{
	
	psi[0] = a[0] * m.uvwxi[no_u][no_v][no_w][no_xi] 
		+ a[1] * m.uvwxi[no_u + 1][no_v][no_w][no_xi] 
		+ a[2] * m.uvwxi[no_u][no_v + 1][no_w][no_xi]
		+ a[3] * m.uvwxi[no_u][no_v][no_w+1][no_xi]
		+ a[4] * 0.5*(m.uvwxi[no_u + 2][no_v][no_w][no_xi] + m.uvwxi[no_u][no_v + 2][no_w][no_xi] + m.uvwxi[no_u][no_v][no_w+2][no_xi] + m.uvwxi[no_u][no_v][no_w][no_xi + 1]);
	
	psi[1] = a[0] * m.uvwxi[no_u+1][no_v][no_w][no_xi]
		+ a[1] * m.uvwxi[no_u + 2][no_v][no_w][no_xi]
		+ a[2] * m.uvwxi[no_u+1][no_v + 1][no_w][no_xi]
		+ a[3] * m.uvwxi[no_u+1][no_v][no_w + 1][no_xi]
		+ a[4] * 0.5*(m.uvwxi[no_u + 3][no_v][no_w][no_xi] + m.uvwxi[no_u+1][no_v + 2][no_w][no_xi] + m.uvwxi[no_u+1][no_v][no_w + 2][no_xi] + m.uvwxi[no_u+1][no_v][no_w][no_xi + 1]);

	psi[2] = a[0] * m.uvwxi[no_u][no_v+1][no_w][no_xi]
		+ a[1] * m.uvwxi[no_u + 1][no_v+1][no_w][no_xi]
		+ a[2] * m.uvwxi[no_u][no_v + 2][no_w][no_xi]
		+ a[3] * m.uvwxi[no_u][no_v+1][no_w + 1][no_xi]
		+ a[4] * 0.5*(m.uvwxi[no_u + 2][no_v+1][no_w][no_xi] + m.uvwxi[no_u][no_v + 3][no_w][no_xi] + m.uvwxi[no_u][no_v+1][no_w + 2][no_xi] + m.uvwxi[no_u][no_v+1][no_w][no_xi + 1]);

	psi[3] = a[0] * m.uvwxi[no_u][no_v][no_w+1][no_xi]
		+ a[1] * m.uvwxi[no_u + 1][no_v][no_w+1][no_xi]
		+ a[2] * m.uvwxi[no_u][no_v + 1][no_w+1][no_xi]
		+ a[3] * m.uvwxi[no_u][no_v][no_w + 2][no_xi]
		+ a[4] * 0.5*(m.uvwxi[no_u + 2][no_v][no_w+1][no_xi] + m.uvwxi[no_u][no_v + 2][no_w+1][no_xi] + m.uvwxi[no_u][no_v][no_w + 3][no_xi] + m.uvwxi[no_u][no_v][no_w+1][no_xi + 1]);

	psi[4] = 0.5*(
		a[0] * (m.uvwxi[no_u + 2][no_v][no_w][no_xi] + m.uvwxi[no_u][no_v + 2][no_w][no_xi] + m.uvwxi[no_u][no_v][no_w+2][no_xi] + m.uvwxi[no_u][no_v][no_w][no_xi + 1]) +
		a[1] * (m.uvwxi[no_u + 3][no_v][no_w][no_xi] + m.uvwxi[no_u+1][no_v + 2][no_w][no_xi] + m.uvwxi[no_u+1][no_v][no_w + 2][no_xi] + m.uvwxi[no_u+1][no_v][no_w][no_xi + 1]) +
		a[2] * (m.uvwxi[no_u + 2][no_v+1][no_w][no_xi] + m.uvwxi[no_u][no_v + 3][no_w][no_xi] + m.uvwxi[no_u][no_v+1][no_w + 2][no_xi] + m.uvwxi[no_u][no_v+1][no_w][no_xi + 1]) +
		a[3] * (m.uvwxi[no_u + 2][no_v][no_w+1][no_xi] + m.uvwxi[no_u][no_v + 2][no_w+1][no_xi] + m.uvwxi[no_u][no_v][no_w + 3][no_xi] + m.uvwxi[no_u][no_v][no_w+1][no_xi + 1]) +
		a[4] * 0.5*(m.uvwxi[no_u + 4][no_v][no_w][no_xi] + 
		            m.uvwxi[no_u][no_v + 4][no_w][no_xi] +
					m.uvwxi[no_u][no_v][no_w + 4][no_xi] +
					m.uvwxi[no_u][no_v][no_w][no_xi + 2] +
				2 * m.uvwxi[no_u + 2][no_v + 2][no_w][no_xi] + 
				2 * m.uvwxi[no_u + 2][no_v][no_w + 2][no_xi] + 
				2 * m.uvwxi[no_u + 2][no_v][no_w][no_xi + 1] +
				2 * m.uvwxi[no_u][no_v + 2][no_w + 2][no_xi] +
				2 * m.uvwxi[no_u][no_v + 2][no_w][no_xi + 1] +
				2 * m.uvwxi[no_u][no_v][no_w + 2][no_xi + 1] )
				);
}

template <class T>
void G_speed(int no_u, int no_v, int no_w, int no_xi, double *psi, double a[5], T &m)
{
	psi[0] = a[0] * m.uwhole[no_u] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[1] * m.uwhole[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[2] * m.uwhole[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[3] * m.uwhole[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[4] * 0.5*(m.uwhole[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 1]);

	psi[1] = a[0] * m.uwhole[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[1] * m.uwhole[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[2] * m.uwhole[no_u + 1] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[3] * m.uwhole[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[4] * 0.5*(m.uwhole[no_u + 3] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uwhole[no_u + 1] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uwhole[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uwhole[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 1]);

	psi[2] = a[0] * m.uwhole[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[1] * m.uwhole[no_u + 1] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[2] * m.uwhole[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[3] * m.uwhole[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[4] * 0.5*(m.uwhole[no_u + 2] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v + 3] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi + 1]);

	psi[3] = a[0] * m.uwhole[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[1] * m.uwhole[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[2] * m.uwhole[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[3] * m.uwhole[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi]
		+ a[4] * 0.5*(m.uwhole[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w + 1] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 3] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi + 1]);

	psi[4] = 0.5*(
		a[0] * (m.uwhole[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 1]) +
		a[1] * (m.uwhole[no_u + 3] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uwhole[no_u + 1] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uwhole[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uwhole[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 1]) +
		a[2] * (m.uwhole[no_u + 2] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v + 3] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi + 1]) +
		a[3] * (m.uwhole[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w + 1] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 3] * m.xi[no_xi]
			+ m.uwhole[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi + 1]) +
		a[4] * 0.5*(m.uwhole[no_u + 4] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi] +
			m.uwhole[no_u] * m.vwhole[no_v + 4] * m.wwhole[no_w] * m.xi[no_xi] +
			m.uwhole[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 4] * m.xi[no_xi] +
			m.uwhole[no_u] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 2] +
			2 * m.uwhole[no_u + 2] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi] +
			2 * m.uwhole[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi] +
			2 * m.uwhole[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 1] +
			2 * m.uwhole[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w + 2] * m.xi[no_xi] +
			2 * m.uwhole[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi + 1] +
			2 * m.uwhole[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi + 1])
		);
}

template <class T>
void GL_speed(int no_u, int no_v, int no_w, int no_xi, double *psi, double a[5], T& m)
{

	psi[0] = a[0] * m.uplus[no_u] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[1] * m.uplus[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[2] * m.uplus[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[3] * m.uplus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[4] * 0.5*(m.uplus[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 1]);

	psi[1] = a[0] * m.uplus[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[1] * m.uplus[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[2] * m.uplus[no_u + 1] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[3] * m.uplus[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[4] * 0.5*(m.uplus[no_u + 3] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uplus[no_u + 1] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uplus[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uplus[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 1]);

	psi[2] = a[0] * m.uplus[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[1] * m.uplus[no_u + 1] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[2] * m.uplus[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[3] * m.uplus[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[4] * 0.5*(m.uplus[no_u + 2] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v + 3] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi + 1]);

	psi[3] = a[0] * m.uplus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[1] * m.uplus[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[2] * m.uplus[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[3] * m.uplus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi]
		+ a[4] * 0.5*(m.uplus[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w + 1] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 3] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi + 1]);

	psi[4] = 0.5*(
		a[0] * (m.uplus[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 1]) +
		a[1] * (m.uplus[no_u + 3] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uplus[no_u + 1] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uplus[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uplus[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 1]) +
		a[2] * (m.uplus[no_u + 2] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v + 3] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi + 1]) +
		a[3] * (m.uplus[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w + 1] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 3] * m.xi[no_xi]
			+ m.uplus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi + 1]) +
		a[4] * 0.5*(m.uplus[no_u + 4] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi] +
			m.uplus[no_u] * m.vwhole[no_v + 4] * m.wwhole[no_w] * m.xi[no_xi] +
			m.uplus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 4] * m.xi[no_xi] +
			m.uplus[no_u] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 2] +
			2 * m.uplus[no_u + 2] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi] +
			2 * m.uplus[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi] +
			2 * m.uplus[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 1] +
			2 * m.uplus[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w + 2] * m.xi[no_xi] +
			2 * m.uplus[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi + 1] +
			2 * m.uplus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi + 1])
		);
}

template <class T>
void GR_speed(int no_u, int no_v, int no_w, int no_xi, double *psi, double a[5], T& m)
{

	psi[0] = a[0] * m.uminus[no_u] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[1] * m.uminus[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[2] * m.uminus[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[3] * m.uminus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[4] * 0.5*(m.uminus[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 1]);

	psi[1] = a[0] * m.uminus[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[1] * m.uminus[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[2] * m.uminus[no_u + 1] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[3] * m.uminus[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[4] * 0.5*(m.uminus[no_u + 3] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uminus[no_u + 1] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uminus[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uminus[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 1]);

	psi[2] = a[0] * m.uminus[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[1] * m.uminus[no_u + 1] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[2] * m.uminus[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi]
		+ a[3] * m.uminus[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[4] * 0.5*(m.uminus[no_u + 2] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v + 3] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi + 1]);

	psi[3] = a[0] * m.uminus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[1] * m.uminus[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[2] * m.uminus[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w + 1] * m.xi[no_xi]
		+ a[3] * m.uminus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi]
		+ a[4] * 0.5*(m.uminus[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w + 1] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 3] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi + 1]);

	psi[4] = 0.5*(
		a[0] * (m.uminus[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 1]) +
		a[1] * (m.uminus[no_u + 3] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uminus[no_u + 1] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uminus[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uminus[no_u + 1] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 1]) +
		a[2] * (m.uminus[no_u + 2] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v + 3] * m.wwhole[no_w] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w + 2] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v + 1] * m.wwhole[no_w] * m.xi[no_xi + 1]) +
		a[3] * (m.uminus[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w + 1] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 3] * m.xi[no_xi]
			+ m.uminus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 1] * m.xi[no_xi + 1]) +
		a[4] * 0.5*(m.uminus[no_u + 4] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi] +
			m.uminus[no_u] * m.vwhole[no_v + 4] * m.wwhole[no_w] * m.xi[no_xi] +
			m.uminus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 4] * m.xi[no_xi] +
			m.uminus[no_u] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 2] +
			2 * m.uminus[no_u + 2] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi] +
			2 * m.uminus[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi] +
			2 * m.uminus[no_u + 2] * m.vwhole[no_v] * m.wwhole[no_w] * m.xi[no_xi + 1] +
			2 * m.uminus[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w + 2] * m.xi[no_xi] +
			2 * m.uminus[no_u] * m.vwhole[no_v + 2] * m.wwhole[no_w] * m.xi[no_xi + 1] +
			2 * m.uminus[no_u] * m.vwhole[no_v] * m.wwhole[no_w + 2] * m.xi[no_xi + 1])
		);
}

double Get_Tau_NS(double density0, double lambda0);

double Get_Tau(double density_left, double density_right, double density0, double lambda_left, double lambda_right, double lambda0, double dt);

double Alpha(double lambda, double u);

double Beta(double lambda, double u);

void Global_to_Local(double* change, double* origin, double* normal);

void Array_zero(double* target, int dim);

bool negative_density_or_pressure(double* primvar);

void Convar_to_Primvar(Fluid2d* fluids, Block2d block);

void Convar_to_primvar(Fluid3d *fluids, Block3d& block);

void Convar_to_primvar_1D(double* primvar, double convar[3]);

void Convar_to_primvar_2D(double *primvar, double convar[4]);

void Convar_to_primvar_3D(double *primvar, double convar[5]);

void Convar_to_ULambda_1d(double* primvar, double convar[3]);

void Convar_to_ULambda_2d(double* primvar, double convar[4]);

void Primvar_to_convar_1D(double* convar, double primvar[3]);

void Primvar_to_convar_2D(double *convar, double primvar[4]);

void Primvar_to_convar_3D(double *convar, double* primvar);

void Convar_to_char1D(double* character, double primvar[3], double convar[3]);

void Char_to_convar1D(double* convar, double primvar[3], double charvar[3]);

void Convar_to_char(double *character, double *primvar, double convar[4]);

void Char_to_convar(double *convar, double *primvar, double character[4]);

void Char_base_3D(double *s, double *primvar);

void Convar_to_char_3D(double *character, double *s, double *convar);

void Char_to_convar_3D(double *convar, double *primvar, double character[5]);

double DensityU(double density, double u);

double DensityE(double density, double u, double pressure);

double Pressure(double density, double densityu, double densityE);

double Q_densityu(double density, double u);

double Q_densityv(double density, double v);

double Q_densityw(double density, double w);

double Q_densityE(double density, double u, double v, double pressure);

double Q_densityE(double density, double u, double v,double w, double pressure);

double U(double density, double q_densityu);

double V(double density, double q_densityv);

double W(double density, double q_densityw);

double Pressure(double density, double q_densityu, double q_densityv, double q_densityE);

double Pressure(double density, double q_densityu, double q_densityv, double q_densityw, double q_densityE);

double Lambda(double density, double u, double densityE);

double Lambda(double density, double u, double v, double densityE);

double Lambda(double density, double u, double v, double w, double densityE);

void Copy_Array(double* target, double* origin, int dim);

void YchangetoX(double* fluidtmp, double* fluid);

void XchangetoY(double* fluidtmp, double* fluid);

void Ydirection(double *change, double *origin);

void Zdirection(double *change, double *origin);

void Local_to_Global(double *change,double *normal);

double Max_three_number(double a1, double a2, double a3);

double Max_five_number(double a1, double a2, double a3, double a4, double a5);

double Max_seven_number(double a1, double a2, double a3, double a4, double a5, double a6, double a7);

Flux1d** Setflux_array(Block1d block);

void SetUniformMesh(Block1d block, Fluid1d* fluids, Interface1d* interfaces, Flux1d** fluxes);

Fluid2d* Setfluid(Block2d& block);

Interface2d* Setinterface_array(Block2d block);

Flux2d_gauss** Setflux_gauss_array(Block2d block);

void SetUniformMesh(Block2d& block, Fluid2d* fluids, Interface2d* xinterfaces, Interface2d* yinterfaces, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes);

void Copy_geo_from_interface_to_line(Interface2d& interface);

void Copy_geo_from_interface_to_flux(Interface2d& interface, Flux2d_gauss* flux, int stages);

void Set_Gauss_Coordinate(Interface2d& xinterface, Interface2d& yinterface);

void Set_Gauss_Intergation_Location_x(Point2d& xgauss, int index, double h);

void Set_Gauss_Intergation_Location_y(Point2d& ygauss, int index, double h);

void CopyFluid_new_to_old(Fluid1d* fluids, Block1d block);

void CopyFluid_new_to_old(Fluid2d* fluids, Block2d block);

void CopyFluid_new_to_old(Fluid3d *fluids, Block3d block);

Fluid3d *Setfluid_array(Block3d &block);

Interface3d *Setinterface_array(Block3d& block);

Flux3d_gauss** Setflux_gauss_array(Block3d& block);

void SetUniformMesh
(Block3d block, Fluid3d* fluids,
	Interface3d *xinterfaces, Interface3d *yinterfaces, Interface3d *zinterfaces);

void Set_node_for_a_cube(Fluid3d& fluid, Block3d block);

void Set_Gauss_loc_for_a_rectangular(double *loc, 
	double *node0, double *node1, double *node2, double *node3);

void Allocate_reconstruction_variable_along_interface(Interface3d & interface);
