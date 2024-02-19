#include "fluid_mesh.h"

int gausspoint = 0; //initilization
double* gauss_loc = new double[gausspoint]; //initilization
double* gauss_weight = new double[gausspoint]; //initilization
double** gauss_loc_3d;
bool is_multi_dimensional = true;

void SetGuassPoint()
{
	if (gausspoint == 0)
	{
		cout << "no gausspoint specify" << endl;
		exit(0);
	}
	if (gausspoint == 1)
	{
		gauss_loc = new double[1]; gauss_loc[0] = 0.0;
		gauss_weight = new double[1]; gauss_weight[0] = 1.0;
	}
	if (gausspoint == 2)
	{
		gauss_loc = new double[2]; gauss_loc[0] = -sqrt(1.0 / 3.0); gauss_loc[1] = -gauss_loc[0];
		gauss_weight = new double[2]; gauss_weight[0] = 0.5; gauss_weight[1] = 0.5;
	}
	if (gausspoint == 3)
	{
		cout << "the usage of 3 gausspoint may be imcompatiable for some reconstructions" << endl;
		//exit(0);
		gauss_loc = new double[3]; gauss_loc[0] = 0.0; gauss_loc[1] = sqrt(3.0 / 5.0); gauss_loc[2] = -gauss_loc[1];
		gauss_weight = new double[3];
		gauss_weight[0] = 4.0 / 9.0; gauss_weight[1] = 0.5 - 0.5 * gauss_weight[0]; gauss_weight[2] = gauss_weight[1];
	}
	if (gausspoint == 4)
	{
		cout << "the usage of 4 gausspoint may be imcompatiable for some reconstructions" << endl;
		gauss_loc = new double[4]; gauss_weight = new double[4];

		gauss_loc[0] = -1.0; gauss_loc[1] = -gauss_loc[0];
		gauss_loc[2] = -sqrt(5) / 5.0; gauss_loc[3] = -gauss_loc[2];
		gauss_weight[0] = 1.0 / 12.0;		gauss_weight[1] = 1.0 / 12.0;
		gauss_weight[2] = 5.0 / 12.0;		gauss_weight[3] = 5.0 / 12.0;

	}
	cout << gausspoint << " gausspoint(s) specify" << endl;
}

void SetGuassPoint_3D()
{
	if (is_multi_dimensional == false)
	{
		gausspoint = 0;
		cout << " if you set multi_demensional_false, gauuspoint shall be zero then" << endl;
	}
	if (gausspoint == 0)
	{
		cout << "no gausspoint specify" << endl;
		gauss_weight = new double[1]; gauss_weight[0] = 1.0;
		return;
	}
	if (gausspoint == 2 || gausspoint == 3 || gausspoint == 5 ||
		gausspoint == 6 || gausspoint == 7 || gausspoint == 8 ||
		gausspoint == 10 || gausspoint == 11 || gausspoint == 12 ||
		gausspoint == 13 || gausspoint == 14 || gausspoint == 15)
	{
		cout << "gausspoint is equal to " << gausspoint << " not resonable gausspoint set in 3-D" << endl;
		exit(0);
	}

	gauss_loc_3d = new double*[gausspoint];
	for (int i = 0; i < gausspoint; i++)
	{
		gauss_loc_3d[i] = new double[2];
	}
	gauss_weight = new double[gausspoint];
	if (gausspoint == 1)
	{
		 gauss_loc_3d[0][0] = 0.0; gauss_loc_3d[0][1] = 0.0;
		 gauss_weight[0] = 1.0;
	}

	if (gausspoint == 4)
	{	
		gauss_loc_3d[0][0] = -sqrt(1.0 / 3.0); gauss_loc_3d[0][1] = gauss_loc_3d[0][0];
		gauss_loc_3d[1][0] = gauss_loc_3d[0][0]; gauss_loc_3d[1][1] = -gauss_loc_3d[0][0];
		gauss_loc_3d[2][0] = -gauss_loc_3d[0][0]; gauss_loc_3d[2][1] = gauss_loc_3d[0][0];
		gauss_loc_3d[3][0] = -gauss_loc_3d[0][0]; gauss_loc_3d[3][1] = -gauss_loc_3d[0][0];

		gauss_weight[0] = 0.25; gauss_weight[1] = 0.25; gauss_weight[2] = 0.25; gauss_weight[3] = 0.25;
	}

	if (gausspoint == 9)
	{
		cout << "the usage of 3 gausspoint may result unkown error" << endl;
		//exit(0);
		double k1 = 2.5 / 9.0; double k2 = 4.0 / 9.0;
		gauss_weight[0] = k1*k1; gauss_weight[1] = k1*k2; gauss_weight[2] = k1*k1;
		gauss_weight[3] = k2 * k1; gauss_weight[4] = k2 * k2; gauss_weight[5] = k2 * k1;
		gauss_weight[6] = k1 * k1; gauss_weight[7] = k1 * k2; gauss_weight[8] = k1 * k1;

		gauss_loc_3d[0][0] = -sqrt(3.0 / 5.0); gauss_loc_3d[0][1] = -sqrt(3.0 / 5.0);
		gauss_loc_3d[1][0] = -sqrt(3.0 / 5.0); gauss_loc_3d[1][1] = 0.0;
		gauss_loc_3d[2][0] = -sqrt(3.0 / 5.0); gauss_loc_3d[2][1] = sqrt(3.0 / 5.0);
		gauss_loc_3d[3][0] = 0.0; gauss_loc_3d[3][1] = -sqrt(3.0 / 5.0);
		gauss_loc_3d[4][0] = 0.0; gauss_loc_3d[4][1] = 0.0;
		gauss_loc_3d[5][0] = 0.0; gauss_loc_3d[5][1] = sqrt(3.0 / 5.0);
		gauss_loc_3d[6][0] = sqrt(3.0 / 5.0); gauss_loc_3d[6][1] = -sqrt(3.0 / 5.0);
		gauss_loc_3d[7][0] = sqrt(3.0 / 5.0); gauss_loc_3d[7][1] = 0.0;
		gauss_loc_3d[8][0] = sqrt(3.0 / 5.0); gauss_loc_3d[8][1] = sqrt(3.0 / 5.0);
		//gauss_loc[0] = 0.0; gauss_loc[1] = -1; gauss_loc[2] = 1;
		//gauss_weight[0] = 2.0/3.0; gauss_weight[1] = 1.0/6.0; gauss_weight[2] = gauss_weight[1];
	}
	cout << gausspoint << " gausspoint(s) specify in 3-D case" << endl;
}