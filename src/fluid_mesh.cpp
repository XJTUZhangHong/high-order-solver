#include "fluid_mesh.h"

int gausspoint = 0; //initilization
double* gauss_loc = new double[gausspoint]; //initilization
double* gauss_weight = new double[gausspoint]; //initilization

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