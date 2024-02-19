#include "output.h"

void output1d(Fluid1d* fluids, Block1d block)
{
	ofstream ResultFile;
	int N = block.nodex, ghost = block.ghost;
	const char* filePath_R = "../data/R.txt";

	ResultFile.open(filePath_R);
	for (int i = ghost; i < N + ghost; i++)
	{
		ResultFile << fluids[i].convar[0] << endl;
	}
	ResultFile.close();
}

void output_error_form(double CFL, double dt_ratio, int mesh_set, int* mesh_number, double** error)
{
	cout << "Accuracy-residual-file-output" << endl;

	ofstream error_out("../data/error.txt");
	error_out << "the CFL number is " << CFL << endl;
	error_out << "the dt ratio over dx is " << dt_ratio << endl;


	double** order = new double* [mesh_set - 1];
	for (int i = 0; i < mesh_set - 1; i++)
	{
		order[i] = new double[3];
		for (int j = 0; j < 3; j++)
		{
			order[i][j] = log(error[i][j] / error[i + 1][j]) / log(2);
		}
	}


	for (int i = 0; i < mesh_set; i++)
	{
		if (i == 0)
		{
			error_out << "1/" << mesh_number[i]
				<< scientific << setprecision(6)
				<< " & " << error[i][0] << " & ~"
				<< " & " << error[i][1] << " & ~"
				<< " & " << error[i][2] << " & ~"
				<< " \\\\ " << endl;
		}
		else
		{
			error_out << "1/" << mesh_number[i]
				<< " & " << scientific << setprecision(6) << error[i][0];
			error_out << " & " << fixed << setprecision(2) << order[i - 1][0];
			error_out << " & " << scientific << setprecision(6) << error[i][1];
			error_out << " & " << fixed << setprecision(2) << order[i - 1][1];
			error_out << " & " << scientific << setprecision(6) << error[i][2];
			error_out << " & " << fixed << setprecision(2) << order[i - 1][2];
			error_out << " \\\\ " << endl;
		}

	}
	error_out.close();
}

void output2d(Fluid2d* fluids, Block2d block)
{
	int N = block.nodex, ghost = block.ghost;
	int n = block.nodex + 2 * block.ghost;

	// Lunix File Output Directory
	ofstream ofs;
	ofs.open("../data/R2d.txt", ios::trunc);
	for (int j = block.ghost; j < block.ghost + block.nodey; j++)
	{
		for (int i = block.ghost; i < block.ghost + block.nodex; i++)
		{
			ofs << fluids[i * block.ny + j].convar[0] << endl;
		}
	}
	ofs.close();
}

void output3d(Fluid3d* fluids, Block3d block)
{
	ofstream ofs;
	ofs.open("../data/R3d.txt", ios::trunc);

	for (int k = block.ghost; k < block.ghost + block.nodez; k++)
	{
		for (int j = block.ghost; j < block.ghost + block.nodey; j++)
		{
			for (int i = block.ghost; i < block.ghost + block.nodex; i++)
			{
				int index = i*(block.ny*block.nz) + j*block.nz + k;

				ofs << fluids[index].convar[0] << endl;
			}
		}
	}
	ofs.close();
}






