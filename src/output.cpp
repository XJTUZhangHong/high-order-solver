#include "output.h"

void output1d(Fluid1d* fluids, Block1d block)
{
	ofstream ResultFile;
	int N = block.nodex, ghost = block.ghost;
	const char* filePath_R = "../../data/R.txt";

	ResultFile.open(filePath_R);
	for (int i = ghost; i < N + ghost; i++)
	{
		ResultFile << fluids[i].convar[0] << endl;
	}
	ResultFile.close();
}

void output2d(Fluid2d* fluids, Block2d block)
{
	int N = block.nodex, ghost = block.ghost;
	int n = block.nodex + 2 * block.ghost;

	// Lunix File Output Directory
	ofstream ofs;
	ofs.open("../../data/R2d.txt", ios::trunc);
	for (int j = block.ghost; j < block.ghost + block.nodey; j++)
	{
		for (int i = block.ghost; i < block.ghost + block.nodex; i++)
		{
			ofs << fluids[i * block.ny + j].convar[0] << endl;
		}
	}
	ofs.close();
}
