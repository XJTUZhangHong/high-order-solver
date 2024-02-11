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