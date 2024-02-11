#include "boundary_condition.h"

// one-dimensional problem
void free_boundary_left(Fluid1d* fluids, Block1d block, Fluid1d bcvalue)
{
	//for this case, the input variable can also use Fluid1d fluids
	//the call is the same with the current method
	for (int i = block.ghost - 1; i >= 0; i--)
	{
		fluids[i].convar[0] = fluids[i + 1].convar[0];
		fluids[i].convar[1] = fluids[i + 1].convar[1];
		fluids[i].convar[2] = fluids[i + 1].convar[2];
	}
}
void free_boundary_right(Fluid1d* fluids, Block1d block, Fluid1d bcvalue)
{
	for (int i = block.nx - block.ghost; i < block.nx; i++)
	{
		fluids[i].convar[0] = fluids[i - 1].convar[0];
		fluids[i].convar[1] = fluids[i - 1].convar[1];
		fluids[i].convar[2] = fluids[i - 1].convar[2];
	}
}