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

void reflection_boundary_left(Fluid1d* fluids, Block1d block, Fluid1d bcvalue)
{
	for (int i = block.ghost - 1; i >= 0; i--)
	{
		int ref = 2 * block.ghost - 1 - i;
		fluids[i].convar[0] = fluids[ref].convar[0];
		fluids[i].convar[1] = -fluids[ref].convar[1];
		fluids[i].convar[2] = fluids[ref].convar[2];
	}
}

void reflection_boundary_right(Fluid1d* fluids, Block1d block, Fluid1d bcvalue)
{
	for (int i = block.nx - block.ghost; i < block.nx; i++)
	{
		int ref = 2 * (block.nx - block.ghost) - 1 - i;
		fluids[i].convar[0] = fluids[ref].convar[0];
		fluids[i].convar[1] = -fluids[ref].convar[1];
		fluids[i].convar[2] = fluids[ref].convar[2];
	}
}

// two-dimensional problem
void free_boundary_left(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int i = block.ghost - 1; i >= 0; i--)
	{
		for (int j = 0; j < block.ny; j++)
		{
			fluids[i * block.ny + j].convar[0] = fluids[(i + 1) * block.ny + j].convar[0];
			fluids[i * block.ny + j].convar[1] = fluids[(i + 1) * block.ny + j].convar[1];
			fluids[i * block.ny + j].convar[2] = fluids[(i + 1) * block.ny + j].convar[2];
			fluids[i * block.ny + j].convar[3] = fluids[(i + 1) * block.ny + j].convar[3];
		}
	}
}

void free_boundary_right(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int i = block.nx - block.ghost; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			fluids[i * block.ny + j].convar[0] = fluids[(i - 1) * block.ny + j].convar[0];
			fluids[i * block.ny + j].convar[1] = fluids[(i - 1) * block.ny + j].convar[1];
			fluids[i * block.ny + j].convar[2] = fluids[(i - 1) * block.ny + j].convar[2];
			fluids[i * block.ny + j].convar[3] = fluids[(i - 1) * block.ny + j].convar[3];
		}

	}

}

void free_boundary_down(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int j = block.ghost - 1; j >= 0; j--)
	{
		for (int i = 0; i < block.nx; i++)
		{
			fluids[i * block.ny + j].convar[0] = fluids[i * block.ny + j + 1].convar[0];
			fluids[i * block.ny + j].convar[1] = fluids[i * block.ny + j + 1].convar[1];
			fluids[i * block.ny + j].convar[2] = fluids[i * block.ny + j + 1].convar[2];
			fluids[i * block.ny + j].convar[3] = fluids[i * block.ny + j + 1].convar[3];
		}
	}

}

void free_boundary_up(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int j = block.ny - block.ghost; j < block.ny; j++)
	{
		for (int i = 0; i < block.nx; i++)
		{
			fluids[i * block.ny + j].convar[0] = fluids[i * block.ny + j - 1].convar[0];
			fluids[i * block.ny + j].convar[1] = fluids[i * block.ny + j - 1].convar[1];
			fluids[i * block.ny + j].convar[2] = fluids[i * block.ny + j - 1].convar[2];
			fluids[i * block.ny + j].convar[3] = fluids[i * block.ny + j - 1].convar[3];
		}
	}
}