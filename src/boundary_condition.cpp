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

void periodic_boundary_left(Fluid1d* fluids, Block1d block, Fluid1d bcvalue)
{
	for (int i = block.nodex_begin - 1; i >= 0; i--)
	{
		int ref = i + block.nodex;
		Copy_Array(fluids[i].convar, fluids[ref].convar, 3);
	}
}

void periodic_boundary_right(Fluid1d* fluids, Block1d block, Fluid1d bcvalue)
{
	for (int i = block.nodex_end + 1; i < block.nx; i++)
	{
		int ref = i - block.nodex;
		Copy_Array(fluids[i].convar, fluids[ref].convar, 3);
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

void periodic_boundary_left(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int i = block.ghost - 1; i >= 0; i--)
	{
		int tar = i + block.nodex;
		for (int j = 0; j < block.ny; j++)
		{
			for (int var = 0; var < 4; var++)
			{
				fluids[i * block.ny + j].convar[var] = fluids[tar * block.ny + j].convar[var];
			}
		}
	}

}

void periodic_boundary_right(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int i = block.nx - block.ghost; i < block.nx; i++)
	{
		int tar = i - block.nodex;
		for (int j = 0; j < block.ny; j++)
		{
			for (int var = 0; var < 4; var++)
			{
				fluids[i * block.ny + j].convar[var] = fluids[tar * block.ny + j].convar[var];
			}
		}
	}

}

void periodic_boundary_down(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int j = block.ghost - 1; j >= 0; j--)
	{
		int tar = j + block.nodey;
		for (int i = 0; i < block.nx; i++)
		{
			for (int var = 0; var < 4; var++)
			{
				fluids[i * block.ny + j].convar[var] = fluids[i * block.ny + tar].convar[var];
			}
		}
	}

}

void periodic_boundary_up(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int j = block.ny - block.ghost; j < block.ny; j++)
	{
		int tar = j - block.nodey;
		for (int i = 0; i < block.nx; i++)
		{
			for (int var = 0; var < 4; var++)
			{
				fluids[i * block.ny + j].convar[var] = fluids[i * block.ny + tar].convar[var];
			}
		}
	}
}

void RT_boundary(Fluid2d* fluids, Block2d block)
{
	// left
	int order;
	order = block.ghost;
	for (int i = order - 1; i >= 0; i--)
	{
		for (int j = 0; j < block.ny; j++)
		{
			int index = i * block.ny + j;
			int ref = (2 * order - 1 - i) * block.ny + j;
			for (int k = 0; k < 4; k++)
			{
				fluids[index].convar[k] = fluids[ref].convar[k];
			}
			fluids[index].convar[1] = -fluids[ref].convar[1];
		}
	}
	// right
	for (int i = block.nx - order; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			int index = i * block.ny + j;
			int ref = (2 * block.nx - 2 * order - 1 - i) * block.ny + j;
			for (int k = 0; k < 4; k++)
			{
				fluids[index].convar[k] = fluids[ref].convar[k];
			}
			fluids[index].convar[1] = -fluids[ref].convar[1];
		}
	}
	// down
	for (int j = block.ghost - 1; j >= 0; j--)
	{
		for (int i = 0; i < block.nx; i++)
		{
			fluids[i * block.ny + j].convar[0] = 2.0;
			fluids[i * block.ny + j].convar[1] = 0.0;
			fluids[i * block.ny + j].convar[2] = 0.0;
			fluids[i * block.ny + j].convar[3] = 1.0 / (Gamma - 1.0);
		}
	}
	// up
	for (int j = block.ny - block.ghost; j < block.ny; j++)
	{
		for (int i = 0; i < block.nx; i++)
		{
			fluids[i * block.ny + j].convar[0] = 1.0;
			fluids[i * block.ny + j].convar[1] = 0.0;
			fluids[i * block.ny + j].convar[2] = 0.0;
			fluids[i * block.ny + j].convar[3] = 2.5 / (Gamma - 1);
		}
	}
}

// three-dimensional problem
void free_boundary_xleft(Fluid3d *fluids, Block3d block, Fluid3d bcvalue)
{
#pragma omp parallel for
	for (int j = 0; j < block.ny; j++)
	{
		for (int i = block.ghost - 1; i >= 0; i--)
		{
			for (int k = 0; k < block.nz; k++)
			{
				int index = i*(block.ny*block.nz) + j*block.nz + k;
				int ref = block.ghost*(block.ny*block.nz) + j*block.nz + k;
				for (int var = 0; var < 5; var++)
				{
					fluids[index].convar[var] = fluids[ref].convar[var];
				}
			}
		}
	}
}

void free_boundary_xright(Fluid3d *fluids, Block3d block, Fluid3d bcvalue)
{
#pragma omp parallel for
	for (int j = 0; j < block.ny; j++)
	{
		for (int i = block.nx - block.ghost; i < block.nx; i++)
		{
			for (int k = 0; k < block.nz; k++)
			{
				int index = i*(block.ny*block.nz) + j*block.nz + k;
				int ref = (block.nx - block.ghost-1)*(block.ny*block.nz) + j*block.nz + k;
				for (int var = 0; var < 5; var++)
				{
					fluids[index].convar[var] = fluids[ref].convar[var];
				}

			}
		}

	}
}

void free_boundary_yleft(Fluid3d *fluids, Block3d block, Fluid3d bcvalue)
{
#pragma omp parallel for
	for (int i = 0; i < block.nx; i++)
	{
		for (int k = 0; k < block.nz; k++)
		{
			for (int j = block.ghost - 1; j >= 0; j--)
			{

				int index = i*(block.ny*block.nz) + j*block.nz + k;
				int ref = i*(block.ny*block.nz) + (block.ghost)*block.nz + k;
				for (int var = 0; var < 5; var++)
				{
					fluids[index].convar[var] = fluids[ref].convar[var];
				}
			}

		}
	}
}

void free_boundary_yright(Fluid3d *fluids, Block3d block, Fluid3d bcvalue)
{
#pragma omp parallel for
	for (int i = 0; i < block.nx; i++)
	{
		for (int k = 0; k < block.nz; k++)
		{
			for (int j = block.ny - block.ghost; j < block.ny; j++)
			{
				int index = i*(block.ny*block.nz) + j*block.nz + k;
				int ref = i*(block.ny*block.nz) + (block.ny - block.ghost-1)*block.nz + k;
				for (int var = 0; var < 5; var++)
				{
					fluids[index].convar[var] = fluids[ref].convar[var];
				}
			}

		}
	}
}

void free_boundary_zleft(Fluid3d *fluids, Block3d block, Fluid3d bcvalue)
{
#pragma omp parallel for
	for (int j = 0; j < block.ny; j++)
	{
		for (int i = 0; i < block.nx; i++)
		{
			for (int k = block.ghost - 1; k >= 0; k--)
			{
				int index = i*(block.ny*block.nz) + j*block.nz + k;
				int ref = i*(block.ny*block.nz) + j*block.nz + block.ghost;
				for (int var = 0; var < 5; var++)
				{
					fluids[index].convar[var] = fluids[ref].convar[var];
				}
			}

		}
	}
}

void free_boundary_zright(Fluid3d *fluids, Block3d block, Fluid3d bcvalue)
{
#pragma omp parallel for
	for (int j = 0; j < block.ny; j++)
	{
		for (int i = 0; i < block.nx; i++)
		{
			for (int k = block.nz - block.ghost; k < block.nz; k++)
			{
				int index = i*(block.ny*block.nz) + j*block.nz + k;
				int ref = i*(block.ny*block.nz) + j*block.nz + block.nz - block.ghost -1;
				for (int var = 0; var < 5; var++)
				{
					fluids[index].convar[var] = fluids[ref].convar[var];
				}
			}
		}

	}
}