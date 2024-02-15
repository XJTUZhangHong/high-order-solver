#pragma once
#include"basic_function.h"
#include"time.h"
#include<string.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#if defined(_WIN32)
#include<direct.h>
#else 
#include<sys/stat.h>
//#include<sys/types.h>
#endif
struct Runtime {
	//for record wall clock time
	clock_t start_initial;
	clock_t finish_initial;
	clock_t start_compute;
	clock_t finish_compute;
	Runtime()
	{
		memset(this, 0, sizeof(Runtime));
	}
};

void output1d(Fluid1d* fluids, Block1d block);

void output_error_form(double CFL, double dt_ratio, int mesh_set, int* mesh_number, double** error);

void output2d(Fluid2d* fluids, Block2d block);
