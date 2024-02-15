#pragma once
#include "flux_function.h"
#include "output.h"
#include "reconstruction.h"
#include "time_advance.h"
#include "boundary_condition.h"

void SodTubeProblem();

void ICforSodTube(Fluid1d* fluids, Block1d block);

void Blastwave();

void ICforBlastwave(Fluid1d* fluids, Block1d block);

void ShuOsher();

void ICforShuOsher(Fluid1d* fluids, Block1d block);

void accuracy_sinwave_1d();

void accuracy_sinwave_1d(double& CFL, double& dt_ratio, int& mesh_number, double* error);

void ICfor_sinwave(Fluid1d* fluids, Block1d block);

void error_for_sinwave(Fluid1d* fluids, Block1d block, double tstop, double* error);

