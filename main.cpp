#include "1Dproblem.h"
#include "2Dproblem.h"

int main()
{
    omp_set_num_threads(6);
    SodTubeProblem();
    return 0;
}