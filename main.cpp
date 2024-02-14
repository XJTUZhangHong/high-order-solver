#include "1Dproblem.h"
#include "2Dproblem.h"

int main()
{
    omp_set_num_threads(6);
    // 1-D problem
    //SodTubeProblem();
    //Blastwave();
    //ShuOsher();
    // 2-D problem
    PlanarShock();
    return 0;
}