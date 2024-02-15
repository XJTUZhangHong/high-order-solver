#include "1Dproblem.h"
#include "2Dproblem.h"

int main()
{
    omp_set_num_threads(6);
    is_using_df_factor = false;
    // 1-D problem
    //SodTubeProblem();
    //Blastwave();
    //accuracy_sinwave_1d();
    //ShuOsher();
    // 2-D problem
    //PlanarShock();
    High_mach_astrophusical_jet();
    return 0;
}