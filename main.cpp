#include "1Dproblem.h"
#include "2Dproblem.h"
#include "3Dproblem.h"

int main()
{
    omp_set_num_threads(32);
    is_using_df_factor = true;
    df_thres = 3.0;
    // 5th-order -- gausspoint = 2
    // 7th-order -- gausspoint = 4
    // 1-D problem
    //SodTubeProblem();
    //Blastwave();
    accuracy_sinwave_1d();
    //ShuOsher();
    // 2-D problem
    //RT_instability();
    //PlanarShock();
    //PlanarSheer();
    //High_mach_astrophusical_jet();
    //doubleMach();
    //viscous_sod_shock_problem();
    //accuracy_sinwave_2d();
    // 3-D problem
    //CubicTube();
    return 0;
}