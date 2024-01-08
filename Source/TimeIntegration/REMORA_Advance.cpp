#include <REMORA.H>

using namespace amrex;

// advance a single level for a single time step
 void
REMORA::Advance (int lev, Real time, Real dt_lev, int /*iteration*/, int /*ncycle*/)
{
    BL_PROFILE("REMORA::Advance()");

    setup_step(lev, time, dt_lev);

    if (solverChoice.use_barotropic)
    {
        int nfast_counter=nfast + 1;

        //***************************************************
        //Compute fast timestep from dt_lev and ratio
        //***************************************************
        Real dtfast_lev=dt_lev/Real(fixed_ndtfast_ratio);

        //***************************************************
        //Advance nfast_counter steps of the 2d integrator
        //***************************************************
        for (int my_iif = 0; my_iif < nfast_counter; my_iif++) {
            advance_2d_onestep(lev, dt_lev, dtfast_lev, my_iif, nfast_counter);
        }
    }

    //***************************************************
    //Advance one step of the 3d integrator
    //***************************************************
    advance_3d_ml(lev, dt_lev);
}
