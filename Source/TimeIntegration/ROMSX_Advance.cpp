#include <ROMSX.H>

using namespace amrex;

// advance a single level for a single time step
 void
ROMSX::Advance (int lev, Real /*time*/, Real dt_lev, int /*iteration*/, int /*ncycle*/)
{

    BL_PROFILE("ROMSX::Advance()");

    setup_step(lev, dt_lev);

    if (solverChoice.use_barotropic)
    {
        bool predictor_2d_step=true;
        int nfast_counter=nfast + 1;
        //Compute fast timestep from dt_lev and ratio
        Real dtfast_lev=dt_lev/Real(fixed_ndtfast_ratio);
        for(int my_iif = 0; my_iif < nfast_counter; my_iif++) {
            advance_2d_onestep(lev, dt_lev, dtfast_lev, my_iif, nfast_counter);
        }
    }

    advance_3d_ml(lev, dt_lev);
}
