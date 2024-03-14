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

    // **************************************************************************************
    // Register old and new coarse data if we are at a level less than the finest level
    // **************************************************************************************
    if (lev < finest_level)
    {
        if (cf_width > 0) {
            // We must fill the ghost cells of these so that the parallel copy works correctly
            cons_old[lev]->FillBoundary(geom[lev].periodicity());
            cons_new[lev]->FillBoundary(geom[lev].periodicity());
            FPr_c[lev].RegisterCoarseData({cons_old[lev], cons_new[lev]}, {time, time + dt_lev});
        }

        if (cf_width >= 0) {
            // We must fill the ghost cells of these so that the parallel copy works correctly
            xvel_old[lev]->FillBoundary(geom[lev].periodicity());
            xvel_new[lev]->FillBoundary(geom[lev].periodicity());
            FPr_u[lev].RegisterCoarseData({xvel_old[lev], xvel_new[lev]}, {time, time + dt_lev});

            yvel_old[lev]->FillBoundary(geom[lev].periodicity());
            yvel_new[lev]->FillBoundary(geom[lev].periodicity());
            FPr_v[lev].RegisterCoarseData({yvel_old[lev], yvel_new[lev]}, {time, time + dt_lev});

            zvel_old[lev]->FillBoundary(geom[lev].periodicity());
            zvel_new[lev]->FillBoundary(geom[lev].periodicity());
            FPr_w[lev].RegisterCoarseData({zvel_old[lev], zvel_new[lev]}, {time, time + dt_lev});
        }
    }
}
