#include <REMORA.H>

using namespace amrex;

//
// Advance all levels by dt
//
void
REMORA::timeStepML (Real time, int /*iteration*/)
{
    // HACK HACK so lev is defined and compiler won't complain, but always say regrid_int=-1
    for (int lev=0; lev <= finest_level;lev++) {
        if (regrid_int > 0)  // We may need to regrid
        {
            // help keep track of whether a level was already regridded
            // from a coarser level call to regrid
            static Vector<int> last_regrid_step(max_level+1, 0);

            // regrid changes level "lev+1" so we don't regrid on max_level
            // also make sure we don't regrid fine levels again if
            // it was taken care of during a coarser regrid
            if (lev < max_level && istep[lev] > last_regrid_step[lev])
            {
                if (istep[lev] % regrid_int == 0)
                {
                    // regrid could add newly refine levels (if finest_level < max_level)
                    // so we save the previous finest level index
                    int old_finest = finest_level;
                    regrid(lev, time);

                    // Mark that we have regridded this level already
                    for (int k = lev; k <= finest_level; ++k) {
                        last_regrid_step[k] = istep[k];
                    }

                    // If there are newly created levels, set the time step
                    for (int k = old_finest+1; k <= finest_level; ++k) {
                        dt[k] = dt[k-1] / MaxRefRatio(k-1);
                    }
                }
            }
        }
    }

    for (int lev=0; lev <= finest_level;lev++)
    {
        // Update what we call "old" and "new" time
        t_old[lev] = t_new[lev];
        t_new[lev] += dt[lev];

        if (Verbose()) {
            amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
            amrex::Print() << "ADVANCE from time = " << t_old[lev] << " to " << t_new[lev]
                           << " with dt = " << dt[lev] << std::endl;
        }

        // We must swap the pointers so the previous step's "new" is now this step's "old"
        std::swap(cons_old[lev], cons_new[lev]);
        std::swap(xvel_old[lev], xvel_new[lev]);
        std::swap(yvel_old[lev], yvel_new[lev]);
        std::swap(zvel_old[lev], zvel_new[lev]);

        setup_step(lev, time, dt[lev]);
    }

    if (solverChoice.use_barotropic)
    {
        int nfast_counter=nfast + 1;

        for (int lev=0; lev <= finest_level; lev++)
        {
            //Compute fast timestep from dt_lev and ratio
            Real dtfast_lev=dt[lev]/Real(fixed_ndtfast_ratio);
            for (int my_iif = 0; my_iif < nfast_counter; my_iif++) {
                advance_2d_onestep(lev, dt[lev], dtfast_lev, my_iif, nfast_counter);
            } // my_iif
        } // lev
    } // use_barotropic

    for (int lev=0; lev <= finest_level; lev++) {
        advance_3d_ml(lev, dt[lev]);
        ++istep[lev];

        if (Verbose())
        {
            amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
            amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
        }
    }

    if (solverChoice.coupling_type == CouplingType::TwoWay) {
        for (int lev=0; lev <= finest_level-1; lev++) {
            AverageDownTo(lev); // average lev+1 down to lev
        }
    }
}
