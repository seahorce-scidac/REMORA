#include <REMORA.H>

using namespace amrex;

//
// Advance all levels by dt
//
void
REMORA::timeStepML (Real time, int /*iteration*/)
{
    if (time == 0.0_rt && solverChoice.init_l1ad_T) {
        average_down(*cons_new[1], *cons_new[0],
                 0, cons_new[0]->nComp(), refRatio(0));
        Print() << "extra average_down " << istep[0] << std::endl;
        WritePlotFile();
    }

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
                        dt[k] = dt[k-1] / nsubsteps[lev];
                    }
                }
            }
        }
    }

    scale_rhs_vars();

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

        // **************************************************************************************
        // Register old and new coarse data if we are at a level less than the finest level
        // **************************************************************************************
        if (lev < finest_level)
        {
            if (cf_width > 0) {
                // We must fill the ghost cells of these so that the parallel copy works correctly
                cons_old[lev]->FillBoundary(geom[lev].periodicity());
                cons_new[lev]->FillBoundary(geom[lev].periodicity());
                FPr_c[lev].RegisterCoarseData({cons_old[lev], cons_new[lev]}, {time, time + dt[lev]});
            }

            if (cf_width >= 0) {
                // We must fill the ghost cells of these so that the parallel copy works correctly
                xvel_old[lev]->FillBoundary(geom[lev].periodicity());
                xvel_new[lev]->FillBoundary(geom[lev].periodicity());
                FPr_u[lev].RegisterCoarseData({xvel_old[lev], xvel_new[lev]}, {time, time + dt[lev]});

                yvel_old[lev]->FillBoundary(geom[lev].periodicity());
                yvel_new[lev]->FillBoundary(geom[lev].periodicity());
                FPr_v[lev].RegisterCoarseData({yvel_old[lev], yvel_new[lev]}, {time, time + dt[lev]});

                zvel_old[lev]->FillBoundary(geom[lev].periodicity());
                zvel_new[lev]->FillBoundary(geom[lev].periodicity());
                FPr_w[lev].RegisterCoarseData({zvel_old[lev], zvel_new[lev]}, {time, time + dt[lev]});
            }
        }
    }

    if (solverChoice.use_barotropic)
    {
        int nfast_counter=nfast + 1;

        for (int my_iif = 0; my_iif < nfast_counter; my_iif++) {
            //Compute fast timestep from dt[lev] and ratio
            for (int lev=0; lev <= finest_level; lev++)
            {
                Real dtfast_lev=dt[lev]/Real(fixed_ndtfast_ratio);
                advance_2d_onestep(lev, dt[lev], dtfast_lev, my_iif, nfast_counter);
                // **************************************************************************************
                // Register old and new coarse data if we are at a level less than the finest level
                // **************************************************************************************
                if (lev < finest_level)
                {
                    if (cf_width >= 0) {
                        Print() << "cf width >= 0  " << dt[lev] << std::endl;
                        // We must fill the ghost cells of these so that the parallel copy works correctly
                        vec_ubar[lev]->FillBoundary(geom[lev].periodicity());
                        FPr_ubar[lev].RegisterCoarseData({vec_ubar[lev].get(), vec_ubar[lev].get()}, {time, time + dt[lev]});

                        vec_vbar[lev]->FillBoundary(geom[lev].periodicity());
                        FPr_vbar[lev].RegisterCoarseData({vec_vbar[lev].get(), vec_vbar[lev].get()}, {time, time + dt[lev]});
                    }
                }
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
        // **************************************************************************************
        // Register old and new coarse data if we are at a level less than the finest level
        // **************************************************************************************
        if (lev < finest_level)
        {
            if (cf_width > 0) {
                // We must fill the ghost cells of these so that the parallel copy works correctly
                cons_old[lev]->FillBoundary(geom[lev].periodicity());
                cons_new[lev]->FillBoundary(geom[lev].periodicity());
                FPr_c[lev].RegisterCoarseData({cons_old[lev], cons_new[lev]}, {time, time + dt[lev]});
            }

            if (cf_width >= 0) {
                Print() << "cf width >= 0  " << dt[lev] << std::endl;
                // We must fill the ghost cells of these so that the parallel copy works correctly
                xvel_old[lev]->FillBoundary(geom[lev].periodicity());
                xvel_new[lev]->FillBoundary(geom[lev].periodicity());
                FPr_u[lev].RegisterCoarseData({xvel_old[lev], xvel_new[lev]}, {time, time + dt[lev]});

                yvel_old[lev]->FillBoundary(geom[lev].periodicity());
                yvel_new[lev]->FillBoundary(geom[lev].periodicity());
                FPr_v[lev].RegisterCoarseData({yvel_old[lev], yvel_new[lev]}, {time, time + dt[lev]});

                zvel_old[lev]->FillBoundary(geom[lev].periodicity());
                zvel_new[lev]->FillBoundary(geom[lev].periodicity());
                FPr_w[lev].RegisterCoarseData({zvel_old[lev], zvel_new[lev]}, {time, time + dt[lev]});
            }
        }
        FillPatch(lev, t_new[lev], *xvel_new[lev], xvel_new, BdyVars::u,0,true,true);
        FillPatch(lev, t_new[lev], *yvel_new[lev], yvel_new, BdyVars::v,0,true,true);
        FillPatch(lev, t_new[lev], *zvel_new[lev], zvel_new, BdyVars::null,0,true,true);
    }

    print_state(*xvel_new[0],IntVect(20,13,0));
    print_state(*xvel_new[1],IntVect(60,42,0));
    print_state(*xvel_new[1],IntVect(59,42,0));
    print_state(*xvel_new[1],IntVect(58,42,0));
    print_state(*xvel_new[1],IntVect(57,42,0));
    scale_rhs_vars_inv();

    if (solverChoice.coupling_type == CouplingType::TwoWay) {
        for (int lev=0; lev <= finest_level-1; lev++) {
            AverageDownTo(lev); // average lev+1 down to lev
        }
    }
}
