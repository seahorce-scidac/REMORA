#include <ROMSX.H>
#include <ROMSX_PhysBCFunct.H>
#include <IndexDefines.H>
#include <TimeInterpolatedData.H>

using namespace amrex;

PhysBCFunctNoOp null_bc;

//
// Fill valid and ghost data in the MultiFab "mf"
// This version fills the MultiFab mf in valid regions with the "state data" at the given time;
// values in mf when it is passed in are *not* used.
//
void
ROMSX::FillPatch (int lev, Real time, MultiFab& mf_to_fill, Vector<MultiFab*> const& mfs,
                  const int bdy_var_type)
{
    BL_PROFILE_VAR("ROMSX::FillPatch()",ROMSX_FillPatch);
    int bccomp;
    amrex::Interpolater* mapper = nullptr;

    const int icomp = 0;
    const int ncomp = mf_to_fill.nComp();

    Box mf_box(mf_to_fill.boxArray()[0]);
    if (mf_box.ixType() == IndexType(IntVect(0,0,0)))
    {
        bccomp = 0;
        mapper = &cell_cons_interp;
    }
    else if (mf_box.ixType() == IndexType(IntVect(1,0,0)))
    {
        bccomp = BCVars::xvel_bc;
        mapper = &face_linear_interp;
    }
    else if (mf_box.ixType() == IndexType(IntVect(0,1,0)))
    {
        bccomp = BCVars::yvel_bc;
        mapper = &face_linear_interp;
    }
    else if (mf_box.ixType() == IndexType(IntVect(0,0,1)))
    {
        bccomp = BCVars::zvel_bc;
        mapper = &face_linear_interp;
    }
    else
    {
        amrex::Abort("Dont recognize this box type in ROMSX_FillPatch");
    }

    if (lev == 0)
    {
        Vector<MultiFab*> fmf = {mfs[lev], mfs[lev]};
        Vector<Real> ftime    = {t_old[lev], t_new[lev]};
        amrex::FillPatchSingleLevel(mf_to_fill, time, fmf, ftime, icomp, icomp, ncomp,
                                    geom[lev], null_bc, bccomp);
    }
    else
    {
        Vector<MultiFab*> fmf = {mfs[lev], mfs[lev]};
        Vector<Real> ftime    = {t_old[lev], t_new[lev]};
        Vector<MultiFab*> cmf = {mfs[lev-1], mfs[lev-1]};
        Vector<Real> ctime    = {t_old[lev-1], t_new[lev-1]};

        amrex::FillPatchTwoLevels(mf_to_fill, time, cmf, ctime, fmf, ftime,
                                  0, icomp, ncomp, geom[lev-1], geom[lev],
                                  null_bc, bccomp, null_bc, bccomp, refRatio(lev-1),
                                  mapper, domain_bcs_type, bccomp);
    } // lev > 0

    // ***************************************************************************
    // Physical bc's at domain boundary
    // ***************************************************************************

#ifdef ROMSX_USE_NETCDF
    // Fill the data which is stored in the boundary data read from netcdf files
    if ( (solverChoice.ic_bc_type == IC_BC_Type::Real) && (lev==0) &&
         (bdy_var_type != BdyVars::null) )
    {
        fill_from_bdyfiles (mf_to_fill,time,bdy_var_type);
    }
#endif

    // Enforce physical boundary conditions
    (*physbcs[lev])(mf_to_fill,0,ncomp,mf_to_fill.nGrowVect(),time,bccomp);

    // Also enforce free-slip at top boundary (on xvel or yvel)
    if ( (mf_box.ixType() == IndexType(IntVect(1,0,0))) ||
         (mf_box.ixType() == IndexType(IntVect(0,1,0))) )
    {
        int khi = geom[lev].Domain().bigEnd(2);
        for (MFIter mfi(mf_to_fill); mfi.isValid(); ++mfi)
        {
            Box gbx  = mfi.growntilebox(); // Note this is face-centered since vel is
            gbx.setSmall(2,khi+1);
            if (gbx.ok()) {
                Array4<Real> vel_arr = mf_to_fill.array(mfi);
                ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    vel_arr(i,j,k) = vel_arr(i,j,khi);
                });
            }
        }
    }
}


// utility to copy in data from old/new data into a struct that holds data for FillPatching
TimeInterpolatedData
ROMSX::GetDataAtTime (int /*lev*/, Real /*time*/)
{
    BL_PROFILE_VAR("GetDataAtTime()",GetDataAtTime);
    TimeInterpolatedData data;

// HACK HACK HACK
#if 0
    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
        for (int i = 0; i < Vars::NumTypes; ++i) {
            data.add_var(&vars_new[lev][i], data.non_owning);
        }
        data.set_time(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
        for (int i = 0; i < Vars::NumTypes; ++i) {
            data.add_var(&vars_old[lev][i], data.non_owning);
        }
        data.set_time(t_old[lev]);
    }
    else if (time > t_old[lev] && time < t_new[lev])
    {
        // do first order interpolation in time between [t_old[lev], t_new[lev]]
        // time interpolation includes the ghost cells
        for (int i = 0; i < Vars::NumTypes; ++i) {
            MultiFab* mf_tmp = new MultiFab(vars_new[lev][i].boxArray(),
                                            vars_new[lev][i].DistributionMap(),
                                            vars_new[lev][i].nComp(), vars_new[lev][i].nGrowVect());
            mf_tmp->setVal(0.0_rt);

            const Real dt_fraction = (time - t_old[lev]) / (t_new[lev] - t_old[lev]);
            MultiFab::Saxpy(*mf_tmp, 1.0_rt - dt_fraction, vars_old[lev][i], 0, 0, mf_tmp->nComp(), mf_tmp->nGrowVect());
            MultiFab::Saxpy(*mf_tmp,          dt_fraction, vars_new[lev][i], 0, 0, mf_tmp->nComp(), mf_tmp->nGrowVect());

            data.add_var(mf_tmp, data.owning);
        }
        data.set_time(time);
    }
    else
    {
        amrex::Error("Requested data at a time outside the interval [t_old, t_new]");
    }

    // We need to make sure to fill these before we compute the viscosity
    for (int i = 0; i < Vars::NumTypes; ++i) {
        data.get_var(Vars::xvel).FillBoundary(geom[lev].periodicity());
        data.get_var(Vars::yvel).FillBoundary(geom[lev].periodicity());
        data.get_var(Vars::zvel).FillBoundary(geom[lev].periodicity());
        data.get_var(Vars::cons).FillBoundary(geom[lev].periodicity());
    }
#endif
    return data;
}

// Fill an entire multifab by interpolating from the coarser level -- this is used
//     only when a new level of refinement is being created during a run (i.e not at initialization)
//     This will never be used with static refinement.
void
ROMSX::FillCoarsePatch (int lev, Real time, MultiFab* mf_to_fill, MultiFab* mf_crse)
{
    BL_PROFILE_VAR("FillCoarsePatch()",FillCoarsePatch);
    AMREX_ASSERT(lev > 0);

    int bccomp = 0;
    int  icomp = 0;
    int  ncomp = 1;
    amrex::Interpolater* mapper = nullptr;

    Box box_mf = ((*mf_to_fill)[0]).box();
    if (box_mf.ixType() == IndexType(IntVect(0,0,0)))
    {
        bccomp = 0;
        mapper = &cell_cons_interp;
        ncomp = NCONS;
    }
    else if (box_mf.ixType() == IndexType(IntVect(1,0,0)))
    {
        bccomp = BCVars::xvel_bc;
        mapper = &face_linear_interp;
    }
    else if (box_mf.ixType() == IndexType(IntVect(0,1,0)))
    {
        bccomp = BCVars::yvel_bc;
        mapper = &face_linear_interp;
    }
    else if (box_mf.ixType() == IndexType(IntVect(0,0,1)))
    {
        bccomp = BCVars::zvel_bc;
        mapper = &face_linear_interp;
    } else {
          amrex::Abort("Dont recognize this box type in ROMSX_FillPatch");
    }

#if 0
    TimeInterpolatedData cdata = GetDataAtTime(lev-1, time);
    TimeInterpolatedData fdata = GetDataAtTime(lev  , time);
    Vector<Real> ctime = {cdata.get_time()};

    Vector<MultiFab*> cmf = {mf_crse};
    Vector<Real> ctime = {time};

    ROMSXPhysBCFunct cphysbc(lev-1,geom[lev-1],
                             domain_bcs_type,domain_bcs_type_d,
                             cdata,
                             m_bc_extdir_vals
#ifdef ROMSX_USE_NETCDF
                            ,ic_bc_type,bdy_data_xlo,bdy_data_xhi,
                             bdy_data_ylo,bdy_data_yhi,bdy_time_interval
#endif
                            );
    ROMSXPhysBCFunct fphysbc(lev,geom[lev],
                             domain_bcs_type,domain_bcs_type_d,
                             fdata,
                             m_bc_extdir_vals
#ifdef ROMSX_USE_NETCDF
                            ,ic_bc_type,bdy_data_xlo,bdy_data_xhi,
                             bdy_data_ylo,bdy_data_yhi,bdy_time_interval
#endif
                            );
#endif

//  amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
//                               cphysbc, 0, fphysbc, 0, refRatio(lev-1),
//                               mapper, domain_bcs_type, bccomp);
    amrex::InterpFromCoarseLevel(*mf_to_fill, time, mf_crse[lev-1], 0, icomp, ncomp, geom[lev-1], geom[lev],
                                 null_bc, 0, null_bc, 0, refRatio(lev-1),
                                 mapper, domain_bcs_type, bccomp);
}
