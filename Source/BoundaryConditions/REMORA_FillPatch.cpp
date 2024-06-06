#include <REMORA.H>
#include <REMORA_PhysBCFunct.H>
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
REMORA::FillPatch (int lev, Real time, MultiFab& mf_to_fill, Vector<MultiFab*> const& mfs,
#ifdef REMORA_USE_NETCDF
                  const int bdy_var_type,
#else
                  const int /*bdy_var_type*/,
#endif
                  const int  icomp,
                  const bool fill_all,
                  const bool fill_set)
{
    BL_PROFILE_VAR("REMORA::FillPatch()",REMORA_FillPatch);
    int bccomp;
    amrex::Interpolater* mapper = nullptr;

    Box mf_box(mf_to_fill.boxArray()[0]);

    //
    // ***************************************************************************
    // The first thing we do is interpolate the momenta on the "valid" faces of
    // the fine grids (where the interface is coarse/fine not fine/fine) -- this
    // will not be over-written below because the FillPatch operators see these as
    // valid faces.
    // ***************************************************************************
    if (lev>0 && fill_set) {
        if (cf_set_width > 0 &&
            mf_box.ixType() == IndexType(IntVect(0,0,0))) {
            FPr_c[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
        } else if (fill_all && cf_set_width >= 0) {
            if (mf_box.ixType() == IndexType(IntVect(1,0,0))) {
                FPr_u[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
            } else if (mf_box.ixType() == IndexType(IntVect(0,1,0))) {
                FPr_v[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
            } else if (mf_box.ixType() == IndexType(IntVect(0,0,1))) {
                FPr_w[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
            }
        }
    }

    int ncomp;
    if (fill_all) {
        ncomp = mf_to_fill.nComp();
    } else {
        ncomp = 1;
    }

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
    else {
        bccomp = BCVars::zvel_bc;
        mapper = &face_linear_interp;
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

    // This is currently unconditionally true, but do_bc flag should actually be passed in
    const bool do_bc=true;
    if (do_bc) {
        // ***************************************************************************
        // Physical bc's at domain boundary
        // ***************************************************************************

        // Enforce physical boundary conditions
        (*physbcs[lev])(mf_to_fill,icomp,ncomp,mf_to_fill.nGrowVect(),time,bccomp);

#ifdef REMORA_USE_NETCDF
        // Fill the data which is stored in the boundary data read from netcdf files
        if ( (solverChoice.ic_bc_type == IC_BC_Type::Real) && (lev==0) &&
             (bdy_var_type != BdyVars::null) )
        {
            fill_from_bdyfiles (mf_to_fill,time,bdy_var_type, icomp);
        }
#endif

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
}
//
// Fill valid and ghost data in the MultiFab "mf"
// This version fills the MultiFab mf in valid regions with the "state data" at the given time;
// values in mf when it is passed in are *not* used.
// Unlike FillPatch, FillPatchNoBC does not apply boundary conditions.
//
void
REMORA::FillPatchNoBC (int lev, Real time, MultiFab& mf_to_fill, Vector<MultiFab*> const& mfs,
#ifdef REMORA_USE_NETCDF
                  const int bdy_var_type,
#else
                  const int /*bdy_var_type*/,
#endif
                  const int  icomp,
                  const bool fill_all,
                  const bool fill_set)
{
    // HACK: Note that this is hacky; should be able to have a single call to FillPatch with a
    // flag for bcs, but for some reason it was acting weird, so we're splitting this out into
    // two functions with repeated code for the time being
    BL_PROFILE_VAR("REMORA::FillPatch()",REMORA_FillPatch);
    int bccomp;
    amrex::Interpolater* mapper = nullptr;

    Box mf_box(mf_to_fill.boxArray()[0]);

    //
    // ***************************************************************************
    // The first thing we do is interpolate the momenta on the "valid" faces of
    // the fine grids (where the interface is coarse/fine not fine/fine) -- this
    // will not be over-written below because the FillPatch operators see these as
    // valid faces.
    // ***************************************************************************
    if (lev>0 && fill_set) {
        if (cf_set_width > 0 &&
            mf_box.ixType() == IndexType(IntVect(0,0,0))) {
            FPr_c[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
        } else if (fill_all && cf_set_width >= 0) {
            if (mf_box.ixType() == IndexType(IntVect(1,0,0))) {
                FPr_u[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
            } else if (mf_box.ixType() == IndexType(IntVect(0,1,0))) {
                FPr_v[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
            } else if (mf_box.ixType() == IndexType(IntVect(0,0,1))) {
                FPr_w[lev-1].FillSet(mf_to_fill, time, null_bc, domain_bcs_type);
            }
        }
    }

    int ncomp;
    if (fill_all) {
        ncomp = mf_to_fill.nComp();
    } else {
        ncomp = 1;
    }

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
    else {
        bccomp = BCVars::zvel_bc;
        mapper = &face_linear_interp;
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
}


// utility to copy in data from old/new data into a struct that holds data for FillPatching
TimeInterpolatedData
REMORA::GetDataAtTime (int /*lev*/, Real /*time*/)
{
    BL_PROFILE_VAR("GetDataAtTime()",GetDataAtTime);
    TimeInterpolatedData data;

// HACK HACK HACK
#if 0
    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3_rt;

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
REMORA::FillCoarsePatch (int lev, Real time, MultiFab* mf_to_fill, MultiFab* mf_crse,
                         const int  icomp,
                         const bool fill_all)
{
    BL_PROFILE_VAR("FillCoarsePatch()",FillCoarsePatch);
    AMREX_ASSERT(lev > 0);

    int ncomp;
    if (fill_all) {
        ncomp = mf_to_fill->nComp();
    } else {
        ncomp = 1;
    }

    int bccomp = 0;
    amrex::Interpolater* mapper = nullptr;

    Box box_mf(mf_to_fill->boxArray()[0]);

    if (box_mf.ixType() == IndexType(IntVect(0,0,0)))
    {
        bccomp = 0;
        mapper = &cell_cons_interp;
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
          amrex::Abort("Dont recognize this box type in REMORA_FillPatch");
    }

#if 0
    TimeInterpolatedData cdata = GetDataAtTime(lev-1, time);
    TimeInterpolatedData fdata = GetDataAtTime(lev  , time);
    Vector<Real> ctime = {cdata.get_time()};

    Vector<MultiFab*> cmf = {mf_crse};
    Vector<Real> ctime = {time};

    REMORAPhysBCFunct cphysbc(lev-1,geom[lev-1],
                             domain_bcs_type,domain_bcs_type_d,
                             cdata,
                             m_bc_extdir_vals
#ifdef REMORA_USE_NETCDF
                            ,ic_bc_type,bdy_data_xlo,bdy_data_xhi,
                             bdy_data_ylo,bdy_data_yhi,bdy_time_interval
#endif
                            );
    REMORAPhysBCFunct fphysbc(lev,geom[lev],
                             domain_bcs_type,domain_bcs_type_d,
                             fdata,
                             m_bc_extdir_vals
#ifdef REMORA_USE_NETCDF
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

void
REMORA::FillBdyCCVels (int lev, MultiFab& mf_cc_vel)
{
    // Impose bc's at domain boundaries
//    for (int lev = 0; lev <= finest_level; ++lev)
//    {
    Box domain(Geom(lev).Domain());

    int ihi = domain.bigEnd(0);
    int jhi = domain.bigEnd(1);
    int khi = domain.bigEnd(2);

    // Impose periodicity first
    mf_cc_vel.FillBoundary(geom[lev].periodicity());

    for (MFIter mfi(mf_cc_vel, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Note that we don't fill corners here -- only the cells that share a face
        //      with interior cells -- this is all that is needed to calculate vorticity
        const Box& bx = mfi.tilebox();
        const Array4<Real>& vel_arr = mf_cc_vel.array(mfi);

        if (!Geom(lev).isPeriodic(0)) {
            // Low-x side
            if (bx.smallEnd(0) <= domain.smallEnd(0)) {
                Real mult = (phys_bc_type[0] == REMORA_BC::no_slip_wall) ? -1. : 1.;
                ParallelFor(makeSlab(bx,0,0), [=] AMREX_GPU_DEVICE(int , int j, int k) noexcept
                {
                    vel_arr(-1,j,k,1) = mult*vel_arr(0,j,k,1); // v
                    vel_arr(-1,j,k,2) = mult*vel_arr(0,j,k,2); // w
                });
            }

            // High-x side
            if (bx.bigEnd(0) >= domain.bigEnd(0)) {
                Real mult = (phys_bc_type[3] == REMORA_BC::no_slip_wall) ? -1. : 1.;
                ParallelFor(makeSlab(bx,0,0), [=] AMREX_GPU_DEVICE(int , int j, int k) noexcept
                {
                    vel_arr(ihi+1,j,k,1) = mult*vel_arr(ihi,j,k,1); // v
                    vel_arr(ihi+1,j,k,2) = mult*vel_arr(ihi,j,k,2); // w
                });
            }
        } // !periodic

        if (!Geom(lev).isPeriodic(1)) {
            // Low-y side
            if (bx.smallEnd(1) <= domain.smallEnd(1)) {
                Real mult = (phys_bc_type[1] == REMORA_BC::no_slip_wall) ? -1. : 1.;
                ParallelFor(makeSlab(bx,1,0), [=] AMREX_GPU_DEVICE(int i, int  , int k) noexcept
                {
                    vel_arr(i,-1,k,0) = mult*vel_arr(i,0,k,0); // u
                    vel_arr(i,-1,k,2) = mult*vel_arr(i,0,k,2); // w
                });
            }

            // High-y side
            if (bx.bigEnd(1) >= domain.bigEnd(1)) {
                Real mult = (phys_bc_type[4] == REMORA_BC::no_slip_wall) ? -1. : 1.;
                ParallelFor(makeSlab(bx,1,0), [=] AMREX_GPU_DEVICE(int i, int , int k) noexcept
                {
                    vel_arr(i,jhi+1,k,0) = mult*vel_arr(i,jhi,k,0); // u
                    vel_arr(i,jhi+1,k,2) = mult*-vel_arr(i,jhi,k,2); // w
                });
            }
        } // !periodic

        if (!Geom(lev).isPeriodic(2)) {
            // Low-z side
            if (bx.smallEnd(2) <= domain.smallEnd(2)) {
                Real mult = (phys_bc_type[2] == REMORA_BC::no_slip_wall) ? -1. : 1.;
                ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
                {
                    vel_arr(i,j,-1,0) = mult*vel_arr(i,j,0,0); // u
                    vel_arr(i,j,-1,1) = mult*vel_arr(i,j,0,1); // v
                });
            }

            // High-z side
            if (bx.bigEnd(2) >= domain.bigEnd(2)) {
                Real mult = (phys_bc_type[5] == REMORA_BC::no_slip_wall) ? -1. : 1.;
                ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept
                {
                    vel_arr(i,j,khi+1,0) = mult*vel_arr(i,j,khi,0); // u
                    vel_arr(i,j,khi+1,1) = mult*vel_arr(i,j,khi,1); // v
                });
            }
        } // !periodic
    } // MFIter

//    } // lev
}
