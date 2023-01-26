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
ROMSX::FillPatch (int lev, Real time, const Vector<MultiFab*>& mfs)
{
    BL_PROFILE_VAR("ROMSX::FillPatch()",ROMSX_FillPatch);
    int bccomp;
    amrex::Interpolater* mapper = nullptr;

    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx) {
        MultiFab& mf = *mfs[var_idx];
        const int icomp = 0;
        const int ncomp = mf.nComp();

        if (var_idx == Vars::cons)
        {
            bccomp = 0;
            mapper = &cell_cons_interp;
        }
        else if (var_idx == Vars::xvel)
        {
            bccomp = BCVars::xvel_bc;
            mapper = &face_linear_interp;
        }
        else if (var_idx == Vars::yvel)
        {
            bccomp = BCVars::yvel_bc;
            mapper = &face_linear_interp;
        }
        else if (var_idx == Vars::zvel)
        {
            bccomp = BCVars::zvel_bc;
            mapper = &face_linear_interp;
        } else {
          amrex::Abort("Dont recognize this variable type in ROMSX_Fillpatch");
        }

        if (lev == 0)
        {
            Vector<MultiFab*> fmf = {&vars_old[lev][var_idx], &vars_new[lev][var_idx]};
            Vector<Real> ftime    = {t_old[lev], t_new[lev]};
            amrex::FillPatchSingleLevel(mf, time, fmf, ftime, icomp, icomp, ncomp,
                                        geom[lev], null_bc, bccomp);
        }
        else
        {
            Vector<MultiFab*> fmf = {&vars_old[lev][var_idx], &vars_new[lev][var_idx]};
            Vector<Real> ftime    = {t_old[lev], t_new[lev]};
            Vector<MultiFab*> cmf = {&vars_old[lev-1][var_idx], &vars_new[lev-1][var_idx]};
            Vector<Real> ctime    = {t_old[lev-1], t_new[lev-1]};

            amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                      0, icomp, ncomp, geom[lev-1], geom[lev],
                                      null_bc, bccomp, null_bc, bccomp, refRatio(lev-1),
                                      mapper, domain_bcs_type, bccomp);
        } // lev > 0
    } // var_idx

    // ***************************************************************************
    // Physical bc's at domain boundary
    // ***************************************************************************
    bool cons_only = false;
    int icomp_cons = 0;
    int ncomp_cons = mfs[Vars::cons]->nComp();

    IntVect ngvect_cons = mfs[Vars::cons]->nGrowVect();
    IntVect ngvect_vels = mfs[Vars::xvel]->nGrowVect();
    //tweaked physbcs
    for(auto& mf : mfs)
    {
	amrex::Abort("Need to initialize physbcs");
    (*physbcs[lev])(*mf,icomp_cons,ncomp_cons,ngvect_cons,time,cons_only);
    }
    /*
    if (m_r2d) amrex::Abort("ReadBoundaryPlanes is not supported");//fill_from_bndryregs(mfs,time);
#ifdef ROMSX_USE_NETCDF
    if (init_type == "real") amrex::Abort("This init type is not supported");//fill_from_wrfbdy(mfs,time);
#endif
    */
}

// utility to copy in data from old/new data into a struct that holds data for FillPatching
TimeInterpolatedData
ROMSX::GetDataAtTime (int lev, Real time)
{
    BL_PROFILE_VAR("GetDataAtTime()",GetDataAtTime);
    TimeInterpolatedData data;

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

    return data;
}
//
// Fill valid and ghost data in the MultiFab "mf"
// This version fills the MultiFab mf in valid regions with the "state data" at the given time;
// values in mf when it is passed in are *not* used.
//
void
ROMSX::FillPatch (int lev, Real time, Real time_mt, Real delta_t, Vector<MultiFab>& mfs)
{
    BL_PROFILE_VAR("ROMSX::FillPatch()",ROMSX_FillPatch);
    int bccomp;
    amrex::Interpolater* mapper = nullptr;

    TimeInterpolatedData fdata = GetDataAtTime(lev, time);
    Vector<Real> ftime         = {fdata.get_time()};

    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx) {
        MultiFab& mf = mfs[var_idx];
        const int icomp = 0;
        const int ncomp = mf.nComp();

        if (var_idx == Vars::cons)
        {
            bccomp = 0;
            mapper = &cell_cons_interp;
        }
        else if (var_idx == Vars::xvel)
        {
            bccomp = BCVars::xvel_bc;
            mapper = &face_linear_interp;
        }
        else if (var_idx == Vars::yvel)
        {
            bccomp = BCVars::yvel_bc;
            mapper = &face_linear_interp;
        }
        else if (var_idx == Vars::zvel)
        {
            bccomp = BCVars::zvel_bc;
            mapper = &face_linear_interp;
        } else {
          amrex::Abort("Dont recognize this variable type in ROMSX_Fillpatch");
        }

        if (lev == 0)
        {
            Vector<MultiFab*> smf = {&fdata.get_var(var_idx)};
            ROMSXPhysBCFunct physbc(lev,time_mt,delta_t,geom[lev],
                                  domain_bcs_type,domain_bcs_type_d,
                                  var_idx,solverChoice.terrain_type,
                                  fdata,m_bc_extdir_vals,z_phys_nd[lev], detJ_cc[lev] 
#ifdef ROMSX_USE_NETCDF
                                 ,init_type,bdy_data_xlo,bdy_data_xhi,
                                  bdy_data_ylo,bdy_data_yhi,bdy_time_interval 
#endif
                                  );
            amrex::FillPatchSingleLevel(mf, time, smf, ftime, 0, icomp, ncomp,
                                        geom[lev], physbc, bccomp);
        }
        else
        {
            TimeInterpolatedData cdata = GetDataAtTime(lev-1, time);
            Vector<Real> ctime = {cdata.get_time()};
            Vector<MultiFab*> cmf = {&cdata.get_var(var_idx)};
            Vector<MultiFab*> fmf = {&fdata.get_var(var_idx)};

            ROMSXPhysBCFunct cphysbc(lev-1,time_mt,delta_t,geom[lev-1],
                                   domain_bcs_type,domain_bcs_type_d,
                                   var_idx,solverChoice.terrain_type,cdata,
                                   m_bc_extdir_vals,z_phys_nd[lev-1],detJ_cc[lev-1] 
#ifdef ROMSX_USE_NETCDF
                                  ,init_type,bdy_data_xlo,bdy_data_xhi,
                                   bdy_data_ylo,bdy_data_yhi,bdy_time_interval 
#endif
                                   );
            ROMSXPhysBCFunct fphysbc(lev,time_mt,delta_t,geom[lev],
                                   domain_bcs_type,domain_bcs_type_d,
                                   var_idx,solverChoice.terrain_type,fdata,
                                   m_bc_extdir_vals,z_phys_nd[lev], detJ_cc[lev] 
#ifdef ROMSX_USE_NETCDF
                                  ,init_type,bdy_data_xlo,bdy_data_xhi,
                                   bdy_data_ylo,bdy_data_yhi,bdy_time_interval 
#endif
                                   );

            amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                      0, icomp, ncomp, geom[lev-1], geom[lev],
                                      cphysbc, bccomp, fphysbc, bccomp, refRatio(lev-1),
                                      mapper, domain_bcs_type, bccomp);
        } // lev > 0
    } // var_idx
}

//
// Fill valid and ghost data in the MultiFabs in "mfs"
// mfs is a Vector<std::reference_wrapper<MultiFab> > containing, in order: cons, xvel, yvel, and zvel data
// This version fills the MultiFabs mfs in valid regions with the values in "mfs" when it is passed in;
// it is used only to compute ghost values for intermediate stages of a time integrator.
//
void
ROMSX::FillIntermediatePatch (int lev, Real time, Real time_mt, Real delta_t,
                            Vector<std::reference_wrapper<MultiFab> > mfs,
                            int ng_cons, int ng_vel, bool cons_only, int scomp_cons, int ncomp_cons)
{
    BL_PROFILE_VAR("FillIntermediatePatch()",FillIntermediatePatch);
    int bccomp;
    amrex::Interpolater* mapper;
    TimeInterpolatedData level_data;

    // We should always pass cons, xvel, yvel, and zvel (in that order) in the mfs vector
    AMREX_ALWAYS_ASSERT(mfs.size() == Vars::NumTypes);

    for (int imf = 0; imf < Vars::NumTypes; ++imf) {
        level_data.add_var(&mfs[imf].get(), level_data.non_owning);
    }
    level_data.set_time(time);

    for (int var_idx = 0; var_idx < Vars::NumTypes; ++var_idx)
    {
        if (cons_only && var_idx != Vars::cons) continue;

        MultiFab& mf = mfs[var_idx].get();

        IntVect ngvect;
        int icomp, ncomp;
        if (var_idx == Vars::cons)
        {
            bccomp = 0;
            mapper = &cell_cons_interp;
            ngvect = IntVect(ng_cons,ng_cons,ng_cons);
            icomp  = scomp_cons;
            ncomp  = ncomp_cons;
        }
        else if (var_idx == Vars::xvel)
        {
            bccomp = NVAR;
            mapper = &face_linear_interp;
            ngvect = IntVect(ng_vel,ng_vel,ng_vel);
            icomp  = 0;
            ncomp  = 1;
        }
        else if (var_idx == Vars::yvel)
        {
            bccomp = NVAR+1;
            mapper = &face_linear_interp;
            ngvect = IntVect(ng_vel,ng_vel,ng_vel);
            icomp  = 0;
            ncomp  = 1;
        }
        else if (var_idx == Vars::zvel)
        {
            bccomp = NVAR+2;
            mapper = &face_linear_interp;
            ngvect = IntVect(ng_vel,ng_vel,0);
            icomp  = 0;
            ncomp  = 1;
        }

        if (lev == 0)
        {
            // on lev, use the mf data and time passed to FillIntermediatePatch().
            Vector<MultiFab*> smf { &mf };
            Vector<Real> stime { time };

            ROMSXPhysBCFunct physbc(lev,time_mt,delta_t,geom[lev],
                                  domain_bcs_type,domain_bcs_type_d,
                                  var_idx,solverChoice.terrain_type,level_data,
                                  m_bc_extdir_vals,z_phys_nd[lev],detJ_cc[lev] 
#ifdef ROMSX_USE_NETCDF
                                 ,init_type,bdy_data_xlo,bdy_data_xhi,
                                  bdy_data_ylo,bdy_data_yhi,bdy_time_interval 
#endif
                                  );

            amrex::FillPatchSingleLevel(mf, ngvect, time, smf, stime, 0, icomp, ncomp,
                                        geom[lev], physbc, bccomp);
        }
        else
        {
            MultiFab mf_tmp(mf.boxArray(), mf.DistributionMap(), ncomp, mf.nGrowVect());

            TimeInterpolatedData cdata = GetDataAtTime(lev-1, time);
            Vector<MultiFab*> cmf = {&cdata.get_var(var_idx)};
            Vector<MultiFab*> fmf = {&mf};
            Vector<Real> ctime = {cdata.get_time()};
            Vector<Real> ftime = {level_data.get_time()};

            ROMSXPhysBCFunct cphysbc(lev-1,time_mt,delta_t,geom[lev-1],
                                   domain_bcs_type,domain_bcs_type_d,
                                   var_idx,solverChoice.terrain_type,cdata,
                                   m_bc_extdir_vals,z_phys_nd[lev-1],detJ_cc[lev-1] 
#ifdef ROMSX_USE_NETCDF
                                  ,init_type,bdy_data_xlo,bdy_data_xhi,
                                   bdy_data_ylo,bdy_data_yhi,bdy_time_interval 
#endif
                                   );
            ROMSXPhysBCFunct fphysbc(lev,time_mt,delta_t,geom[lev],
                                   domain_bcs_type,domain_bcs_type_d,
                                   var_idx,solverChoice.terrain_type,level_data,
                                   m_bc_extdir_vals,z_phys_nd[lev],detJ_cc[lev] 
#ifdef ROMSX_USE_NETCDF
                                  ,init_type,bdy_data_xlo,bdy_data_xhi,
                                   bdy_data_ylo,bdy_data_yhi,bdy_time_interval 
#endif
                                   );

            amrex::FillPatchTwoLevels(mf_tmp, ngvect, time, cmf, ctime, fmf, ftime,
                                      0, icomp, ncomp, geom[lev-1], geom[lev],
                                      cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                      mapper, domain_bcs_type, bccomp);

            // Replace mf with mf_tmp
            if (ncomp == mf.nComp())
                std::swap(mf_tmp, mf);
            else
                MultiFab::Copy(mf,mf_tmp,0,0,1,mf.nGrowVect());
        }
    }
}

// Fill an entire multifab by interpolating from the coarser level -- this is used
//     only when a new level of refinement is being created during a run (i.e not at initialization)
//     This will never be used with static refinement.
void
ROMSX::FillCoarsePatch (int lev, Real time, Real time_mt, Real delta_t,
                      MultiFab& mf, int icomp, int ncomp, int var_idx)
{
    BL_PROFILE_VAR("FillCoarsePatch()",FillCoarsePatch);
    AMREX_ASSERT(lev > 0);

    int bccomp;
    amrex::Interpolater* mapper;

    if (var_idx == Vars::cons)
    {
        bccomp = 0;
        mapper = &cell_cons_interp;
    }
    else if (var_idx == Vars::xvel || var_idx == Vars::xmom)
    {
        bccomp = NVAR;
        mapper = &face_linear_interp;
    }
    else if (var_idx == Vars::yvel || var_idx == Vars::ymom)
    {
        bccomp = NVAR+1;
        mapper = &face_linear_interp;
    }
    else if (var_idx == Vars::zvel || var_idx == Vars::zmom)
    {
        bccomp = NVAR+2;
        mapper = &face_linear_interp;
    }

    TimeInterpolatedData cdata = GetDataAtTime(lev-1, time);
    TimeInterpolatedData fdata = GetDataAtTime(lev  , time);
    Vector<MultiFab*> cmf = {&cdata.get_var(var_idx)};
    Vector<MultiFab*> fmf = {&fdata.get_var(var_idx)};
    Vector<Real> ctime = {cdata.get_time()};
    Vector<Real> ftime = {fdata.get_time()};

    ROMSXPhysBCFunct cphysbc(lev-1,time_mt,delta_t,geom[lev-1],
                           domain_bcs_type,domain_bcs_type_d,
                           var_idx,solverChoice.terrain_type,cdata,
                           m_bc_extdir_vals,z_phys_nd[lev-1],detJ_cc[lev-1] 
#ifdef ROMSX_USE_NETCDF
                          ,init_type,bdy_data_xlo,bdy_data_xhi,
                           bdy_data_ylo,bdy_data_yhi,bdy_time_interval
#endif
                           );
    ROMSXPhysBCFunct fphysbc(lev,time_mt,delta_t,geom[lev],
                           domain_bcs_type,domain_bcs_type_d,
                           var_idx,solverChoice.terrain_type,fdata,
                           m_bc_extdir_vals,z_phys_nd[lev],detJ_cc[lev] 
#ifdef ROMSX_USE_NETCDF
                          ,init_type,bdy_data_xlo,bdy_data_xhi,
                           bdy_data_ylo,bdy_data_yhi,bdy_time_interval 
#endif
                           );

    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
                                    cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                    mapper, domain_bcs_type, bccomp);
}

void
ROMSX::FillCoarsePatchAllVars (int lev, Real time, Real time_mt, Real delta_t, Vector<MultiFab>& vmf)
{
    for (int var_idx = 0; var_idx < vmf.size(); ++var_idx) {
        FillCoarsePatch(lev, time, time_mt, delta_t, vmf[var_idx], 0, vmf[var_idx].nComp(), var_idx);
    }
}
