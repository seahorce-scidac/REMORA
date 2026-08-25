/**
 * \file REMORA_init.cpp
 */

#include <REMORA.H>
#include <REMORA_Constants.H>
#include <REMORA_prob_common.H>

using namespace amrex;

/**
 * @param[in   ] lev     level to initialize on
 */
void
REMORA::init_analytic(int lev)
{
    prob->init_analytic_prob(lev, geom[lev], solverChoice, *this, *cons_new[lev], *xvel_new[lev], *yvel_new[lev]);

    set_grid_scale(lev);
}

/**
 * \brief Initialize the biology tracers from whichever source
 *        remora.biology_ic_type selects.
 *
 * Deliberately decoupled from remora.ic_type. A NetCDF initial file supplies
 * only the tracers it happens to contain -- the BioToy file has `oxygen` but
 * no `PO4` or `ODU` -- so validating an option set the file does not cover
 * would otherwise mean regenerating binary input that no fresh clone can
 * reproduce. With biology_ic_type = analytic the physical fields still come
 * from NetCDF while biology comes from the problem's analytic profile.
 *
 * Must be called after the physical fields are initialized: the analytic
 * biology profiles are functions of temperature.
 *
 * @param[in   ] lev     level to initialize on
 */
void
REMORA::init_biology_ic (int lev)
{
    if (!REMORABiology::has_biology(biology_model)) {
        return;
    }

    const bool use_analytic =
        (biology_ic_type == REMORABiology::BiologyICType::analytic) ||
        (biology_ic_type == REMORABiology::BiologyICType::follow_ic_type &&
         solverChoice.ic_type != IC_Type::netcdf);

    if (use_analytic) {
        // The base-class hook in REMORA_prob_common.H is an amrex::Error, so a
        // problem that enables biology without providing an analytic biology
        // IC fails loudly rather than starting from uninitialized tracers.
        prob->init_analytic_biology(lev, geom[lev], solverChoice, fennel_params,
                                    *this, *cons_new[lev]);
    } else {
#ifdef REMORA_USE_NETCDF
        init_biology_from_netcdf(lev);
#else
        amrex::Abort("remora.biology_ic_type = netcdf requires a NetCDF build");
#endif
    }
}

/**
 * \brief Full-domain counterpart of init_biology_ic, for the hires_init_level
 *        path where initial data is specified on a high-resolution level and
 *        averaged down.
 *
 * Fills the biology components of vec_cons_full_domain[hires_init_level].
 * Must be called after the physical fields on that level are set (the
 * analytic profiles are functions of temperature) and before whatever
 * averages vec_cons_full_domain down.
 */
void
REMORA::init_biology_ic_full_domain ()
{
    if (!REMORABiology::has_biology(biology_model)) {
        return;
    }

    const bool use_analytic =
        (biology_ic_type == REMORABiology::BiologyICType::analytic) ||
        (biology_ic_type == REMORABiology::BiologyICType::follow_ic_type &&
         solverChoice.ic_type != IC_Type::netcdf);

    if (use_analytic) {
        prob->init_analytic_biology(hires_init_level, geom[hires_init_level],
                                    solverChoice, fennel_params, *this,
                                    *vec_cons_full_domain[hires_init_level]);
    } else {
#ifdef REMORA_USE_NETCDF
        init_biology_full_domain_from_netcdf();
#else
        amrex::Abort("remora.biology_ic_type = netcdf requires a NetCDF build");
#endif
    }

    // Average down here, after the branch, rather than inside one of the two sources. When
    // this loop lived in the NetCDF reader, the analytic branch filled hires_init_level and
    // nothing else, so every level below it -- including level 0, the one the run actually
    // integrates -- kept the zeros that init_data_full_domain_from_netcdf had written.
    for (int lev = hires_init_level-1; lev >= 0; lev--) {
        average_down_with_grow_cells(lev, vec_cons_full_domain);
    }
}

/**
 * @param[in   ] lev     level to initialize on
 */
void
REMORA::init_beta_plane_coriolis (int lev)
{
    std::unique_ptr<MultiFab>& mf_fcor = vec_fcor[lev];
    auto geomdata  = Geom(lev).data();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*cons_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Valid region only, since vec_yr's ghost cells aren't all filled. set_coriolis
        // FillPatches vec_fcor with foextrap right after this, as for coriolis_type = netcdf.
        Box bx = mfi.tilebox();
        auto fcor_arr = (mf_fcor)->array(mfi);
        auto yr_arr = vec_yr[lev]->const_array(mfi);
        Real coriolis_f0 = solverChoice.coriolis_f0;
        Real coriolis_beta = solverChoice.coriolis_beta;
        // Units: yr is metric (1/pn, or y_rho from the grid file), so beta_plane assumes
        // prob_lo/prob_hi are too. A domain given in degrees would silently mix the two.
        Real Esize = geomdata.ProbHi()[1] - geomdata.ProbLo()[1];

        ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            // yr is measured from the southern edge of the domain (ROMS ana_grid convention), so
            // f is coriolis_f0 at mid-domain. Adding prob_lo back would shift f by beta*prob_lo.
            fcor_arr(i,j,0) = coriolis_f0 + coriolis_beta * (yr_arr(i,j,0) - Real(0.5) * Esize);
        });
    } //mfi

    vec_fcor[lev]->FillBoundary(geom[lev].periodicity());
}

/**
 * @param[in   ] lev     level to operate on
 */
void
REMORA::set_zeta_average (int lev)
{
    std::unique_ptr<MultiFab>& mf_zeta = vec_zeta[lev];
    std::unique_ptr<MultiFab>& mf_Zt_avg1  = vec_Zt_avg1[lev];
    for ( MFIter mfi(*cons_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        int nstp = 0;
        Array4<Real> const& Zt_avg1 = (mf_Zt_avg1)->array(mfi);
        Array4<const Real> const& zeta     = mf_zeta->const_array(mfi);

        Box  bx3 = mfi.tilebox()      ;  bx3.grow(IntVect(NGROW+1,NGROW+1,0)); //   cell-centered, grown by 3

        ParallelFor(makeSlab(bx3,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Zt_avg1(i,j,0) = zeta(i,j,0,nstp);
        });

    }
}

/**
 * @param[in   ] lev     level to operate on
 */
void
REMORA::set_2darrays (int lev)
{
    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

    vec_ubar[lev]->setVal(zero);
    vec_vbar[lev]->setVal(zero);

    MultiFab* U_old = xvel_new[lev];
    MultiFab* V_old = yvel_new[lev];
    std::unique_ptr<MultiFab>& mf_ubar = vec_ubar[lev];
    std::unique_ptr<MultiFab>& mf_vbar = vec_vbar[lev];
    std::unique_ptr<MultiFab>& mf_Hz  = vec_Hz[lev];
    int nstp = 0;

    for ( MFIter mfi(*cons_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& ubar = (mf_ubar)->array(mfi);
        Array4<Real> const& vbar = (mf_vbar)->array(mfi);

        Array4<const Real> const& Hz       = mf_Hz->const_array(mfi);
        Array4<const Real> const& u        = U_old->const_array(mfi);
        Array4<const Real> const& v        = V_old->const_array(mfi);

        Box  bx2 = mfi.tilebox()      ;  bx2.grow(IntVect(NGROW  ,NGROW  ,0)); //   cell-centered, grown by 2
        Box ubx2 = mfi.nodaltilebox(0); ubx2.grow(IntVect(NGROW  ,NGROW  ,0)); // x-face-centered, grown by 2
        Box vbx2 = mfi.nodaltilebox(1); vbx2.grow(IntVect(NGROW  ,NGROW  ,0)); // y-face-centered, grown by 2

        ParallelFor(makeSlab(ubx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Real CF = zero;
            Real sum_of_hz = zero;

            for (int k=0; k<=N; k++) {
                Real avg_hz = Real(0.5)*(Hz(i,j,k)+Hz(i-1,j,k));
                sum_of_hz += avg_hz;
                CF += avg_hz*u(i,j,k,nstp);
            }
            ubar(i,j,0,0) = CF / sum_of_hz;
        });

        ParallelFor(makeSlab(vbx2,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Real CF = zero;
            Real sum_of_hz = zero;

            for(int k=0; k<=N; k++) {
                Real avg_hz = Real(0.5)*(Hz(i,j,k)+Hz(i,j-1,k));
                sum_of_hz += avg_hz;
                CF += avg_hz*v(i,j,k,nstp);
            }
            vbar(i,j,0,0) = CF / sum_of_hz;
        });
    }

    FillPatch(lev, t_new[lev], *vec_ubar[lev], GetVecOfPtrs(vec_ubar), ubar_bc(), bdy_ubar(),0,false,false,0,0,zero,*vec_ubar[lev]);
    FillPatch(lev, t_new[lev], *vec_vbar[lev], GetVecOfPtrs(vec_vbar), vbar_bc(), bdy_vbar(),0,false,false,0,0,zero,*vec_vbar[lev]);
}

/**
 * @param[in   ] lev            level to operate on
 * @param[in   ] solver_choice  algorithmic choices
 */
void
REMORA::init_gls_vmix (int lev, SolverChoice solver_choice)
{
    vec_tke[lev]->setVal(solver_choice.gls_Kmin);
    vec_gls[lev]->setVal(solver_choice.gls_Pmin);
    vec_Lscale[lev]->setVal(zero);
    vec_Akk[lev]->setVal(solver_choice.Akk_bak);
    vec_Akp[lev]->setVal(solver_choice.Akp_bak);
    vec_Akv[lev]->setVal(solver_choice.Akv_bak);
    for (int n = 0; n < NAT; n++) {
        vec_Akt[lev]->setVal(solver_choice.Akt_bak[n], n, 1);
    }

    auto N = Geom(lev).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vec_Akk[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box bx = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        Array4<Real> const& Akk = vec_Akk[lev]->array(mfi);
        Array4<Real> const& Akp = vec_Akp[lev]->array(mfi);
        Array4<Real> const& Akt = vec_Akt[lev]->array(mfi);
        Array4<Real> const& Akv = vec_Akv[lev]->array(mfi);

        ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Akk(i,j, 0) = zero;
            Akk(i,j, N+1) = zero;

            Akp(i,j, 0) = zero;
            Akp(i,j, N+1) = zero;

            Akv(i,j, 0) = zero;
            Akv(i,j, N+1) = zero;

            Akt(i,j, 0) = zero;
            Akt(i,j, N+1) = zero;
        });
    }
}

/**
 * @param[in   ] lev     level to operate on
 */
void
REMORA::init_clim_nudg_coeff (int lev) {
    // Fill nudging coefficients from constant values, then overwrite
    // with coeffs read from file if using
    vec_nudg_coeff[BdyVars::u][lev]->setVal(solverChoice.nudg_coeff[BdyVars::u]);
    vec_nudg_coeff[BdyVars::v][lev]->setVal(solverChoice.nudg_coeff[BdyVars::v]);
    for (int icomp = 0; icomp < ncons; ++icomp) {
        vec_nudg_coeff[BdyVars::cons(icomp)][lev]->setVal(solverChoice.nudg_coeff[BdyVars::cons(icomp)]);
    }
    vec_nudg_coeff[bdy_ubar()][lev]->setVal(solverChoice.nudg_coeff[bdy_ubar()]);
    vec_nudg_coeff[bdy_vbar()][lev]->setVal(solverChoice.nudg_coeff[bdy_vbar()]);
    vec_nudg_coeff[bdy_zeta()][lev]->setVal(solverChoice.nudg_coeff[bdy_zeta()]);
#ifdef REMORA_USE_NETCDF
    if (solverChoice.do_any_clim_nudg) {
        // A coefficient file is optional. Without one every variable keeps the constant
        // timescale set above, which is a reasonable setup for tracers nudged on the
        // single remora.tnudg timescale.
        if (nc_clim_coeff_file.empty()) {
            amrex::Print() << "No remora.nc_clim_coeff_file given; climatology nudging will use "
                              "the constant timescales from remora.tnudg, m2nudg, and m3nudg"
                           << std::endl;
        } else {
            amrex::Print() << "Calling init_clim_nudg_coeff_from_netcdf \n " << std::endl;
            init_clim_nudg_coeff_from_netcdf(lev);
            amrex::Print() << "Climatology weights loaded from netcdf file \n " << std::endl;
        }
    }
#endif
}

void
REMORA::init_stretch_coeffs () {
    int nz = geom[0].Domain().length(2);
    s_r.resize(nz);
    s_w.resize(nz+1);
    Cs_r.resize(nz);
    Cs_w.resize(nz+1);

    calc_stretch_coeffs();
}

void REMORA::allocate_bathymetry_grid_vars_full_domain () {
    // Make fake boxArray that covers the whole domain on level 0
    BoxArray ba;
    ba.define(makeSlab(geom[0].Domain(),2,0));
    Box refined_domain = makeSlab(geom[0].Domain(),2,0);

    DistributionMapping dm(ba);
    vec_h_full_domain[0].reset(new MultiFab(ba, dm, 1, IntVect(1,1,0)));
    vec_pm_full_domain[0].reset(new MultiFab(ba, dm, 1, IntVect(1,1,0)));
    vec_pn_full_domain[0].reset(new MultiFab(ba, dm, 1, IntVect(1,1,0)));

    auto h_growvect = vec_h[0]->nGrowVect();
    auto pm_growvect = vec_pm[0]->nGrowVect();
    auto pn_growvect = vec_pn[0]->nGrowVect();
    for (int lev=1; lev <= hires_grid_level; lev++) {
        ba = ba.refine(refRatio(lev-1));
        refined_domain.refine(refRatio(lev-1));

        // Always allocate at least as many grow cells as there are in the level's normal variable multifab
        // This makes copying and boundary filling much easier
        vec_h_full_domain[lev].reset(new MultiFab(ba, dm, 1, max(cum_ref_ratios[lev],h_growvect)));
        vec_pm_full_domain[lev].reset(new MultiFab(ba, dm, 1, max(cum_ref_ratios[lev],pm_growvect)));
        vec_pn_full_domain[lev].reset(new MultiFab(ba, dm, 1, max(cum_ref_ratios[lev],pn_growvect)));
    }
    // A NetCDF read covers only as many grow cells as the file carries, and the analytic
    // bathymetry hook fills h alone, so parts of these arrays can reach the average-down
    // unwritten. Zero them so what lands at level 0 does not depend on the arena.
    for (int lev = 0; lev <= hires_grid_level; lev++) {
        vec_h_full_domain[lev]->setVal(0.0);
        vec_pm_full_domain[lev]->setVal(0.0);
        vec_pn_full_domain[lev]->setVal(0.0);
    }
    nc_hires_grid_box = refined_domain;
}

void
REMORA::init_bathymetry_full_domain_from_analytic ()
{
    // init_analytic_bathymetry needs to be able to handle the full number of grow cells that vec_h_full_domain has
    prob->init_analytic_bathymetry(hires_grid_level, Geom(hires_grid_level), solverChoice, *this, *vec_h_full_domain[hires_grid_level]);
    // Average down to fill levels below hires_grid_level. Use a special average_down so grow cells
    // get populated by averaged down fine data
    for (int lev=hires_grid_level-1; lev >= 0; lev--) {
        average_down_with_grow_cells(lev, vec_h_full_domain);
    }
}

void REMORA::allocate_init_full_domain () {
    // Make fake boxArray that covers the whole domain on level 0
    BoxArray ba;
    ba.define(geom[0].Domain());
    BoxArray ba2d;
    ba2d.define(makeSlab(geom[0].Domain(),2,0));
    Box refined_domain = geom[0].Domain();

    DistributionMapping dm(ba);
    vec_cons_full_domain[0].reset(new MultiFab(ba, dm, ncons, IntVect(1,1,0)));
    vec_xvel_full_domain[0].reset(new MultiFab(convert(ba,IntVect(1,0,0)), dm, 1, IntVect(0,1,0)));
    vec_yvel_full_domain[0].reset(new MultiFab(convert(ba,IntVect(0,1,0)), dm, 1, IntVect(1,0,0)));
    vec_zeta_full_domain[0].reset(new MultiFab(ba2d, dm, 1, IntVect(1,1,0)));


    auto cons_growvect = cons_new[0]->nGrowVect();
    auto xvel_growvect = xvel_new[0]->nGrowVect();
    auto yvel_growvect = yvel_new[0]->nGrowVect();
    auto zeta_growvect = vec_zeta[0]->nGrowVect();
    for (int lev=1; lev <= hires_init_level; lev++) {
        ba = ba.refine(refRatio(lev-1));
        ba2d = ba2d.refine(refRatio(lev-1));
        refined_domain.refine(refRatio(lev-1));

        // Always allocate at least as many grow cells as there are in the level's normal variable multifab
        // This makes copying and boundary filling much easier
        vec_cons_full_domain[lev].reset(new MultiFab(ba, dm, ncons, max(cum_ref_ratios[lev],cons_growvect)));
        vec_xvel_full_domain[lev].reset(new MultiFab(convert(ba,IntVect(1,0,0)), dm, 1, max(cum_ref_ratios[lev],xvel_growvect) - IntVect(1,0,0)));
        vec_yvel_full_domain[lev].reset(new MultiFab(convert(ba,IntVect(0,1,0)), dm, 1, max(cum_ref_ratios[lev],yvel_growvect) - IntVect(0,1,0)));
        vec_zeta_full_domain[lev].reset(new MultiFab(ba2d, dm, 1, max(cum_ref_ratios[lev],zeta_growvect)));
    }
    // See the note in allocate_bathymetry_grid_vars_full_domain: the vertical grow cells in
    // particular are never covered by the read, since cum_ref_ratios has no z component.
    for (int lev = 0; lev <= hires_init_level; lev++) {
        vec_cons_full_domain[lev]->setVal(0.0);
        vec_xvel_full_domain[lev]->setVal(0.0);
        vec_yvel_full_domain[lev]->setVal(0.0);
        vec_zeta_full_domain[lev]->setVal(0.0);
    }
    nc_hires_init_box = refined_domain;
}

void
REMORA::init_full_domain_from_analytic ()
{
    prob->init_analytic_prob(hires_init_level, geom[hires_init_level], solverChoice, *this, *vec_cons_full_domain[hires_init_level], *vec_xvel_full_domain[hires_init_level], *vec_yvel_full_domain[hires_init_level]);

    // Biology must be filled on the hires level before the average-down loop
    // below, or the biology components of vec_cons_full_domain are averaged
    // down uninitialized.
    init_biology_ic_full_domain();

    for (int lev=hires_init_level-1; lev >= 0; lev--) {
        average_down_with_grow_cells(lev, vec_cons_full_domain);
        average_down_with_grow_cells(lev, vec_xvel_full_domain);
        average_down_with_grow_cells(lev, vec_yvel_full_domain);
    }
}

void
REMORA::init_full_domain_zeta_from_analytic ()
{
    prob->init_analytic_zeta(hires_init_level, geom[hires_init_level], solverChoice, *this, *vec_zeta_full_domain[hires_init_level]);

    for (int lev=hires_init_level-1; lev >= 0; lev--) {
        average_down_with_grow_cells(lev, vec_zeta_full_domain);
    }
}
