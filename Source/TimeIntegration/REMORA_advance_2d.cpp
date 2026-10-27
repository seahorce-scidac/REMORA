#include <REMORA.H>

using namespace amrex;

namespace {
/** \brief Copy the freshly filled leapfrog component into the other two, in the coarse-fine
 *         ghost band only.
 *
 * FillPatch writes only component knew. The solver updates every component inside the fine
 * grid but nothing writes the coarse-fine ghost band, so its other components stay a full
 * parent step stale -- and the 2D solver reads zeta(krhs) and ubar(kstp) from exactly those
 * cells. ROMS's put_refine2d sets the contact points it is about to read, so make every
 * component agree with the parent state just filled.
 *
 * The band is found by marking the valid region and calling FillBoundary: ghosts another box
 * on this level covers pick up a 1, so the remaining zeros are the coarse-fine and
 * domain-boundary ghosts. Domain ghosts are excluded so the physical boundary conditions,
 * applied after the FillPatch, are not overwritten.
 */
void fill_ghost_kcomps (MultiFab& mf, int knew, const Geometry& geom)
{
    MultiFab valid(mf.boxArray(), mf.DistributionMap(), 1, mf.nGrowVect());
    valid.setVal(Real(0.0));
    valid.setVal(Real(1.0), 0, 1, 0);
    valid.FillBoundary(geom.periodicity());

    const Box dom = amrex::convert(geom.Domain(), mf.boxArray().ixType());

    const bool per_x = geom.isPeriodic(0);
    const bool per_y = geom.isPeriodic(1);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box gbx = mfi.growntilebox();
        const auto& a = mf.array(mfi);
        const auto& v = valid.const_array(mfi);
        const auto dlo = amrex::lbound(dom);
        const auto dhi = amrex::ubound(dom);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (v(i,j,k) > Real(0.5)) { return; }
            // Skip ghosts outside the domain only where the domain really ends. A ghost beyond
            // a periodic edge is not a physical-boundary cell: if the level wraps onto itself
            // the FillBoundary above already marked it valid, and if it does not -- a patch
            // covering part of the periodic width -- it is an ordinary coarse-fine ghost that
            // FillPatch just filled from the parent, needing the same synchronisation as any
            // other.
            if (!per_x && (i < dlo.x || i > dhi.x)) { return; }
            if (!per_y && (j < dlo.y || j > dhi.y)) { return; }
            const Real val = a(i,j,k,knew);
            a(i,j,k,0) = val;
            a(i,j,k,1) = val;
            a(i,j,k,2) = val;
        });
    }
}
} // namespace
/** Nonlinear shallow-water rpimitive equations predictor (Leap-frog) and
 * corrector (Adams-Moulton) time-stepping engine. Corresponds to Nonlinear/step2d_LF_AM3.h
 * in ROMS.
 *
 * @param[in   ] lev        level of refinement (coarsest level is 0)
 * @param[in   ] mf_rhoS    density perturbation
 * @param[in   ] mf_rhoA    vertically-averaged density
 * @param[inout] mf_ru2d    RHS contributions to 2D u-momentum
 * @param[inout] mf_rv2d    RHS contribtuions to 2D v-momentum
 * @param[inout] mf_rufrc   before first predictor, vertical integral of 3D RHS for uvel, converted to forcing terms
 * @param[inout] mf_rvfrc   before first predictor, vertical integral of 3D RHS for vvel, converted to forcing term
 * @param[inout] mf_Zt_avg1 average of sea surface height over all fast steps
 * @param[inout] mf_DU_avg1 time-averaged u-flux for 2D equations
 * @param[inout] mf_DU_avg2 time-averaged u-flux for 3D equation coupling
 * @param[inout] mf_DV_avg1 time-averaged v-flux for 2D equations
 * @param[inout] mf_DV_avg2 time-averaged v-flux for 3D equation coupling
 * @param[inout] mf_rubar   RHS of vertically integrated u-momentum
 * @param[inout] mf_rvbar   RHS of vertically integrated v-momentum
 * @param[inout] mf_rzeta   RHS of sea surface height
 * @param[inout] mf_ubar    vertically integrated u-momentum
 * @param[inout] mf_vbar    vertically integrated v-momentum
 * @param[inout] mf_zeta    Sea-surface height
 * @param[in   ] mf_h       Bathymetry
 * @param[in   ] mf_pm      1 / dx
 * @param[in   ] mf_pn      1 / dy
 * @param[in   ] mf_fcor    Coriolis factor
 * @param[inout] mf_visc2_p Harmonic viscosity at psi points
 * @param[inout] mf_visc2_r Harmoic viscosity at rho points
 * @param[in   ] mf_mskr    Land-sea mask at rho-points
 * @param[in   ] mf_msku    Land-sea mask at u-points
 * @param[in   ] mf_mskv    Land-sea mask at v-points
 * @param[in   ] mf_mskp    Land-sea mask at psi-points
 * @param[in   ] dtfast_lev Length of current barotropic step
 * @param[in   ] predictor_2d_step Is this a predictor step?
 * @param[in   ] first_2d_step Is this the first barotropic step?
 * @param[in   ] my_iif     Which barotropic predictor-corrector pair?
 * @param[inout] next_indx1 Cached index for
 */

void
REMORA::advance_2d (int lev,
                   MultiFab const* mf_rhoS,
                   MultiFab const* mf_rhoA,
                   MultiFab      * mf_ru2d,
                   MultiFab      * mf_rv2d,
                   MultiFab      * mf_rufrc,
                   MultiFab      * mf_rvfrc,
                   MultiFab      * mf_Zt_avg1,
                   std::unique_ptr<MultiFab>& mf_DU_avg1,
                   std::unique_ptr<MultiFab>& mf_DU_avg2,
                   std::unique_ptr<MultiFab>& mf_DV_avg1,
                   std::unique_ptr<MultiFab>& mf_DV_avg2,
                   std::unique_ptr<MultiFab>& mf_rubar,
                   std::unique_ptr<MultiFab>& mf_rvbar,
                   std::unique_ptr<MultiFab>& mf_rzeta,
                   std::unique_ptr<MultiFab>& mf_ubar,
                   std::unique_ptr<MultiFab>& mf_vbar,
                   MultiFab      * mf_zeta,
                   MultiFab const* mf_h,
                   MultiFab const* mf_pm,
                   MultiFab const* mf_pn,
                   MultiFab const* mf_fcor,
                   MultiFab const* mf_visc2_p,
                   MultiFab const* mf_visc2_r,
                   MultiFab const* mf_mskr,
                   MultiFab const* mf_msku,
                   MultiFab const* mf_mskv,
                   MultiFab const* mf_mskp,
                   Real dtfast_lev,
                   bool predictor_2d_step,
                   bool first_2d_step, int my_iif,
                   int & next_indx1)
{
    BL_PROFILE("REMORA::advance2d()");
    int iic = istep[lev];
    const int nstp  = 0;
    int ntfirst = 0;

    int knew = 3;
    int krhs = (my_iif + iic) % 2 + 1;
    int kstp = my_iif <=1 ? iic % 2 + 1 : (iic % 2 + my_iif % 2 + 1) % 2 + 1;
    int indx1 = krhs;
    if (predictor_2d_step) {
        next_indx1 = 3 - indx1;
    } else {
        knew = next_indx1;
        kstp = 3 - knew;
        krhs = 3;
        //If it's not the auxiliary time step, set indx1 to next_indx1
        // NOTE: should this ever not execute?
        // Include indx1 updates for diagnostic purposes?
        //        if (my_iif<nfast+1)
        //            indx1=next_indx1;
    }
    int ptsk = 3-kstp;
    knew-=1;
    krhs-=1;
    kstp-=1;
    // Include indx1 updates for diagnostic purposes?
    //indx1-=1;
    ptsk-=1;
    // Per-level scratch (see MakeNewLevel), every cell read below is written first in
    // this call, so nothing carries over between fast steps.
    MultiFab& mf_DUon = *vec_DUon[lev];
    MultiFab& mf_DVom = *vec_DVom[lev];

    // What one pass hands the next: the cell fields the face pass reads from both sides of
    // a face, and the momentum right-hand sides (face-centred, stored at their face index).
    // The rest of what the step used to stage here (Drhs, Dstp, fomn, rhs_zeta, zeta_new)
    // is recomputed where it is read, with the expression and operand order of the kernels
    // it replaces, so the results are unchanged.
    int ncomp = 0;
    int Dnew_comp = ncomp++;
    int zwrk_comp = ncomp++;
    int gzeta_comp = ncomp++;
    int gzeta2_comp = ncomp++;
    int gzetaSA_comp = ncomp++;
    int rhs_ubar_comp = ncomp++;
    int rhs_vbar_comp = ncomp++;

    AMREX_ALWAYS_ASSERT(ncomp == n_scratch_2d);
    MultiFab& mf = *vec_scratch_2d[lev];

    // Barotropic mass fluxes, both directions in one launch, with the depth
    // Drhs = zeta(krhs) + h on either side of the face
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*mf_rhoS, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real const> const& ubar = mf_ubar->const_array(mfi);
        Array4<Real const> const& vbar = mf_vbar->const_array(mfi);
        Array4<Real const> const& zeta = mf_zeta->const_array(mfi);
        Array4<Real const> const& h    = mf_h->const_array(mfi);

        Array4<Real const> const& on_u = vec_metric_u[lev]->const_array(mfi,MetricU::on_u);
        Array4<Real const> const& om_v = vec_metric_v[lev]->const_array(mfi,MetricV::om_v);

        auto DUon = mf_DUon.array(mfi);
        auto DVom = mf_DVom.array(mfi);

        ParallelFor(makeSlab(mfi.grownnodaltilebox(0, IntVect(NGROW,NGROW,0)),2,0),
                    makeSlab(mfi.grownnodaltilebox(1, IntVect(NGROW,NGROW,0)),2,0),
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Real cff1= Real(0.5) * on_u(i,j,0) *((zeta(i,j,0,krhs)+h(i,j,0))+
                                                (zeta(i-1,j,0,krhs)+h(i-1,j,0)));
            DUon(i,j,0)=ubar(i,j,0,krhs)*cff1;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Real cff1= Real(0.5) * om_v(i,j,0) * ((zeta(i,j,0,krhs)+h(i,j,0))+
                                                 (zeta(i,j-1,0,krhs)+h(i,j-1,0)));
            DVom(i,j,0)=vbar(i,j,0,krhs)*cff1;
        });
    }

    // These are needed to pass the tests with bathymetry but I don't quite see why
    if (do_substep && cf_impose_flux) {
        set_2d_cf_flux(lev, t_old[lev], mf_DUon, mf_DVom);
    }

    // One exchange for both so their messages are in flight together
    amrex::FillBoundary(Vector<MultiFab*>{&mf_DUon, &mf_DVom}, geom[lev].periodicity());

#ifdef REMORA_USE_NETCDF
    if (solverChoice.do_m2_clim_nudg) {
        ubar_clim_data_from_file->update_interpolated_to_time(model_time(t_new[lev]), lev, vec_ubar[lev].get(), geom, ref_ratio);
        vbar_clim_data_from_file->update_interpolated_to_time(model_time(t_new[lev]), lev, vec_vbar[lev].get(), geom, ref_ratio);
    }
#endif

    // The passes below each run over one box with a branch-free body: run-time choices
    // (time-stepping scheme, advection scheme, optional terms) are compile-time options,
    // so every loop vectorizes on CPU and each pass is one kernel on GPU. Passes are not
    // fused further where that slows the CPU: advection and viscosity fluxes, and the two
    // halves of the momentum right-hand side, ran 10-40% slower in one body than apart
    // (each loop then reads 20-35 arrays at once, and a body with per-cell region tests
    // does not vectorize at all).

    //
    //  Do not perform the actual time stepping during the auxiliary
    //  (nfast(ng)+1) time step: only the time averages are accumulated.
    //
    const bool do_step = (my_iif < nfast);

    //  Time-averaging weights
    auto weight1 = vec_weight1.dataPtr();
    auto weight2 = vec_weight2.dataPtr();
    enum AvgMode : int { avg_first_pred, avg_pred, avg_corr };
    const int avg_mode = predictor_2d_step ? (first_2d_step ? avg_first_pred : avg_pred)
                                           : avg_corr;
    Real cff1_wt1 = zero;
    Real cff2_wt  = zero;
    if (avg_mode == avg_first_pred) {
        cff2_wt = (Real(-1.0)/Real(12.0))*weight2[my_iif+1];
    } else if (avg_mode == avg_pred) {
        cff1_wt1 = weight1[my_iif-1];
        cff2_wt  = (Real(8.0)/Real(12.0))*weight2[my_iif]-
                   (one/Real(12.0))*weight2[my_iif+1];
    } else if (first_2d_step) {
        cff2_wt = weight2[my_iif];
    } else {
        cff2_wt = Real(5.0)/Real(12.0)*weight2[my_iif];
    }

    //
    //  During the first time-step, the predictor step is Forward-Euler and the corrector
    //  step is Backward-Euler. Otherwise, the predictor step is Leap-frog and the corrector
    //  step is Adams-Moulton.
    //
    enum Scheme : int { first_step, leapfrog, adams_moulton };
    const int scheme = (my_iif == 0) ? first_step
                                     : (predictor_2d_step ? leapfrog : adams_moulton);
    Real zcff1 = zero, zcff2 = zero, zcff3 = zero, zcff4 = zero, zcff5 = zero;
    Real mcff1 = zero, mcff2 = zero, mcff3 = zero;
    if (scheme == first_step) {
        zcff1 = dtfast_lev;
        mcff1 = Real(0.5)*dtfast_lev;
    } else if (scheme == leapfrog) {
        zcff1 = two * dtfast_lev;
        zcff4 = Real(4.0) / Real(25.0);
        zcff5 = one - two*zcff4;
        mcff1 = dtfast_lev;
    } else {
        zcff1 = dtfast_lev * Real(5.0)/Real(12.0);
        zcff2 = dtfast_lev * Real(8.0)/Real(12.0);
        zcff3 = dtfast_lev * one/Real(12.0);
        zcff4 = two/Real(5.0);
        zcff5 = one-zcff4;
        mcff1 = Real(0.5)*dtfast_lev*Real(5.0)/Real(12.0);
        mcff2 = Real(0.5)*dtfast_lev*Real(8.0)/Real(12.0);
        mcff3 = Real(0.5)*dtfast_lev*one/Real(12.0);
    }
    const Real fac = Real(1000.0)/solverChoice.rho0;
    const int  pred = predictor_2d_step ? 1 : 0;

    //  Coupling from 3D to 2D on the first predictor step: the first 3D step, the second,
    //  or AB3 after that. Every other fast step just adds the 3D forcing.
    enum Couple : int { couple_none, couple_first, couple_second, couple_ab3 };
    int couple = couple_none;
    if (first_2d_step && predictor_2d_step) {
        couple = (iic==ntfirst) ? couple_first
                                : ((iic==(ntfirst+1)) ? couple_second : couple_ab3);
    }
    const int add_frc = (couple == couple_none) ? 1 : 0;

    const Real g_half = Real(0.5) * solverChoice.g;
    const Real third  = one / Real(3.0);
    const Real sixth  = one / Real(6.0);

    const int upstream3 = (solverChoice.uv_Hadv_scheme == AdvectionScheme::upstream3) ? 1 : 0;
    const int use_cor   = solverChoice.use_coriolis ? 1 : 0;
    const int use_curv  = solverChoice.use_curvilinear_grid ? 1 : 0;

    // Boundary handling of the third-order advection
    const Box& domain = geom[lev].Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);
    const bool is_periodic_in_x = geom[0].isPeriodic(0);
    const bool is_periodic_in_y = geom[0].isPeriodic(1);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*mf_rhoS, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real const> const& rhoS = mf_rhoS->const_array(mfi);
        Array4<Real const> const& rhoA = mf_rhoA->const_array(mfi);
        Array4<Real const> const& h    = mf_h->const_array(mfi);

        Array4<Real      > const& rufrc   = mf_rufrc->array(mfi);
        Array4<Real      > const& rvfrc   = mf_rvfrc->array(mfi);
        Array4<Real      > const& Zt_avg1 = mf_Zt_avg1->array(mfi);
        Array4<Real      > const& ubar    = mf_ubar->array(mfi);
        Array4<Real      > const& vbar    = mf_vbar->array(mfi);
        Array4<Real      > const& zeta = mf_zeta->array(mfi);
        Array4<Real      > const& DU_avg1 = (mf_DU_avg1)->array(mfi);
        Array4<Real      > const& DU_avg2 = (mf_DU_avg2)->array(mfi);
        Array4<Real      > const& DV_avg1 = (mf_DV_avg1)->array(mfi);
        Array4<Real      > const& DV_avg2 = (mf_DV_avg2)->array(mfi);
        Array4<Real      > const& ru2d = (mf_ru2d)->array(mfi);
        Array4<Real      > const& rv2d = (mf_rv2d)->array(mfi);
        Array4<Real      > const& rubar = (mf_rubar)->array(mfi);
        Array4<Real      > const& rvbar = (mf_rvbar)->array(mfi);
        Array4<Real      > const& rzeta = (mf_rzeta)->array(mfi);
        Array4<Real const> const& visc2_p = mf_visc2_p->const_array(mfi);
        Array4<Real const> const& visc2_r = mf_visc2_r->const_array(mfi);

        Array4<Real const> const& pm   = mf_pm->const_array(mfi);
        Array4<Real const> const& pn   = mf_pn->const_array(mfi);
        Array4<Real const> const& fcor = mf_fcor->const_array(mfi);

        Array4<Real const> const& pmn    = vec_metric_r[lev]->const_array(mfi,MetricR::pmn);
        Array4<Real const> const& omn    = vec_metric_r[lev]->const_array(mfi,MetricR::omn);
        Array4<Real const> const& metr   = vec_metric_r[lev]->const_array(mfi);
        Array4<Real const> const& metp   = vec_metric_p[lev]->const_array(mfi);
        Array4<Real const> const& on_u   = vec_metric_u[lev]->const_array(mfi,MetricU::on_u);
        Array4<Real const> const& pmn4_u = vec_metric_u[lev]->const_array(mfi,MetricU::pmn4_u);
        Array4<Real const> const& om_v   = vec_metric_v[lev]->const_array(mfi,MetricV::om_v);
        Array4<Real const> const& pmn4_v = vec_metric_v[lev]->const_array(mfi,MetricV::pmn4_v);

        Array4<Real const> const& mskr = mf_mskr->const_array(mfi);
        Array4<Real const> const& msku = mf_msku->const_array(mfi);
        Array4<Real const> const& mskv = mf_mskv->const_array(mfi);
        Array4<Real const> const& mskp = mf_mskp->const_array(mfi);

        Array4<Real const> const& DUon = mf_DUon.const_array(mfi);
        Array4<Real const> const& DVom = mf_DVom.const_array(mfi);

        auto Dnew     = mf.array(mfi,Dnew_comp);
        auto zwrk     = mf.array(mfi,zwrk_comp);
        auto gzeta    = mf.array(mfi,gzeta_comp);
        auto gzeta2   = mf.array(mfi,gzeta2_comp);
        auto gzetaSA  = mf.array(mfi,gzetaSA_comp);
        auto rhs_ubar = mf.array(mfi,rhs_ubar_comp);
        auto rhs_vbar = mf.array(mfi,rhs_vbar_comp);

        Box bxD  = makeSlab(mfi.tilebox(),2,0);
        Box xbxD = makeSlab(mfi.nodaltilebox(0),2,0);
        Box ybxD = makeSlab(mfi.nodaltilebox(1),2,0);

        //
        // Time averages of the free surface and the barotropic fluxes
        //
        AnyCTO(TypeList<CompileTimeOptions<avg_first_pred,avg_pred,avg_corr>>{}, {avg_mode},
        [&] (auto const& fu, auto const& fv) {
            ParallelFor(makeSlab(mfi.grownnodaltilebox(0, IntVect(NGROW,NGROW,0)),2,0),
                        makeSlab(mfi.grownnodaltilebox(1, IntVect(NGROW,NGROW,0)),2,0), fu, fv);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int, auto mode)
        {
            if (mode.value == avg_first_pred) {
                DU_avg1(i,j,0)=zero;
                DU_avg2(i,j,0)=cff2_wt*DUon(i,j,0);
            } else if (mode.value == avg_pred) {
                DU_avg1(i,j,0) += cff1_wt1*DUon(i,j,0);
                DU_avg2(i,j,0) += cff2_wt*DUon(i,j,0);
            } else {
                DU_avg2(i,j,0)=DU_avg2(i,j,0)+cff2_wt*DUon(i,j,0);
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int, auto mode)
        {
            if (mode.value == avg_first_pred) {
                DV_avg1(i,j,0)=zero;
                DV_avg2(i,j,0)=cff2_wt*DVom(i,j,0);
            } else if (mode.value == avg_pred) {
                DV_avg1(i,j,0) += cff1_wt1*DVom(i,j,0);
                DV_avg2(i,j,0) += cff2_wt*DVom(i,j,0);
            } else {
                DV_avg2(i,j,0)=DV_avg2(i,j,0)+cff2_wt*DVom(i,j,0);
            }
        });

        if (predictor_2d_step) {
            Box gbx3D = makeSlab(mfi.growntilebox(IntVect(NGROW+1,NGROW+1,0)),2,0);
            if (first_2d_step) {
                ParallelFor(gbx3D, [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    Zt_avg1(i,j,0)=zero;
                });
            } else {
                ParallelFor(gbx3D, [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    Zt_avg1(i,j,0) += cff1_wt1*zeta(i,j,0,krhs);
                });
            }
        }

        if (!do_step) { continue; }

        //
        //  Free-surface step and the pressure-gradient fields, on the tile grown by one.
        //  zeta(knew) and rzeta(krhs) are stored on the whole grown tile rather than the
        //  growntilebox the solver needs: where two tiles of a box overlap, both compute the
        //  same value from the same exchanged fluxes, and nothing in this loop reads them.
        //
        Box tbxp1D = makeSlab(grow(mfi.tilebox(),IntVect(NGROW-1,NGROW-1,0)),2,0);
        AnyCTO(TypeList<CompileTimeOptions<first_step,leapfrog,adams_moulton>,
                        CompileTimeOptions<0,1>>{}, {scheme, pred},
        [&] (auto const& f) { ParallelFor(tbxp1D, f); },
        [=] AMREX_GPU_DEVICE (int i, int j, int, auto sch, auto is_pred)
        {
            Real rhs_zeta = zero;
            Real zeta_new, zw;
            if (sch.value == first_step) {
                rhs_zeta = (DUon(i,j,0)-DUon(i+1,j,0))+
                           (DVom(i,j,0)-DVom(i,j+1,0));
                zeta_new = (zeta(i,j,0,kstp)+ pmn(i,j,0)*zcff1*rhs_zeta) * mskr(i,j,0);
                zw = Real(0.5)*(zeta(i,j,0,kstp)+zeta_new);
            } else if (sch.value == leapfrog) {
                rhs_zeta = (DUon(i,j,0)-DUon(i+1,j,0))+
                           (DVom(i,j,0)-DVom(i,j+1,0));
                zeta_new = (zeta(i,j,0,kstp)+
                            pmn(i,j,0)*zcff1*rhs_zeta) * mskr(i,j,0);
                zw = zcff5*zeta(i,j,0,krhs)+
                     zcff4*(zeta(i,j,0,kstp)+zeta_new);
            } else {
                Real cff=zcff1*((DUon(i,j,0)-DUon(i+1,j,0))+
                                (DVom(i,j,0)-DVom(i,j+1,0)));
                zeta_new = zeta(i,j,0,kstp)+
                    pmn(i,j,0)*(cff+
                                zcff2*rzeta(i,j,0,kstp)-
                                zcff3*rzeta(i,j,0,ptsk));
                zeta_new *= mskr(i,j,0);
                zw = zcff5*zeta_new+zcff4*zeta(i,j,0,krhs);
            }
            Dnew(i,j,0) = zeta_new+h(i,j,0);
            zwrk(i,j,0) = zw;
            Real gz = (fac+rhoS(i,j,0))*zw;
            gzeta(i,j,0) = gz;
            gzeta2(i,j,0) = gz*zw;
            gzetaSA(i,j,0) = zw*(rhoS(i,j,0)-rhoA(i,j,0));

            //  Load new free-surface values into shared array at both predictor and
            //  corrector steps, and the predictor's right-hand side.
            zeta(i,j,0,knew) = zeta_new;
            if (is_pred.value) {
                rzeta(i,j,0,krhs) = rhs_zeta;
            }
        });

        //
        // Momentum fluxes at the rho and psi points around the tile: horizontal advection
        // and harmonic viscosity (uv3dmix with Hz = Drhs and no time step), and the
        // curvilinear terms. Per-tile, thread-private.
        //
        int nflx = 0;
        const int aUFx = nflx++, aUFe = nflx++, aVFx = nflx++, aVFe = nflx++; // advection
        const int hUFx = nflx++, hUFe = nflx++, hVFx = nflx++, hVFe = nflx++; // viscosity
        const int cUFx = nflx++, cVFe = nflx++;                               // curvilinear
        FArrayBox fab_flx(grow(bxD,IntVect(1,1,0)), nflx, The_Async_Arena());
        auto flx = fab_flx.array();

        // Cells: the x-fluxes of u on [lo-1:hi+1] x [lo:hi], the y-fluxes of v on
        // [lo:hi] x [lo-1:hi+1]
        AnyCTO(TypeList<CompileTimeOptions<0,1>>{}, {upstream3},
        [&] (auto const& fu, auto const& fv) {
            ParallelFor(growLo(xbxD,0,1), growLo(ybxD,1,1), fu, fv);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int, auto up3)
        {
            if (up3.value) {
                Real uxx_i   = ubar(i-1,j,0,krhs)-two*ubar(i  ,j,0,krhs)+ubar(i+1,j,0,krhs);
                Real uxx_ip1 = ubar(i  ,j,0,krhs)-two*ubar(i+1,j,0,krhs)+ubar(i+2,j,0,krhs);

                Real Huxx_i   = DUon(i-1,j,0)-two*DUon(i  ,j,0)+DUon(i+1,j,0);
                Real Huxx_ip1 = DUon(i  ,j,0)-two*DUon(i+1,j,0)+DUon(i+2,j,0);

                if (i == dlo.x && !is_periodic_in_x) {
                    uxx_i = uxx_ip1;
                    Huxx_i = Huxx_ip1;
                }
                else if (i == dhi.x && !is_periodic_in_x) {
                    uxx_ip1 = uxx_i;
                    Huxx_ip1 = Huxx_i;
                }

                Real ubar_avg = ubar(i  ,j,0,krhs)+ubar(i+1,j,0,krhs);

                flx(i,j,0,aUFx)=Real(0.25)*(ubar_avg-sixth*(uxx_i+uxx_ip1)) * (DUon(i,j,0)+ DUon(i+1,j,0)-sixth*(Huxx_i+ Huxx_ip1));
            } else {
                flx(i,j,0,aUFx) = Real(0.25) * (DUon(i,j,0)+DUon(i+1,j,0)) * (ubar(i,j,0,krhs)+ubar(i+1,j,0,krhs));
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int, auto up3)
        {
            if (up3.value) {
                Real vee_j    = vbar(i,j-1,0,krhs)-two*vbar(i,j  ,0,krhs)+vbar(i,j+1,0,krhs);
                Real vee_jp1  = vbar(i,j  ,0,krhs)-two*vbar(i,j+1,0,krhs)+vbar(i,j+2,0,krhs);

                Real Hvee_j   = DVom(i,j-1,0)-two*DVom(i,j  ,0)+DVom(i,j+1,0);
                Real Hvee_jp1 = DVom(i,j  ,0)-two*DVom(i,j+1,0)+DVom(i,j+2,0);

                if (j == dlo.y and !is_periodic_in_y) {
                    vee_j = vee_jp1;
                    Hvee_j = Hvee_jp1;
                }
                else if (j == dhi.y and !is_periodic_in_y) {
                    vee_jp1 = vee_j;
                    Hvee_jp1 = Hvee_j;
                }

                Real cff1=vbar(i,j  ,0,krhs)+vbar(i,j+1,0,krhs);

                flx(i,j,0,aVFe) = Real(0.25) * (cff1-(vee_j + vee_jp1)*sixth) * (DVom(i,j  ,0)+ DVom(i,j+1,0) -
                                                   sixth  * (Hvee_j+ Hvee_jp1));
            } else {
                flx(i,j,0,aVFe) = Real(0.25) * (DVom(i,j,0) + DVom(i,j+1,0))*(vbar(i,j,0,krhs)+vbar(i,j+1,0,krhs));
            }
        });

        ParallelFor(growLo(xbxD,0,1), growLo(ybxD,1,1),
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real cff = Real(0.5)*(zeta(i,j,0,krhs)+h(i,j,0)) * ( metr(i,j,0,MetricR::pmon_r) *
                                        ( (pn(i  ,j,0) + pn(i+1,j,0)) * ubar(i+1,j,0,krhs)-
                                          (pn(i-1,j,0) + pn(i  ,j,0)) * ubar(i  ,j,0,krhs) )-
                                                    metr(i,j,0,MetricR::pnom_r) *
                                        ( (pm(i,j  ,0) + pm(i,j+1,0)) * vbar(i,j+1,0,krhs)-
                                          (pm(i,j-1,0) + pm(i,j  ,0)) * vbar(i,j  ,0,krhs) ) );
            Real on_r = metr(i,j,0,MetricR::on_r);
            flx(i,j,0,hUFx) = on_r * on_r * visc2_r(i,j,0) * cff;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {

            const Real cff = Real(0.5)*(zeta(i,j,0,krhs)+h(i,j,0)) * ( metr(i,j,0,MetricR::pmon_r) *
                                        ( (pn(i  ,j,0) + pn(i+1,j,0)) * ubar(i+1,j,0,krhs)-
                                          (pn(i-1,j,0) + pn(i  ,j,0)) * ubar(i  ,j,0,krhs) )-
                                                    metr(i,j,0,MetricR::pnom_r) *
                                        ( (pm(i,j  ,0) + pm(i,j+1,0)) * vbar(i,j+1,0,krhs)-
                                          (pm(i,j-1,0) + pm(i,j  ,0)) * vbar(i,j  ,0,krhs) ) );
            Real om_r = metr(i,j,0,MetricR::om_r);
            flx(i,j,0,hVFe) = om_r * om_r * visc2_r(i,j,0) * cff;
        });

        // psi points: the y-fluxes of u and the x-fluxes of v, both on [lo:hi+1]^2
        AnyCTO(TypeList<CompileTimeOptions<0,1>>{}, {upstream3},
        [&] (auto const& fu, auto const& fv) {
            ParallelFor(growHi(xbxD,1,1), growHi(ybxD,0,1), fu, fv);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int, auto up3)
        {
            if (up3.value) {
                Real uee_j   = ubar(i,j-1,0,krhs)-two*ubar(i,j  ,0,krhs)+ubar(i,j+1,0,krhs);
                Real uee_jm1 = ubar(i,j-2,0,krhs)-two*ubar(i,j-1,0,krhs)+ubar(i,j  ,0,krhs);

                Real Hvxx_i   = DVom(i-1,j,0)-two*DVom(i  ,j,0)+DVom(i+1,j,0);
                Real Hvxx_im1 = DVom(i-2,j,0)-two*DVom(i-1,j,0)+DVom(i  ,j,0);

                if (j == dlo.y and !is_periodic_in_y) {
                    uee_jm1 = uee_j;
                } else if (j == dhi.y+1 and !is_periodic_in_y) {
                    uee_j = uee_jm1;
                }

                Real cff1=ubar(i,j  ,0,krhs)+ubar(i,j-1,0,krhs);
                Real cff2=DVom(i,j,0)+DVom(i-1,j,0);

                flx(i,j,0,aUFe)=Real(0.25)*(cff1-(uee_j+uee_jm1)*sixth)*
                  (cff2-sixth*(Hvxx_i+Hvxx_im1));
            } else {
                flx(i,j,0,aUFe) = Real(0.25) * (DVom(i,j,0)+DVom(i-1,j,0)) * (ubar(i,j,0,krhs)+ubar(i,j-1,0,krhs));
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int, auto up3)
        {
            if (up3.value) {
                Real vxx_i   = vbar(i-1,j,0,krhs)-two*vbar(i  ,j,0,krhs)+vbar(i+1,j,0,krhs);
                Real vxx_im1 = vbar(i-2,j,0,krhs)-two*vbar(i-1,j,0,krhs)+vbar(i  ,j,0,krhs);

                Real Huee_j   = DUon(i,j-1,0)-two*DUon(i,j  ,0)+DUon(i,j+1,0);
                Real Huee_jm1 = DUon(i,j-2,0)-two*DUon(i,j-1,0)+DUon(i,j  ,0);

                if (i == dlo.x and !is_periodic_in_x) {
                    vxx_im1 = vxx_i;
                } else if (i == dhi.x + 1 and !is_periodic_in_x) {
                    vxx_i = vxx_im1;
                }

                Real cff1=vbar(i  ,j,0,krhs)+vbar(i-1,j,0,krhs);
                Real cff2=DUon(i,j,0)+DUon(i,j-1,0);

                flx(i,j,0,aVFx)=Real(0.25)*(cff1-(vxx_i + vxx_im1)*sixth)* (cff2-sixth*(Huee_j+ Huee_jm1));
            } else {
                flx(i,j,0,aVFx) = Real(0.25) * (DUon(i,j,0) + DUon(i,j-1,0))*(vbar(i,j,0,krhs)+vbar(i-1,j,0,krhs));
            }
        });

        ParallelFor(growHi(xbxD,1,1), growHi(ybxD,0,1),
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real cff = mskp(i,j,0) * Real(0.125) *
                        ((zeta(i-1,j  ,0,krhs)+h(i-1,j  ,0)) + (zeta(i,j  ,0,krhs)+h(i,j  ,0))+
                         (zeta(i-1,j-1,0,krhs)+h(i-1,j-1,0)) + (zeta(i,j-1,0,krhs)+h(i,j-1,0)))*
                        (metp(i,j,0,MetricP::pmon_p)*
                         ((pn(i  ,j-1,0)+pn(i  ,j,0))*vbar(i  ,j,0,krhs)-
                          (pn(i-1,j-1,0)+pn(i-1,j,0))*vbar(i-1,j,0,krhs))+
                         metp(i,j,0,MetricP::pnom_p)*
                         ((pm(i-1,j  ,0)+pm(i,j  ,0))*ubar(i,j  ,0,krhs)-
                          (pm(i-1,j-1,0)+pm(i,j-1,0))*ubar(i,j-1,0,krhs)));
            const Real om_p = metp(i,j,0,MetricP::om_p);
            flx(i,j,0,hUFe) = om_p*om_p*visc2_p(i,j,0)*cff;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {

            const Real cff = mskp(i,j,0) * Real(0.125) *
                        ((zeta(i-1,j  ,0,krhs)+h(i-1,j  ,0)) + (zeta(i,j  ,0,krhs)+h(i,j  ,0))+
                         (zeta(i-1,j-1,0,krhs)+h(i-1,j-1,0)) + (zeta(i,j-1,0,krhs)+h(i,j-1,0)))*
                        (metp(i,j,0,MetricP::pmon_p)*
                         ((pn(i  ,j-1,0)+pn(i  ,j,0))*vbar(i  ,j,0,krhs)-
                          (pn(i-1,j-1,0)+pn(i-1,j,0))*vbar(i-1,j,0,krhs))+
                         metp(i,j,0,MetricP::pnom_p)*
                         ((pm(i-1,j  ,0)+pm(i,j  ,0))*ubar(i,j  ,0,krhs)-
                          (pm(i-1,j-1,0)+pm(i,j-1,0))*ubar(i,j-1,0,krhs)));
            const Real on_p = metp(i,j,0,MetricP::on_p);
            flx(i,j,0,hVFx) = on_p*on_p*visc2_p(i,j,0)*cff;
        });

        if (use_curv) {
            Array4<Real const> const& dndx = vec_dndx[lev]->const_array(mfi);
            Array4<Real const> const& dmde = vec_dmde[lev]->const_array(mfi);
            ParallelFor(grow(bxD,IntVect(1,1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                Real cff1 = Real(0.5) * (vbar(i,j,0,krhs) + vbar(i  ,j+1,0,krhs));
                Real cff2 = Real(0.5) * (ubar(i,j,0,krhs) + ubar(i+1,j  ,0,krhs));
                Real cff3 = cff1 * dndx(i,j,0);
                Real cff4 = cff2 * dmde(i,j,0);
                Real cff = (zeta(i,j,0,krhs)+h(i,j,0)) * (cff3 - cff4);
                flx(i,j,0,cUFx) = cff * cff1;
                flx(i,j,0,cVFe) = cff * cff2;
            });
        }

        //
        // Right-hand side of the 2D momentum equations, the terms in the order the
        // separate kernels added them: pressure gradient, advection, Coriolis (with
        // fomn = fcor/(pm*pn)), curvilinear, viscosity.
        //
        ParallelFor(xbxD, ybxD,
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Real rhs = g_half * on_u(i,j,0) *
                          ((    h(i-1,j,0) +     h(i,j,0))*
                           (gzeta(i-1,j,0) - gzeta(i,j,0))+
                           (    h(i-1,j,0) -     h(i,j,0))*
                           (      gzetaSA(i-1,j,0) + gzetaSA(i,j,0)+
                            third*(   rhoA(i-1,j,0) -    rhoA(i,j,0))*
                                 (   zwrk(i-1,j,0) -    zwrk(i,j,0)))+
                           (gzeta2(i-1,j,0)- gzeta2(i  ,j,0)));

            rhs -= ((flx(i,j  ,0,aUFx)-flx(i-1,j,0,aUFx)) +
                    (flx(i,j+1,0,aUFe)-flx(i  ,j,0,aUFe)));
            rhs_ubar(i,j,0) = rhs;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Real rhs = g_half*om_v(i,j,0) *
                          ((    h(i,j-1,0) +     h(i,j,0))*
                           (gzeta(i,j-1,0) - gzeta(i,j,0))+
                           (    h(i,j-1,0) -     h(i,j,0))*
                           (gzetaSA(i,j-1,0)+ gzetaSA(i,j  ,0)+
                            third*(rhoA(i,j-1,0)- rhoA(i,j  ,0))*
                                 (zwrk(i,j-1,0)- zwrk(i,j  ,0)))+
                           (gzeta2(i,j-1,0)- gzeta2(i,j  ,0)));

            rhs -= ((flx(i+1,j,0,aVFx)-flx(i  ,j,0,aVFx)) +
                    (flx(i  ,j,0,aVFe)-flx(i,j-1,0,aVFe)));
            rhs_vbar(i,j,0) = rhs;
        });

        AnyCTO(TypeList<CompileTimeOptions<0,1>,CompileTimeOptions<0,1>>{}, {use_cor, use_curv},
        [&] (auto const& fu, auto const& fv) { ParallelFor(xbxD, ybxD, fu, fv); },
        [=] AMREX_GPU_DEVICE (int i, int j, int, auto cor, auto curv)
        {
            Real rhs = rhs_ubar(i,j,0);
            if (cor.value) {
                Real UFx_i   = Real(0.5) * (zeta(i  ,j,0,krhs)+h(i  ,j,0)) * (fcor(i  ,j,0)*omn(i  ,j,0)) * (vbar(i  ,j,0,krhs)+vbar(i  ,j+1,0,krhs));
                Real UFx_im1 = Real(0.5) * (zeta(i-1,j,0,krhs)+h(i-1,j,0)) * (fcor(i-1,j,0)*omn(i-1,j,0)) * (vbar(i-1,j,0,krhs)+vbar(i-1,j+1,0,krhs));
                rhs += Real(0.5)*(UFx_i + UFx_im1);
            }

            if (curv.value) {
                rhs += Real(0.5) * (flx(i,j,0,cUFx) + flx(i-1,j,0,cUFx));
            }

            const Real cff1=Real(0.5)*(pn(i-1,j,0)+pn(i,j,0))*(flx(i,j  ,0,hUFx)-flx(i-1,j,0,hUFx));
            const Real cff2=Real(0.5)*(pm(i-1,j,0)+pm(i,j,0))*(flx(i,j+1,0,hUFe)-flx(i  ,j,0,hUFe));
            rhs_ubar(i,j,0) = rhs + (cff1+cff2);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int, auto cor, auto curv)
        {
            Real rhs = rhs_vbar(i,j,0);

            if (cor.value) {
                Real VFe_j   = Real(0.5) * (zeta(i,j  ,0,krhs)+h(i,j  ,0)) * (fcor(i,j  ,0)*omn(i,j  ,0)) * (ubar(i,j  ,0,krhs)+ubar(i+1,j  ,0,krhs));
                Real VFe_jm1 = Real(0.5) * (zeta(i,j-1,0,krhs)+h(i,j-1,0)) * (fcor(i,j-1,0)*omn(i,j-1,0)) * (ubar(i,j-1,0,krhs)+ubar(i+1,j-1,0,krhs));
                rhs -= Real(0.5)*(VFe_j + VFe_jm1);
            }

            if (curv.value) {
                rhs -= Real(0.5) * (flx(i,j,0,cVFe) + flx(i,j-1,0,cVFe));
            }

            const Real cff1=Real(0.5)*(pn(i,j-1,0)+pn(i,j,0))*(flx(i+1,j,0,hVFx)-flx(i,j  ,0,hVFx));
            const Real cff2=Real(0.5)*(pm(i,j-1,0)+pm(i,j,0))*(flx(i  ,j,0,hVFe)-flx(i,j-1,0,hVFe));
            rhs_vbar(i,j,0) = rhs + (cff1-cff2);
        });

#ifdef REMORA_USE_NETCDF
        if (solverChoice.do_m2_clim_nudg) {
            Array4<Real const> const& ubar_krhs = mf_ubar->const_array(mfi, krhs);
            Array4<Real const> const& vbar_krhs = mf_vbar->const_array(mfi, krhs);
            Array4<const Real> const& ubar_clim = ubar_clim_data_from_file->get_interpolated_mf(lev)->const_array(mfi);
            Array4<const Real> const& vbar_clim = vbar_clim_data_from_file->get_interpolated_mf(lev)->const_array(mfi);
            Array4<const Real> const& ubar_nudg_coeff = vec_nudg_coeff[bdy_ubar()][lev]->const_array(mfi);
            Array4<const Real> const& vbar_nudg_coeff = vec_nudg_coeff[bdy_vbar()][lev]->const_array(mfi);
            // Drhs, the depth the nudging is weighted with, at rho points around the tile
            FArrayBox fab_Drhs(grow(bxD,IntVect(1,1,0)),1,The_Async_Arena());
            auto Drhs = fab_Drhs.array();
            ParallelFor(grow(bxD,IntVect(1,1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                Drhs(i,j,0) = zeta(i,j,0,krhs)+h(i,j,0);
            });
            // Boxes are like this to match the ROMS loops i=IstrU..Iend, j=JstrV..Jend, which
            // exclude the u-face on the west/east domain edges and the v-face on south/north
            Box xbxD_adj = clim_nudg_momentum_box(mfi.nodaltilebox(0), 0,
                                                  Geom(lev).Domain(), Geom(lev).isPeriodic(0));
            xbxD_adj.makeSlab(2,0);
            Box ybxD_adj = clim_nudg_momentum_box(mfi.nodaltilebox(1), 1,
                                                  Geom(lev).Domain(), Geom(lev).isPeriodic(1));
            ybxD_adj.makeSlab(2,0);
            apply_clim_nudg(xbxD_adj, 1, 0, rhs_ubar, ubar_krhs, ubar_clim, ubar_nudg_coeff, Drhs,
                            vec_metric_u[lev]->const_array(mfi,MetricU::om_u), on_u);
            apply_clim_nudg(ybxD_adj, 0, 1, rhs_vbar, vbar_krhs, vbar_clim, vbar_nudg_coeff, Drhs,
                            om_v, vec_metric_v[lev]->const_array(mfi,MetricV::on_v));
        }
#endif

        //
        // Coupling from 3d to 2d on the first predictor step
        //
        if (couple != couple_none) {
            const Real ccff1 = Real(23.0)/Real(12.0);
            const Real ccff2 = Real(16.0)/Real(12.0);
            const Real ccff3 = Real(5.0)/Real(12.0);
            AnyCTO(TypeList<CompileTimeOptions<couple_first,couple_second,couple_ab3>>{}, {couple},
            [&] (auto const& fu, auto const& fv) { ParallelFor(xbxD, ybxD, fu, fv); },
            [=] AMREX_GPU_DEVICE (int i, int j, int, auto cpl)
            {
                if (cpl.value == couple_first) {
                    rufrc(i,j,0)    -= rhs_ubar(i,j,0);
                    rhs_ubar(i,j,0) += rufrc(i,j,0);
                    ru2d(i,j,0,nstp)  = rufrc(i,j,0);
                } else if (cpl.value == couple_second) {
                    rufrc(i,j,0)=rufrc(i,j,0)-rhs_ubar(i,j,0);
                    rhs_ubar(i,j,0)=rhs_ubar(i,j,0)+Real(1.5)*rufrc(i,j,0)-Real(0.5)*ru2d(i,j,0,0);
                    ru2d(i,j,0,1)=rufrc(i,j,0);
                    Real r_swap= ru2d(i,j,0,1);
                    ru2d(i,j,0,1) = ru2d(i,j,0,0);
                    ru2d(i,j,0,0) = r_swap;
                } else {
                    rufrc(i,j,0)=rufrc(i,j,0)-rhs_ubar(i,j,0);
                    rhs_ubar(i,j,0)=rhs_ubar(i,j,0)+
                        ccff1*rufrc(i,j,0)-
                        ccff2*ru2d(i,j,0,0)+
                        ccff3*ru2d(i,j,0,1);
                    ru2d(i,j,0,1)=rufrc(i,j,0);
                    Real r_swap= ru2d(i,j,0,1);
                    ru2d(i,j,0,1) = ru2d(i,j,0,0);
                    ru2d(i,j,0,0) = r_swap;
                }
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int, auto cpl)
            {
                if (cpl.value == couple_first) {
                    rvfrc(i,j,0)    -= rhs_vbar(i,j,0);
                    rhs_vbar(i,j,0) += rvfrc(i,j,0);
                    rv2d(i,j,0,nstp)  = rvfrc(i,j,0);
                } else if (cpl.value == couple_second) {
                    rvfrc(i,j,0)=rvfrc(i,j,0)-rhs_vbar(i,j,0);
                    rhs_vbar(i,j,0)=rhs_vbar(i,j,0)+Real(1.5)*rvfrc(i,j,0)-Real(0.5)*rv2d(i,j,0,0);
                    rv2d(i,j,0,1)=rvfrc(i,j,0);
                    Real r_swap= rv2d(i,j,0,1);
                    rv2d(i,j,0,1) = rv2d(i,j,0,0);
                    rv2d(i,j,0,0) = r_swap;
                } else {
                    rvfrc(i,j,0)=rvfrc(i,j,0)-rhs_vbar(i,j,0);
                    rhs_vbar(i,j,0)=rhs_vbar(i,j,0)+
                          ccff1*rvfrc(i,j,0)-
                          ccff2*rv2d(i,j,0,0)+
                          ccff3*rv2d(i,j,0,1);
                    rv2d(i,j,0,1)=rvfrc(i,j,0);
                    Real r_swap= rv2d(i,j,0,1);
                    rv2d(i,j,0,1) = rv2d(i,j,0,0);
                    rv2d(i,j,0,0) = r_swap;
                }
            });
        }

        //
        // Time step of the 2D momentum equations, with the 3D forcing added first on every
        // fast step but the coupling one, and Dstp = zeta(kstp) + h. If predictor step,
        // load right-side-term into shared arrays for the subsequent corrector step.
        //
        AnyCTO(TypeList<CompileTimeOptions<0,1>,CompileTimeOptions<0,1>,CompileTimeOptions<0,1>>{},
               {scheme == adams_moulton ? 1 : 0, pred, add_frc},
        [&] (auto const& fu, auto const& fv) { ParallelFor(xbxD, ybxD, fu, fv); },
        [=] AMREX_GPU_DEVICE (int i, int j, int, auto am3, auto is_pred, auto frc)
        {
            Real rhs = rhs_ubar(i,j,0);
            if (frc.value) { rhs += rufrc(i,j,0); }
            Real cff=pmn4_u(i,j,0);
            Real Dnew_avg=one/(Dnew(i,j,0)+Dnew(i-1,j,0));
            Real Dstp_sum=(zeta(i,j,0,kstp)+h(i,j,0))+(zeta(i-1,j,0,kstp)+h(i-1,j,0));
            if (am3.value) {
                ubar(i,j,0,knew)=(ubar(i,j,0,kstp)*Dstp_sum+
                                 cff*(mcff1*rhs+
                                      mcff2*rubar(i,j,0,kstp)-
                                      mcff3*rubar(i,j,0,ptsk)))*Dnew_avg * msku(i,j,0);
            } else {
                ubar(i,j,0,knew)=(ubar(i,j,0,kstp)*Dstp_sum+
                                  cff*mcff1*rhs)*Dnew_avg * msku(i,j,0);
            }
            if (is_pred.value) {
                rubar(i,j,0,krhs)=rhs;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int, auto am3, auto is_pred, auto frc)
        {
            Real rhs = rhs_vbar(i,j,0);
            if (frc.value) { rhs += rvfrc(i,j,0); }
            Real cff=pmn4_v(i,j,0);
            Real Dnew_avg=one/(Dnew(i,j,0)+Dnew(i,j-1,0));
            Real Dstp_sum=(zeta(i,j,0,kstp)+h(i,j,0))+(zeta(i,j-1,0,kstp)+h(i,j-1,0));
            if (am3.value) {
                vbar(i,j,0,knew)=(vbar(i,j,0,kstp)*Dstp_sum+
                                 cff*(mcff1*rhs+
                                      mcff2*rvbar(i,j,0,kstp)-
                                      mcff3*rvbar(i,j,0,ptsk)))*Dnew_avg * mskv(i,j,0);
            } else {
                vbar(i,j,0,knew)=(vbar(i,j,0,kstp)*Dstp_sum+
                                  cff*mcff1*rhs)*Dnew_avg * mskv(i,j,0);
            }
            if (is_pred.value) {
                rvbar(i,j,0,krhs)=rhs;
            }
        });
    }

    // Don't do the FillPatch or rivers at the last truncated predictor step.
    // We may need to move the zeta FillPatch further up
    if (my_iif<nfast) {
        int know;
        Real dt2d;
        if (my_iif==0) {
            know = krhs;
            dt2d = dtfast_lev;
        } else if (predictor_2d_step) {
            know = krhs;
            dt2d = two * dtfast_lev;
        } else {
            know = kstp;
            dt2d = dtfast_lev;
        }

        MultiFab ubar_know(*vec_ubar[lev], make_alias, know, 1);
        MultiFab vbar_know(*vec_vbar[lev], make_alias, know, 1);
        MultiFab zeta_know(*vec_zeta[lev], make_alias, know, 1);

        if (lev == 0) {
            // At level 0 FillPatch is one exchange of component knew followed by the
            // boundary conditions, and those read only `know` components (the field's own,
            // and zeta's for Flather and Chapman), h and the masks -- nothing the exchange
            // of knew changes. Post the three exchanges together so their messages are in
            // flight at once -- this runs every fast step -- then apply the conditions in
            // the same order as before.
            const Periodicity& per = geom[lev].periodicity();
            amrex::FillBoundary(Vector<MultiFab*>{vec_ubar[lev].get(), vec_vbar[lev].get(),
                                                  vec_zeta[lev].get()},
                                Vector<int>{knew, knew, knew}, Vector<int>{1, 1, 1},
                                Vector<IntVect>{vec_ubar[lev]->nGrowVect(),
                                                vec_vbar[lev]->nGrowVect(),
                                                vec_zeta[lev]->nGrowVect()},
                                Vector<Periodicity>{per, per, per});
            FillPatchBCs(lev, t_old[lev], *vec_ubar[lev], ubar_bc(), bdy_ubar(),
                         knew, 1, 0, know, dt2d, ubar_know);
            FillPatchBCs(lev, t_old[lev], *vec_vbar[lev], vbar_bc(), bdy_vbar(),
                         knew, 1, 0, know, dt2d, vbar_know);
            FillPatchBCs(lev, t_old[lev], *vec_zeta[lev], zeta_bc(), bdy_zeta(),
                         knew, 1, 0, know, dt2d, zeta_know);
        } else {
            FillPatch(lev, t_old[lev], *vec_ubar[lev], GetVecOfPtrs(vec_ubar), ubar_bc(), bdy_ubar(),
                      knew, false,true, 0,know, dt2d, ubar_know);
            FillPatch(lev, t_old[lev], *vec_vbar[lev], GetVecOfPtrs(vec_vbar), vbar_bc(), bdy_vbar(),
                      knew, false,true, 0,know, dt2d, vbar_know);
            // With the snapshot pair the coarse contribution is interpolated onto t_old[lev],
            // which timeStep guarantees lies inside the parent's step; without it FillPatch
            // falls back to the parent frozen at the end of that step. See roll_2d_snapshot.
            if (cf_time_interp_zeta && int(vec_zeta_crse_old.size()) >= lev
                && vec_zeta_crse_old[lev-1] && vec_zeta_crse_new[lev-1]) {
                FillPatch(lev, t_old[lev], *vec_zeta[lev], GetVecOfPtrs(vec_zeta), zeta_bc(), bdy_zeta(),
                          knew, false,false, 0,know, dt2d, zeta_know,
                          GetVecOfPtrs(vec_zeta_crse_old), GetVecOfPtrs(vec_zeta_crse_new));
            } else {
                FillPatch(lev, t_old[lev], *vec_zeta[lev], GetVecOfPtrs(vec_zeta), zeta_bc(), bdy_zeta(),
                          knew, false,false, 0,know, dt2d, zeta_know);
            }
        }

        // See fill_ghost_kcomps: without this the components the solver reads in the ghost
        // band lag the one FillPatch just wrote by a whole parent step.
        if (lev > 0 && cf_fill_all_kcomp) {
            fill_ghost_kcomps(*vec_zeta[lev], knew, geom[lev]);
            fill_ghost_kcomps(*vec_ubar[lev], knew, geom[lev]);
            fill_ghost_kcomps(*vec_vbar[lev], knew, geom[lev]);
        }

        // Replace the interface faces the FillPatchers just set from the parent's ubar with
        // the parent's mass flux, which conserves mass. Must follow the FillPatch.
        //
        // cf_set_2d_bcs = 2 imposes it only on the first fast step, where ROMS calls
        // put_refine2d (main3d.F, before the barotropic loop; inside it only composite grids
        // and the NESTING_DEBUG check run). Every fast step instead re-applies the parent's
        // step-averaged transport twenty times per step rather than letting the barotropic
        // solver carry the interface.
        const bool set_2d_now = (cf_set_2d_bcs == 1) || (cf_set_2d_bcs == 2 && my_iif == 0);
        if (do_substep && set_2d_now) {
            set_2d_cf_bcs(lev, t_old[lev], know, knew);
            // set_2d_cf_bcs writes each box's own interface faces. Where a box boundary meets
            // the interface, the neighbouring box's ghost copy of that face still holds the
            // FillPatch value, and the shared face between the two boxes would be advanced
            // from different neighbours -- the answer then depends on the box layout.
            const Periodicity& per = geom[lev].periodicity();
            amrex::FillBoundary(Vector<MultiFab*>{vec_ubar[lev].get(), vec_vbar[lev].get()},
                                Vector<int>{knew, knew}, Vector<int>{1, 1},
                                Vector<IntVect>{vec_ubar[lev]->nGrowVect(),
                                                vec_vbar[lev]->nGrowVect()},
                                Vector<Periodicity>{per, per});
        }

#ifdef REMORA_USE_NETCDF
        if (solverChoice.do_rivers) {
            river_source_transportbar->update_interpolated_to_time(model_time(t_old[lev]));
            int* river_direction_d = river_direction.data();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(*mf_rhoS, TilingIfNotGPU()); mfi.isValid(); ++mfi )
            {
                Array4<const int > const& river_pos = vec_river_position[lev]->const_array(mfi);
                Array4<const Real> const& river_transportbar = river_source_transportbar->fab_interp->array();
                Array4<Real      > const& ubar    = mf_ubar->array(mfi);
                Array4<Real      > const& vbar    = mf_vbar->array(mfi);
                Array4<Real const> const& zeta    = mf_zeta->const_array(mfi);
                Array4<Real const> const& h       = mf_h->const_array(mfi);
                Array4<Real const> const& on_u    = vec_metric_u[lev]->const_array(mfi,MetricU::on_u);
                Array4<Real const> const& om_v    = vec_metric_v[lev]->const_array(mfi,MetricV::om_v);

                Box gbx1D = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
                gbx1D.makeSlab(2,0);

                ParallelFor(gbx1D, [=] AMREX_GPU_DEVICE (int i, int j, int )
                {
                    int iriver = river_pos(i,j,0);
                    if (iriver >= 0) {
                        if (river_direction_d[iriver] == 0) {
                            Real cff = one / (on_u(i,j,0) * Real(0.5) * (zeta(i-1,j,0,knew) + h(i-1,j,0) +
                                        zeta(i,j,0,knew) + h(i,j,0)));
                            ubar(i,j,0,knew) = river_transportbar(iriver,0,0) * cff;
                        } else {
                            Real cff = one / (om_v(i,j,0) * Real(0.5) * (zeta(i,j-1,0,knew) + h(i,j-1,0) +
                                        zeta(i,j,0,knew) + h(i,j,0)));
                            vbar(i,j,0,knew) = river_transportbar(iriver,0,0) * cff;
                        }
                    }
                });
            }
            if (lev == 0) {
                // FillPatchNoBC at level 0 is just the exchange: do both in one
                const Periodicity& per = geom[lev].periodicity();
                amrex::FillBoundary(Vector<MultiFab*>{vec_ubar[lev].get(), vec_vbar[lev].get()},
                                    Vector<int>{knew, knew}, Vector<int>{1, 1},
                                    Vector<IntVect>{vec_ubar[lev]->nGrowVect(),
                                                    vec_vbar[lev]->nGrowVect()},
                                    Vector<Periodicity>{per, per});
            } else {
                FillPatchNoBC(lev, t_old[lev], *vec_ubar[lev], GetVecOfPtrs(vec_ubar), bdy_ubar(),
                          knew, false,false);
                FillPatchNoBC(lev, t_old[lev], *vec_vbar[lev], GetVecOfPtrs(vec_vbar), bdy_vbar(),
                          knew, false,false);
            }
        }
#endif
    }
}
