#include <iomanip>
#include <REMORA.H>
#include <REMORA_PhysBCFunct.H>

using namespace amrex;

namespace {
    PhysBCFunctNoOp no_bc;

/** \brief Replace each group of nr faces along the interface by their mean.
 *
 * ROMS hands every fine face under one parent face the same parent flux: get_persisted2d
 * copies DU_avg2 at the donor face picked by integer division of the fine index, and
 * put_refine2d and u2dbc_im then scale it by the edge-length ratio. AMReX's face
 * interpolator instead gives the flux a linear variation along the interface. Both make the
 * fine fluxes sum to the parent's; they differ in how the total is shared out. Because that
 * interpolation is conservative, the group mean is the parent value, so averaging turns the
 * linear profile back into ROMS's constant one.
 *
 * Only groups lying wholly inside a box are touched, and only faces the set mask covers, so
 * faces the parent never wrote keep their zero.
 */
void group_average_faces (MultiFab& mf, const iMultiFab& mask, int set_mask, int tdir, int nr)
{
    if (nr <= 1) { return; }
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        const auto& a = mf.array(mfi);
        const auto& m = mask.const_array(mfi);
        const int lo_t = bx.smallEnd(tdir);
        const int hi_t = bx.bigEnd(tdir);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            const int t  = (tdir == 0) ? i : j;
            const int t0 = amrex::coarsen(t, nr) * nr;
            // one thread per group, and only where the whole group is in this box
            if (t != t0 || t0 < lo_t || t0 + nr - 1 > hi_t) { return; }
            Real sum = Real(0.0);
            int  cnt = 0;
            for (int s = t0; s < t0 + nr; ++s) {
                const int si = (tdir == 0) ? s : i;
                const int sj = (tdir == 0) ? j : s;
                if (m(si,sj,k) != set_mask) { continue; }
                sum += a(si,sj,k); ++cnt;
            }
            if (cnt == 0) { return; }
            const Real avg = sum / Real(cnt);
            for (int s = t0; s < t0 + nr; ++s) {
                const int si = (tdir == 0) ? s : i;
                const int sj = (tdir == 0) ? j : s;
                if (m(si,sj,k) != set_mask) { continue; }
                a(si,sj,k) = avg;
            }
        });
    }
}
}

/**
 * Store this level's barotropic mass flux per unit cell edge length, for a finer level to
 * interpolate.
 *
 * @param[in] lev            level of refinement
 */
void
REMORA::store_2d_flux (int lev)
{
    // The previous step's end-of-step value is this step's start-of-step value: DU_avg2 is
    // accumulated during a step and only meaningful once it is over.
    MultiFab::Copy(*vec_Dubar_old[lev], *vec_Dubar_new[lev], 0, 0, 1,
                   vec_Dubar_new[lev]->nGrowVect());
    MultiFab::Copy(*vec_Dvbar_old[lev], *vec_Dvbar_new[lev], 0, 0, 1,
                   vec_Dvbar_new[lev]->nGrowVect());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*vec_Dubar_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real const> const& pn      = vec_pn[lev]->const_array(mfi);
        Array4<Real const> const& DU_avg2 = vec_DU_avg2[lev]->const_array(mfi);
        Array4<Real      > const& Dubar   = vec_Dubar_new[lev]->array(mfi);

        Box ubx = mfi.grownnodaltilebox(0,IntVect(NGROW,NGROW,0));
        ubx.makeSlab(2,0);

        // DU_avg2 is the flux through a whole face, m^3/s. Dividing by on_u leaves D*ubar,
        // so a finer face can multiply its own edge length back in and the fine fluxes sum
        // to the coarse one. Passing ubar across, or the raw flux, loses that.
        ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Real on_u = two / (pn(i,j,0) + pn(i-1,j,0));
            Dubar(i,j,0) = DU_avg2(i,j,0) / on_u;
        });
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*vec_Dvbar_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real const> const& pm      = vec_pm[lev]->const_array(mfi);
        Array4<Real const> const& DV_avg2 = vec_DV_avg2[lev]->const_array(mfi);
        Array4<Real      > const& Dvbar   = vec_Dvbar_new[lev]->array(mfi);

        Box vbx = mfi.grownnodaltilebox(1,IntVect(NGROW,NGROW,0));
        vbx.makeSlab(2,0);

        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            Real om_v = two / (pm(i,j,0) + pm(i,j-1,0));
            Dvbar(i,j,0) = DV_avg2(i,j,0) / om_v;
        });
    }

    // There is no previous step on the first one, so extrapolate rather than leave old at
    // the zero it was initialised with. Only remora.time_interp_flux reads old at all.
    if (istep[lev] == 0) {
        MultiFab::Copy(*vec_Dubar_old[lev], *vec_Dubar_new[lev], 0, 0, 1,
                       vec_Dubar_new[lev]->nGrowVect());
        MultiFab::Copy(*vec_Dvbar_old[lev], *vec_Dvbar_new[lev], 0, 0, 1,
                       vec_Dvbar_new[lev]->nGrowVect());
    }
}

/**
 * Set the normal barotropic velocity on this level's coarse-fine interface from the parent's
 * time-averaged mass flux, as ROMS does in put_refine2d (nesting.F):
 *
 *     ubar_f = Dubar_c / D_f,    D_f = 0.5*(h + zeta)_{i-1} + 0.5*(h + zeta)_i
 *
 * D comes from know by default, but ROMS uses one time index for both sides: the nested
 * branch of u2dbc_im.F builds D from zeta(kout) and writes ubar(kout), and put_refine2d uses
 * indx1 for both. Since the next half-step forms DUon = ubar(krhs)*D(krhs) with krhs equal to
 * this knew, ROMS recovers Dubar_parent exactly while this carries an extra D(knew)/D(know).
 * remora.cf_d_knew = 1 matches ROMS. Measured negligible -- at most 1.2e-07 on the step-1
 * Channel_Test interface velocity against a 2.4e-03 artifact, and no change to Dogbone volume
 * drift or to the Channel_Test blow-up -- so it is off by default pending gold regeneration.
 *
 * Momentum only. setup_step resets all three zeta components to Zt_avg1, so what a finer
 * level interpolates for the free surface is already the parent's fast-time average -- as in
 * ROMS, where set_zeta runs ahead of put_refine2d.
 *
 * @param[in] lev            level of refinement
 * @param[in] time           simulation time to interpolate the parent's flux to
 * @param[in] know           zeta time component to take the sea surface height from
 * @param[in] knew           ubar/vbar time component to set
 */
void
REMORA::set_2d_cf_bcs (int lev, Real time, int know, int knew)
{
    if (lev == 0 || cf_set_width < 0) { return; }

    BL_PROFILE("REMORA::set_2d_cf_bcs()");

    // Overwrites what the fill patchers just set from the parent's ubar, which would not
    // conserve mass. Every fast step: unlike a ROMS contact point on a physical perimeter,
    // these faces are interior and the barotropic solver rewrites them.
    const int set_mask = FPr_Dubar[lev-1].GetSetMaskVal();

    // ROMS builds D from the same zeta component it writes the velocity to, so the next
    // half-step's DUon = ubar(krhs)*D(krhs) reproduces the parent's transport exactly.
    const int kD = cf_d_knew ? knew : know;

    MultiFab Dubar_cf(vec_Dubar_new[lev]->boxArray(), vec_Dubar_new[lev]->DistributionMap(),
                      1, vec_Dubar_new[lev]->nGrowVect());
    MultiFab Dvbar_cf(vec_Dvbar_new[lev]->boxArray(), vec_Dvbar_new[lev]->DistributionMap(),
                      1, vec_Dvbar_new[lev]->nGrowVect());
    Dubar_cf.setVal(zero);
    Dvbar_cf.setVal(zero);

    FPr_Dubar[lev-1].FillSet(Dubar_cf, time, no_bc, domain_bcs_type);
    FPr_Dvbar[lev-1].FillSet(Dvbar_cf, time, no_bc, domain_bcs_type);

    // an x-face varies along y across the interface, and a y-face along x
    if (cf_flux_pc) {
        const IntVect rr = refRatio(lev-1);
        group_average_faces(Dubar_cf, *FPr_Dubar[lev-1].GetMask(), set_mask, 1, rr[1]);
        group_average_faces(Dvbar_cf, *FPr_Dvbar[lev-1].GetMask(), set_mask, 0, rr[0]);
    }


#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*vec_ubar[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real      > const& ubar  = vec_ubar[lev]->array(mfi);
        Array4<Real const> const& Dubar = Dubar_cf.const_array(mfi);
        Array4<int  const> const& cmask = FPr_Dubar[lev-1].GetMask()->const_array(mfi);
        Array4<Real const> const& zeta  = vec_zeta[lev]->const_array(mfi);
        Array4<Real const> const& h     = vec_h[lev]->const_array(mfi);
        Array4<Real const> const& msku  = vec_msku[lev]->const_array(mfi);

        // The mask carries no ghost cells, so this is exactly the region it covers.
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            if (cmask(i,j,0) != set_mask) { return; }

            Real D = Real(0.5) * (h(i-1,j,0,0) + zeta(i-1,j,0,kD) +
                                  h(i  ,j,0,0) + zeta(i  ,j,0,kD));
            if (D <= zero) { return; }

            ubar(i,j,0,knew) = Dubar(i,j,0) / D * msku(i,j,0);
        });
    }

    // Diagnostic: the fine level is not written to NetCDF, so dump the imposed u-face values
    // here. The last block printed for a step is the one the step ends on.
    // one block per baroclinic step: only the first fast step, from cf_print_iface onward.
    // cf_print_iface = N lines up with ROMS child history record N.
    if (cf_print_iface >= 0 && istep[0] >= cf_print_iface) {
        const auto dx_p = geom[lev].CellSizeArray();
        const auto lo_p = geom[lev].ProbLoArray();
        for (MFIter mfi(*vec_ubar[lev]); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto const& ubar  = vec_ubar[lev]->const_array(mfi);
            auto const& cmask = FPr_Dubar[lev-1].GetMask()->const_array(mfi);
            const auto lo = lbound(bx); const auto hi = ubound(bx);
            for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (cmask(i,j,0) != set_mask) { continue; }
                    amrex::AllPrint() << "[IFACE] lev " << lev
                        << " t " << std::setprecision(12) << time
                        << " told0 " << t_old[0] << " tnew0 " << t_new[0]
                        << " i " << i << " j " << j
                        << " x " << lo_p[0] + i * dx_p[0]
                        << " y " << lo_p[1] + (j + Real(0.5)) * dx_p[1]
                        << " ubar " << std::setprecision(12) << ubar(i,j,0,knew) << "\n";
                }
            }
        }

        // The neighbourhood of the western interface at mid-channel: what the fine level
        // actually holds in and around its ghost region, which no output file carries.
        for (MFIter mfi(*vec_zeta[lev]); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto lo = lbound(bx); const auto hi = ubound(bx);
            auto const& zeta  = vec_zeta[lev]->const_array(mfi);
            auto const& ubar  = vec_ubar[lev]->const_array(mfi);
            auto const& cmask = FPr_Dubar[lev-1].GetMask()->const_array(mfi);
            const int jmid = (lo.y + hi.y) / 2;
            int iface = -1;
            for (int i = lo.x; i <= hi.x; ++i) {
                if (cmask(i,jmid,0) == set_mask) { iface = i; break; }
            }
            if (iface < 0) { continue; }
            for (int off = -3; off <= 3; ++off) {
                const int i = iface + off;
                amrex::AllPrint() << "[HALO] lev " << lev
                    << " t " << std::setprecision(12) << time
                    << " j " << jmid << " off " << off << " i " << i
                    << " x " << lo_p[0] + i * dx_p[0]
                    << " zeta_know " << zeta(i,jmid,0,know)
                    << " zeta_knew " << zeta(i,jmid,0,knew)
                    << " ubar_know " << ubar(i,jmid,0,know)
                    << " ubar_knew " << ubar(i,jmid,0,knew)
                    << " mask " << cmask(i,jmid,0) << "\n";
            }
            break;   // one box is enough for this diagnostic
        }
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*vec_vbar[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real      > const& vbar  = vec_vbar[lev]->array(mfi);
        Array4<Real const> const& Dvbar = Dvbar_cf.const_array(mfi);
        Array4<int  const> const& cmask = FPr_Dvbar[lev-1].GetMask()->const_array(mfi);
        Array4<Real const> const& zeta  = vec_zeta[lev]->const_array(mfi);
        Array4<Real const> const& h     = vec_h[lev]->const_array(mfi);
        Array4<Real const> const& mskv  = vec_mskv[lev]->const_array(mfi);

        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            if (cmask(i,j,0) != set_mask) { return; }

            Real D = Real(0.5) * (h(i,j-1,0,0) + zeta(i,j-1,0,kD) +
                                  h(i,j  ,0,0) + zeta(i,j  ,0,kD));
            if (D <= zero) { return; }

            vbar(i,j,0,knew) = Dvbar(i,j,0) / D * mskv(i,j,0);
        });
    }
}

/**
 * Impose the parent's barotropic mass flux directly on the coarse-fine interface faces.
 *
 * set_2d_cf_bcs writes a velocity, which the solver then turns back into a flux using its own
 * depth at a different time index, so the flux it actually carries is not the one imposed.
 * Writing the flux itself makes the transport exact whatever the depth does.
 *
 * @param[in]     lev       level of refinement
 * @param[in]     time      simulation time to interpolate the parent's flux to
 * @param[inout]  mf_DUon   barotropic u-flux
 * @param[inout]  mf_DVom   barotropic v-flux
 */
void
REMORA::set_2d_cf_flux (int lev, Real time, MultiFab& mf_DUon, MultiFab& mf_DVom)
{
    if (lev == 0 || cf_set_width < 0) { return; }

    BL_PROFILE("REMORA::set_2d_cf_flux()");

    const int set_mask = FPr_Dubar[lev-1].GetSetMaskVal();

    MultiFab Dubar_cf(vec_Dubar_new[lev]->boxArray(), vec_Dubar_new[lev]->DistributionMap(),
                      1, vec_Dubar_new[lev]->nGrowVect());
    MultiFab Dvbar_cf(vec_Dvbar_new[lev]->boxArray(), vec_Dvbar_new[lev]->DistributionMap(),
                      1, vec_Dvbar_new[lev]->nGrowVect());
    Dubar_cf.setVal(zero);
    Dvbar_cf.setVal(zero);

    FPr_Dubar[lev-1].FillSet(Dubar_cf, time, no_bc, domain_bcs_type);
    FPr_Dvbar[lev-1].FillSet(Dvbar_cf, time, no_bc, domain_bcs_type);

    // an x-face varies along y across the interface, and a y-face along x
    if (cf_flux_pc) {
        const IntVect rr = refRatio(lev-1);
        group_average_faces(Dubar_cf, *FPr_Dubar[lev-1].GetMask(), set_mask, 1, rr[1]);
        group_average_faces(Dvbar_cf, *FPr_Dvbar[lev-1].GetMask(), set_mask, 0, rr[0]);
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(mf_DUon, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real      > const& DUon  = mf_DUon.array(mfi);
        Array4<Real const> const& Dubar = Dubar_cf.const_array(mfi);
        Array4<int  const> const& cmask = FPr_Dubar[lev-1].GetMask()->const_array(mfi);
        Array4<Real const> const& pn    = vec_pn[lev]->const_array(mfi);
        Array4<Real const> const& msku  = vec_msku[lev]->const_array(mfi);

        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            if (cmask(i,j,0) != set_mask) { return; }

            Real on_u = two / (pn(i,j,0) + pn(i-1,j,0));
            DUon(i,j,0) = Dubar(i,j,0) * on_u * msku(i,j,0);
        });
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(mf_DVom, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real      > const& DVom  = mf_DVom.array(mfi);
        Array4<Real const> const& Dvbar = Dvbar_cf.const_array(mfi);
        Array4<int  const> const& cmask = FPr_Dvbar[lev-1].GetMask()->const_array(mfi);
        Array4<Real const> const& pm    = vec_pm[lev]->const_array(mfi);
        Array4<Real const> const& mskv  = vec_mskv[lev]->const_array(mfi);

        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            if (cmask(i,j,0) != set_mask) { return; }

            Real om_v = two / (pm(i,j,0) + pm(i,j-1,0));
            DVom(i,j,0) = Dvbar(i,j,0) * om_v * mskv(i,j,0);
        });
    }
}
