#include <REMORA.H>
#include <REMORA_PhysBCFunct.H>

using namespace amrex;

namespace {
    PhysBCFunctNoOp no_bc;
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
 * D comes from know, the current state, as ROMS builds it from indx1; knew is only the slot
 * being written.
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

    MultiFab Dubar_cf(vec_Dubar_new[lev]->boxArray(), vec_Dubar_new[lev]->DistributionMap(),
                      1, vec_Dubar_new[lev]->nGrowVect());
    MultiFab Dvbar_cf(vec_Dvbar_new[lev]->boxArray(), vec_Dvbar_new[lev]->DistributionMap(),
                      1, vec_Dvbar_new[lev]->nGrowVect());
    Dubar_cf.setVal(zero);
    Dvbar_cf.setVal(zero);

    FPr_Dubar[lev-1].FillSet(Dubar_cf, time, no_bc, domain_bcs_type);
    FPr_Dvbar[lev-1].FillSet(Dvbar_cf, time, no_bc, domain_bcs_type);


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

            Real D = Real(0.5) * (h(i-1,j,0,0) + zeta(i-1,j,0,know) +
                                  h(i  ,j,0,0) + zeta(i  ,j,0,know));
            if (D <= zero) { return; }

            ubar(i,j,0,knew) = Dubar(i,j,0) / D * msku(i,j,0);
        });
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

            Real D = Real(0.5) * (h(i,j-1,0,0) + zeta(i,j-1,0,know) +
                                  h(i,j  ,0,0) + zeta(i,j  ,0,know));
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
