#include <REMORA.H>

using namespace amrex;

/**
 * Vertical volume flux across the w-levels (ROMS omega times the cell area), from the
 * vertical integral of the horizontal mass flux divergence. W is zero at the bed and at
 * the free surface.
 *
 * @param[in   ] bx2d   k = 0 slab of the columns to compute
 * @param[  out] W      vertical volume flux
 * @param[in   ] Huon   u-volume flux
 * @param[in   ] Hvom   v-volume flux
 * @param[in   ] z_w    z coordinates at w points
 * @param[in   ] h      bathymetry
 * @param[in   ] N      number of vertical levels minus one
 */
void
REMORA::compute_W (const Box& bx2d,
                   const Array4<Real      >& W,
                   const Array4<Real const>& Huon,
                   const Array4<Real const>& Hvom,
                   const Array4<Real const>& z_w,
                   const Array4<Real const>& h,
                   const int N)
{
    ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        W(i,j,0) = zero;
    });
    ParallelFor(bx2d, N+1,
    [=] AMREX_GPU_DEVICE (int i, int j, int , int kk)
    {
        //  Starting with zero vertical velocity at the bottom, integrate
        //  from the bottom (k=0) to the free-surface (k=N).  The w(:,:,N(ng))
        //  contains the vertical velocity at the free-surface, d(zeta)/d(t).
        //  Notice that barotropic mass flux divergence is not used directly.
        //
        int k = kk + 1;
        W(i,j,k) = W(i,j,k-1) - (Huon(i+1,j,k-1)-Huon(i,j,k-1)) - (Hvom(i,j+1,k-1)-Hvom(i,j,k-1));
    });
    ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        W(i,j,N+1)=W(i,j,N+1)/(z_w(i,j,N+1)+h(i,j,0,0)); // wrk_i
    });
    ParallelFor(bx2d, N,
    [=] AMREX_GPU_DEVICE (int i, int j, int , int kk)
    {
        int k = kk + 1;
        W(i,j,k) = W(i,j,k)- W(i,j,N+1)*(z_w(i,j,k)+h(i,j,0,0));
    });
    ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        W(i,j,N+1) = zero;
    });
}

/**
 * W converted from a volume flux to a velocity, as ROMS scale_omega. It is zero at the bed
 * and the free surface.
 *
 * @param[in ] lev       level
 * @param[out] mf_omega  z-nodal omega (m/s); ghost cells are filled across grids and
 *                       periodic boundaries, and zero outside the domain
 * @param[in ] mf_Huon   u-volume flux on at least the valid faces
 * @param[in ] mf_Hvom   v-volume flux on at least the valid faces
 */
void
REMORA::compute_omega (int lev, MultiFab& mf_omega,
                       const MultiFab& mf_Huon, const MultiFab& mf_Hvom)
{
    const int N = Geom(lev).Domain().size()[2]-1;


    // Not tiled: compute_W integrates whole columns
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf_omega); mfi.isValid(); ++mfi)
    {
        Array4<Real> const& W = mf_omega.array(mfi);

        Array4<Real const> const& pm = vec_pm[lev]->const_array(mfi);
        Array4<Real const> const& pn = vec_pn[lev]->const_array(mfi);

        Box bx2d = makeSlab(mfi.validbox(),2,0);
        compute_W(bx2d, W, mf_Huon.const_array(mfi), mf_Hvom.const_array(mfi),
                  vec_z_w[lev]->const_array(mfi), vec_h[lev]->const_array(mfi), N);

        ParallelFor(bx2d, N, [=] AMREX_GPU_DEVICE (int i, int j, int, int kk)
        {
            W(i,j,kk+1) *= pm(i,j,0) * pn(i,j,0);
        });
    }

    mf_omega.FillBoundary(geom[lev].periodicity());
}

/**
 * omega as ROMS writes it to its history file: main3d calls set_massflux on the current
 * velocities and then omega, before output. The fluxes left in vec_Huon and vec_Hvom by
 * update_massflux_3d are instead time-averaged over the step (the ones the tracers move
 * with), so they are not used here.
 *
 * @param[in ] lev       level
 * @param[out] mf_omega  z-nodal omega (m/s)
 */
void
REMORA::compute_diagnostic_omega (int lev, MultiFab& mf_omega)
{
    MultiFab mf_Huon(convert(grids[lev],IntVect(1,0,0)), dmap[lev], 1, IntVect(0));
    MultiFab mf_Hvom(convert(grids[lev],IntVect(0,1,0)), dmap[lev], 1, IntVect(0));
    set_massflux_3d(lev, mf_Huon, mf_Hvom);
    compute_omega(lev, mf_omega, mf_Huon, mf_Hvom);
}

/**
 * The "true" vertical velocity w (m/s) at w-points, a port of ROMS wvelocity. In
 * terrain-following coordinates
 *
 *     Hz * omega = w - d(z)/d(t) - div(z),   div(z) = pm * u * d(z)/d(xi) + pn * v * d(z)/d(eta)
 *
 * so w adds back to omega the stretching of the levels by the moving free surface
 * (from the barotropic flux divergence) and the motion along the sloping levels, which
 * is computed at rho-points and shifted to w-points with ROMS's cubic interpolation.
 * Like ROMS, this is for output only: nothing in the dynamics reads it, and omega comes
 * from the mass fluxes of the current velocities (compute_diagnostic_omega).
 *
 * @param[in ] lev     level
 * @param[out] mf_w    z-nodal w on the valid region
 */
void
REMORA::compute_wvel (int lev, MultiFab& mf_w)
{
    MultiFab mf_omega(convert(grids[lev],IntVect(0,0,1)), dmap[lev], 1, IntVect(0));
    compute_diagnostic_omega(lev, mf_omega);
    compute_wvel(lev, mf_w, mf_omega);
}

/**
 * @param[in ] lev       level
 * @param[out] mf_w      z-nodal w on the valid region
 * @param[in ] mf_omega  z-nodal omega from compute_diagnostic_omega, on the valid region
 */
void
REMORA::compute_wvel (int lev, MultiFab& mf_w, const MultiFab& mf_omega)
{
    const int nz = Geom(lev).Domain().length(2);
    // The cubic stencil reaches three rho-levels in from each end
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nz >= 3, "compute_wvel needs at least 3 vertical levels");

    const Real cff1 = Real(3.0)/Real(8.0);
    const Real cff2 = Real(3.0)/Real(4.0);
    const Real cff3 = Real(1.0)/Real(8.0);
    const Real cff4 = Real(9.0)/Real(16.0);
    const Real cff5 = Real(1.0)/Real(16.0);

    // Not tiled: each column of w reads the whole column of vert
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf_omega); mfi.isValid(); ++mfi)
    {
        Array4<Real      > const& wvel  = mf_w.array(mfi);
        Array4<Real const> const& omega = mf_omega.const_array(mfi);

        Array4<Real const> const& u    = xvel_new[lev]->const_array(mfi);
        Array4<Real const> const& v    = yvel_new[lev]->const_array(mfi);
        Array4<Real const> const& z_r  = vec_z_r[lev]->const_array(mfi);
        Array4<Real const> const& z_w  = vec_z_w[lev]->const_array(mfi);
        Array4<Real const> const& pm   = vec_pm[lev]->const_array(mfi);
        Array4<Real const> const& pn   = vec_pn[lev]->const_array(mfi);
        Array4<Real const> const& DU_avg1 = vec_DU_avg1[lev]->const_array(mfi);
        Array4<Real const> const& DV_avg1 = vec_DV_avg1[lev]->const_array(mfi);

        const Box bx = mfi.validbox().enclosedCells(2);
        FArrayBox fab_vert(bx,1,amrex::The_Async_Arena());
        auto vert = fab_vert.array();

        // Motion along the s-levels, (Ui + Vj) * GRADs(z), at rho-points
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            const Real wrk_ulo = u(i  ,j,k) * (z_r(i  ,j,k)-z_r(i-1,j,k)) * (pm(i-1,j,0)+pm(i  ,j,0));
            const Real wrk_uhi = u(i+1,j,k) * (z_r(i+1,j,k)-z_r(i  ,j,k)) * (pm(i  ,j,0)+pm(i+1,j,0));
            const Real wrk_vlo = v(i,j  ,k) * (z_r(i,j  ,k)-z_r(i,j-1,k)) * (pn(i,j-1,0)+pn(i,j  ,0));
            const Real wrk_vhi = v(i,j+1,k) * (z_r(i,j+1,k)-z_r(i,j  ,k)) * (pn(i,j  ,0)+pn(i,j+1,0));
            vert(i,j,k) = Real(0.25)*(wrk_ulo+wrk_uhi) + Real(0.25)*(wrk_vlo+wrk_vhi);
        });

        // ROMS k (w-levels 0..N, rho-levels 1..N) maps to REMORA k (w-levels 0..nz,
        // rho-levels 0..nz-1), so ROMS vert(k) is vert(k-1) here
        ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            // Free-surface stretching: d(zeta)/d(t) spread linearly through the column
            const Real wrk = (DU_avg1(i,j,0)-DU_avg1(i+1,j,0) + DV_avg1(i,j,0)-DV_avg1(i,j+1,0))
                           / (z_w(i,j,nz)-z_w(i,j,0));
            const Real pmpn = pm(i,j,0)*pn(i,j,0);
            auto interior = [&] (int k) {
                return omega(i,j,k) + pmpn*wrk*(z_w(i,j,k)-z_w(i,j,0));
            };

            Real slope = (z_r(i,j,0)-z_w(i,j,0)) / (z_r(i,j,1)-z_r(i,j,0)); // extrapolation slope
            wvel(i,j,0) = cff1*(vert(i,j,0) - slope*(vert(i,j,1)-vert(i,j,0)))
                        + cff2*vert(i,j,0) - cff3*vert(i,j,1);
            wvel(i,j,1) = interior(1) + cff1*vert(i,j,0) + cff2*vert(i,j,1) - cff3*vert(i,j,2);
            for (int k = 2; k <= nz-2; ++k) {
                wvel(i,j,k) = interior(k) + cff4*(vert(i,j,k-1)+vert(i,j,k  ))
                                          - cff5*(vert(i,j,k-2)+vert(i,j,k+1));
            }
            slope = (z_w(i,j,nz)-z_r(i,j,nz-1)) / (z_r(i,j,nz-1)-z_r(i,j,nz-2)); // extrapolation slope
            wvel(i,j,nz) = pmpn*wrk*(z_w(i,j,nz)-z_w(i,j,0))
                         + cff1*(vert(i,j,nz-1) + slope*(vert(i,j,nz-1)-vert(i,j,nz-2)))
                         + cff2*vert(i,j,nz-1) - cff3*vert(i,j,nz-2);
            wvel(i,j,nz-1) = interior(nz-1) + cff1*vert(i,j,nz-1) + cff2*vert(i,j,nz-2)
                                            - cff3*vert(i,j,nz-3);
        });
    }
}

/**
 * Horizontal mass fluxes Hz*u/n and Hz*v/m of the current velocities on the valid faces,
 * as ROMS set_massflux.
 *
 * @param[in ] lev       level
 * @param[out] mf_Huon   u-volume flux
 * @param[out] mf_Hvom   v-volume flux
 */
void
REMORA::set_massflux_3d (int lev, MultiFab& mf_Huon, MultiFab& mf_Hvom)
{
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf_Huon, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real      > const& Huon = mf_Huon.array(mfi);
        Array4<Real const> const& u    = xvel_new[lev]->const_array(mfi);
        Array4<Real const> const& Hz   = vec_Hz[lev]->const_array(mfi);
        Array4<Real const> const& pn   = vec_pn[lev]->const_array(mfi);

        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real on_u = two / (pn(i-1,j,0)+pn(i,j,0));
            Huon(i,j,k)=Real(0.5)*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k)* on_u;
        });
    }
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf_Hvom, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real      > const& Hvom = mf_Hvom.array(mfi);
        Array4<Real const> const& v    = yvel_new[lev]->const_array(mfi);
        Array4<Real const> const& Hz   = vec_Hz[lev]->const_array(mfi);
        Array4<Real const> const& pm   = vec_pm[lev]->const_array(mfi);

        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real om_v= two / (pm(i,j-1,0)+pm(i,j,0));
            Hvom(i,j,k)=Real(0.5)*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k)* om_v;
        });
    }
}
