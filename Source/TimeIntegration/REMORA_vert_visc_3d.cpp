#include <REMORA.H>

using namespace amrex;

namespace {

/** Average of a cell-centered field onto the u-face (IOFF = 1), the v-face (JOFF = 1), or
 *  its value at the cell center for tracers (both 0), where the average of a value with
 *  itself would return it unchanged. */
template <int IOFF, int JOFF>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real face_avg (Array4<Real const> const& a, int i, int j, int k) noexcept
{
    if constexpr (IOFF == 0 && JOFF == 0) {
        return a(i,j,k);
    } else {
        return Real(0.5) * (a(i-IOFF,j-JOFF,k) + a(i,j,k));
    }
}

}

/**
 * Time-steps vertical viscosity (or diffusivity, for tracers) implicitly with a
 * tridiagonal solve in each column. Corresponds to the vertical viscosity part of
 * step3d_uv.F and the vertical diffusion part of step3d_t.F in ROMS.
 *
 * @param[in   ] phi_bx     box to update on
 * @param[in   ] ioff       x-direction offset: 1 for u, 0 otherwise
 * @param[in   ] joff       y-direction offset: 1 for v, 0 otherwise
 * @param[inout] phi        velocity (u or v) or tracer
 * @param[in   ] Hz         vertical height of cells
 * @param[in   ] Akv        vertical viscosity or diffusivity coefficient
 * @param[in   ] nnew       index of time step to update
 * @param[in   ] N          number of vertical levels
 * @param[in   ] dt_lev     time step at this refinement level
 */
void
REMORA::vert_visc_3d (const Box& phi_bx, const int ioff, const int joff,
                     const Array4<Real      >& phi,
                     const Array4<Real const>& Hz,
                     const Array4<Real const>& Akv,
                     const int nnew, const int N, const Real dt_lev)
{
    BL_PROFILE("REMORA::vert_visc_3d()");
    if (verbose > 1)
        amrex::Print() << "updating on box in vert_visc_3d: " << phi_bx << std::endl;

    AMREX_ASSERT(phi_bx.smallEnd(2) == 0 && phi_bx.bigEnd(2) == N);

    // Tile-private scratch for the elimination: CF and DC on the interior w-levels 1..N
    FArrayBox fab(phi_bx, 2, The_Async_Arena());
    auto const& CF = fab.array(0);
    auto const& DC = fab.array(1);

    const Box bx2d = makeSlab(phi_bx,2,0);
    const Real sixth = one / Real(6.0);

    // Where the face lies: 0 for a tracer at the cell center, 1 for u, 2 for v
    AMREX_ASSERT(ioff + joff <= 1);
    const int face = (ioff == 1) ? 1 : ((joff == 1) ? 2 : 0);

    // The face values of Hz and Akv are formed where they are used rather than stored.
    // Both kernels walk each column in order: on the GPU one thread per column loops over
    // the levels, on the CPU the level loop is outermost, so the loop over i vectorizes.

    //
    //  Use conservative, parabolic spline reconstruction of vertical
    //  viscosity derivatives. Then, time step vertical viscosity term
    //  implicitly by solving a tridiagonal system: forward elimination over the
    //  interior w-levels k = 1..N, with zero flux at the bottom (k = 0).
    //
    AnyCTO(TypeList<CompileTimeOptions<0,1,2>>{}, {face},
    [&] (auto const& f) { ParallelFor(bx2d, N, f); },
    [=] AMREX_GPU_DEVICE (int i, int j, int, int kk, auto face_c)
    {
        constexpr int IO = (face_c.value == 1) ? 1 : 0;
        constexpr int JO = (face_c.value == 2) ? 1 : 0;
        const int k = kk+1;
        const Real Hzkm1 = face_avg<IO,JO>(Hz,i,j,k-1);
        const Real Hzk   = face_avg<IO,JO>(Hz,i,j,k  );
        const Real AKm1  = face_avg<IO,JO>(Akv,i,j,k-1);
        const Real AK    = face_avg<IO,JO>(Akv,i,j,k  );
        const Real AKp1  = face_avg<IO,JO>(Akv,i,j,k+1);

        const Real oHzkm1 = one / Hzkm1;
        const Real oHz    = one / Hzk;

        const Real FC  = sixth * Hzkm1 - dt_lev * AKm1 / Hzkm1;
        const Real CFk = sixth * Hzk   - dt_lev * AKp1 / Hzk;
        const Real BC  = third * (Hzkm1 + Hzk) + dt_lev * AK * (oHzkm1 + oHz);

        const Real CFm1 = (k == 1) ? zero : CF(i,j,k-1);
        const Real DCm1 = (k == 1) ? zero : DC(i,j,k-1);

        const Real cff = one / (BC - FC * CFm1);
        CF(i,j,k) = CFk * cff;
        DC(i,j,k) = cff * (phi(i,j,k,nnew) - phi(i,j,k-1,nnew) - FC * DCm1);
    });

    //
    //  Backward substitution from the top (zero flux at k = N+1), then the update of
    //  cell k from the fluxes AK*DC at its two faces, as soon as both are final.
    //
    AnyCTO(TypeList<CompileTimeOptions<0,1,2>>{}, {face},
    [&] (auto const& f) { ParallelFor(bx2d, N+1, f); },
    [=] AMREX_GPU_DEVICE (int i, int j, int, int kk, auto face_c)
    {
        constexpr int IO = (face_c.value == 1) ? 1 : 0;
        constexpr int JO = (face_c.value == 2) ? 1 : 0;
        const int k = N - kk;
        const Real DCp1 = (k == N) ? zero : DC(i,j,k+1);
        Real DCk = zero;
        if (k >= 1) {
            DCk = DC(i,j,k) - CF(i,j,k) * DCp1;
            DC(i,j,k) = DCk;
        }
        const Real flux_top = DCp1 * face_avg<IO,JO>(Akv,i,j,k+1);
        const Real flux_bot = DCk  * face_avg<IO,JO>(Akv,i,j,k  );

        const Real oHz = one / face_avg<IO,JO>(Hz,i,j,k);
        phi(i,j,k) += dt_lev * oHz * (flux_top - flux_bot);
    });
}
