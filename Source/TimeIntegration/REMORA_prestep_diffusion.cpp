#include <REMORA.H>

using namespace amrex;

/**
 * Called from prestep. The tracer update is a bit different from the u,v updates so we test
 * for it, but checking if ioff=0 and joff=0. In some cases, though, we
 * can recover the tracer update from the generic one by setting those indices.
 * Setting icc and ntfirst identically for the tracers should be equivalent
 * to setting ioff=0 and joff=0
 *
 * @param[in   ] vel_bx   tile box
 * @param[in   ] gbx      grown tile box
 * @param[in   ] ioff     offset in x direction
 * @param[in   ] joff     offset in y direction
 * @param[  out] vel      velocity or scalar to update
 * @param[in   ] vel_old  velocity or scalar at last time
 * @param[inout] rvel     velocity or scalar RHS
 * @param[in   ] Hz       vertical cell height
 * @param[in   ] Akv      vertical viscosity coefficient
 * @param[in   ] sstr     surface flux
 * @param[in   ] bstr     bottom flux
 * @param[in   ] z_r      z coordinates at rho points
 * @param[in   ] pm       1/dx
 * @param[in   ] pn       1/dy
 * @param[in   ] iic      which time step we're on
 * @param[in   ] ntfirst  what is the first time step?
 * @param[in   ] nnew     index of time step to update
 * @param[in   ] nstp     index of last time step
 * @param[in   ] nrhs     index of RHS component
 * @param[in   ] N        number of vertical levels
 * @param[in   ] lambda   weighting coefficient for the newest (implicit) time step derivatives
 * @param[in   ] dt_lev   time step at this level
 */

void
REMORA::prestep_diffusion (const Box& vel_bx, const Box& gbx,
                          const int ioff, const int joff,
                          const Array4<Real      >& vel,
                          const Array4<Real const>& vel_old,
                          const Array4<Real      >& rvel,
                          const Array4<Real const>& Hz,
                          const Array4<Real const>& Akv,
                          const Array4<Real const>& sstr,
                          const Array4<Real const>& bstr,
                          const Array4<Real const>& z_r,
                          const Array4<Real const>& pm,
                          const Array4<Real const>& pn,
                          const int iic, const int ntfirst, const int nnew, int nstp, int nrhs, int N,
                          const Real lambda, const Real dt_lev)
{
    BL_PROFILE("REMORA::prestep_diffusion()");

    BoxArray ba_gbxvel = intersect(BoxArray(vel_bx), gbx);
    AMREX_ASSERT((ba_gbxvel.size() == 1));
    Box gbxvel = ba_gbxvel[0];
    AMREX_ASSERT(gbxvel.smallEnd(2) == 0 && gbxvel.bigEnd(2) == N);

    //
    //  Weighting coefficient for the newest (implicit) time step derivatives
    //  using either a Crank-Nicolson implicit scheme (lambda=Real(0.5)) or a
    //  backward implicit scheme (lambda=one).
    //

    Real oml_dt = dt_lev*(one-lambda);
    //N is one less than ROMS

    //  Except the commented out part means lambda is always one
    if (verbose > 1) {
        amrex::Print() << "in update_vel_3d with box " << vel_bx << std::endl;
        Print() << "vel old " << Box(vel_old) << std::endl;
        Print() << "Akv " << Box(Akv) << std::endl;
    }

    // Where the variable lies: 0 at the cell center (tracers), 1 on the u-face, 2 on the v-face
    AMREX_ASSERT(ioff + joff <= 1);
    const int face = (ioff == 1) ? 1 : ((joff == 1) ? 2 : 0);
    // Which time-stepping update: the first step (and every tracer step), the second, or the
    // later ones
    const int scheme = (iic==ntfirst || face == 0) ? 0 : ((iic==ntfirst+1) ? 1 : 2);

    // The vertical flux FC at the w-level below the cell being updated, carried up the column
    // from one level to the next: FC(k) = dt*bstr at the bottom (k = 0), the diffusive flux
    // in the interior and dt*sstr at the surface (k = N+1).
    const Box bx2d = makeSlab(gbxvel,2,0);
    FArrayBox fab_FC(bx2d,1,The_Async_Arena());
    auto const& FCb = fab_FC.array();

    AnyCTO(TypeList<CompileTimeOptions<0,1,2>,CompileTimeOptions<0,1,2>>{}, {face, scheme},
    [&] (auto const& f) { ParallelFor(bx2d, N+1, f); },
    [=] AMREX_GPU_DEVICE (int i, int j, int, int k, auto face_c, auto scheme_c)
    {
        constexpr int io = (face_c.value == 1) ? 1 : 0;
        constexpr int jo = (face_c.value == 2) ? 1 : 0;

        const Real FC_lo = (k == 0) ? dt_lev*bstr(i,j,0) : FCb(i,j,0);
        Real FC_hi;
        if (k == N) {
            FC_hi = dt_lev*sstr(i,j,0);
        } else {
            const int kp = k+1;
            Real cff = one / ( z_r(i,j,kp)+z_r(i-io,j-jo,kp)
                              -z_r(i,j,kp-1)-z_r(i-io,j-jo,kp-1));
            FC_hi = oml_dt * cff * (vel_old(i,j,kp,nstp)-vel_old(i,j,kp-1,nstp)) *
                                   (Akv(i,j,kp)     +Akv(i-io,j-jo,kp));
            FCb(i,j,0) = FC_hi;
        }

        if (scheme_c.value == 0) {
            Real cff1=vel_old(i,j,k,nstp)*Real(0.5)*(Hz(i,j,k)+Hz(i-io,j-jo,k));
            Real cff2=FC_hi-FC_lo;
            vel(i,j,k,nnew)=cff1+cff2;
        } else if (scheme_c.value == 1) {
            Real cff1=vel_old(i,j,k,nstp)*Real(0.5)*(Hz(i,j,k)+Hz(i-io,j-jo,k));
            Real cff2=FC_hi-FC_lo;
            Real DC = Real(0.25) * dt_lev * (pm(i,j,0)+pm(i-io,j-jo,0))
                                  * (pn(i,j,0)+pn(i-io,j-jo,0));
            Real cff3 = DC * Real(0.5);
            int indx=nrhs ? 0 : 1;
            Real r_swap= rvel(i,j,k,indx);
            rvel(i,j,k,indx) = rvel(i,j,k,nrhs);
            rvel(i,j,k,nrhs) = r_swap;
            vel(i,j,k,nnew)=cff1- cff3*rvel(i,j,k,indx)+ cff2;
        } else {
            Real cff1 =  Real(5.0)/Real(12.0);
            Real cff2 = Real(16.0)/Real(12.0);
            Real cff3=vel_old(i,j,k,nstp)*Real(0.5)*(Hz(i,j,k)+Hz(i-io,j-jo,k));
            Real cff4=FC_hi-FC_lo;
            Real DC = Real(0.25) * dt_lev * (pm(i,j,0)+pm(i-io,j-jo,0))
                                  * (pn(i,j,0)+pn(i-io,j-jo,0));

            int indx=nrhs ? 0 : 1;
            Real r_swap= rvel(i,j,k,indx);
            rvel(i,j,k,indx) = rvel(i,j,k,nrhs);
            rvel(i,j,k,nrhs) = r_swap;

            vel(i,j,k,nnew) = cff3 + DC*(cff1*rvel(i,j,k,nrhs)-
                                                cff2*rvel(i,j,k,indx))+cff4;
            rvel(i,j,k,nrhs) = zero;
        }
    });
}
