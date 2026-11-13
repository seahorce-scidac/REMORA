#include <REMORA.H>

using namespace amrex;

namespace {

/** Strain term at a rho point, shared by the x-flux of u and the y-flux of v there */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real strain_r (int i, int j, int k,
               Array4<Real const> const& Hz, Array4<Real const> const& metr,
               Array4<Real const> const& pm, Array4<Real const> const& pn,
               Array4<Real const> const& uold, Array4<Real const> const& vold, int nrhs) noexcept
{
    return Real(0.5)*Hz(i,j,k) * ( metr(i,j,0,MetricR::pmon_r) *
                                ( (pn(i  ,j,0) + pn(i+1,j,0)) * uold(i+1,j,k,nrhs)-
                                  (pn(i-1,j,0) + pn(i  ,j,0)) * uold(i  ,j,k,nrhs) )-
                                            metr(i,j,0,MetricR::pnom_r) *
                                ( (pm(i,j  ,0) + pm(i,j+1,0)) * vold(i,j+1,k,nrhs)-
                                  (pm(i,j-1,0) + pm(i,j  ,0)) * vold(i,j  ,k,nrhs) ) );
}

/** Strain term at a psi point, shared by the y-flux of u and the x-flux of v there */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real strain_p (int i, int j, int k,
               Array4<Real const> const& Hz, Array4<Real const> const& metp,
               Array4<Real const> const& mskp,
               Array4<Real const> const& pm, Array4<Real const> const& pn,
               Array4<Real const> const& uold, Array4<Real const> const& vold, int nrhs) noexcept
{
    const Real pmon_p = metp(i,j,0,MetricP::pmon_p);
    const Real pnom_p = metp(i,j,0,MetricP::pnom_p);
    return mskp(i,j,0) * Real(0.125) *
                (Hz(i-1,j  ,k) + Hz(i,j ,k)+ Hz(i-1,j-1,k) + Hz(i,j-1,k))*
                (pmon_p*
                 ((pn(i  ,j-1,0)+pn(i  ,j,0))*vold(i  ,j,k,nrhs)-
                  (pn(i-1,j-1,0)+pn(i-1,j,0))*vold(i-1,j,k,nrhs))+
                 pnom_p*
                 ((pm(i-1,j  ,0)+pm(i,j  ,0))*uold(i,j  ,k,nrhs)-
                  (pm(i-1,j-1,0)+pm(i,j-1,0))*uold(i,j-1,k,nrhs)));
}

}

/**
 * @param[in   ] xbx       x-nodal tile box
 * @param[in   ] ybx       y-nodal tile box
 * @param[inout] u         u-direction velocity
 * @param[inout] v         v-direction velocity
 * @param[in   ] uold      u-direction velocity at last step
 * @param[in   ] vold      v-direction velocity at last step
 * @param[inout] rufrc     forcing term for u-velocity RHS
 * @param[inout] rvfrc     forcing term for v-velocity RHS
 * @param[in   ] visc2_p   harmonic viscosity on psi points
 * @param[in   ] visc2_r   harmonic viscosity on rho points
 * @param[in   ] Hz        vertical cell height
 * @param[in   ] pm        1/dx
 * @param[in   ] pn        1/dy
 * @param[in   ] metr      metric products at rho points (MetricR)
 * @param[in   ] metp      metric products at psi points (MetricP)
 * @param[in   ] mskp      land-sea mask on psi-points
 * @param[in   ] nrhs      index of RHS component
 * @param[in   ] nnew      index of time step to update
 * @param[in   ] dt_lev    time step at this level
 */
void
REMORA::uv3dmix  (const Box& xbx, const Box& ybx,
                 const Array4<Real      >& u,
                 const Array4<Real      >& v,
                 const Array4<Real const>& uold,
                 const Array4<Real const>& vold,
                 const Array4<Real      >& rufrc,
                 const Array4<Real      >& rvfrc,
                 const Array4<Real const>& visc2_p,
                 const Array4<Real const>& visc2_r,
                 const Array4<Real const>& Hz,
                 const Array4<Real const>& pm,
                 const Array4<Real const>& pn,
                 const Array4<Real const>& metr,
                 const Array4<Real const>& metp,
                 const Array4<Real const>& mskp,
                 int nrhs, int nnew,
                 const Real dt_lev)
{
    BL_PROFILE("REMORA::uv3dmix()");
    //-----------------------------------------------------------------------
    //  Add in harmonic viscosity s terms.
    //-----------------------------------------------------------------------

    int ncomp = 0;
    int UFx_comp = ncomp++;
    int UFe_comp = ncomp++;
    int VFe_comp = ncomp++;
    int VFx_comp = ncomp++;

    Box bx = enclosedCells(xbx);
    // Each flux is read only on the box its kernel below writes
    FArrayBox fab(grow(bx,IntVect(2,2,0)),ncomp,amrex::The_Async_Arena());

    auto UFx=fab.array(UFx_comp);
    auto UFe=fab.array(UFe_comp);
    auto VFx=fab.array(VFx_comp);
    auto VFe=fab.array(VFe_comp);

    auto N = xbx.hiVect()[2] - ybx.loVect()[2];

    // Think of the cell-centered box as                              [ 0:nx-1, 0:ny-1] (0,0,0) cc
    //
    // xbx is the x-face-centered box on which we update u  [ 0:nx  , 0:ny-1] (1,0,0) x-faces
    //   to do so requires                        UFx on    [-1:nx  , 0:ny-1] (0,0,0) cc
    //      which requires                        uold on   [-1:nx+1, 0:ny-1] (1,0,0) x-faces
    //        and requires                        vold on   [-1:nx  , 0:ny  ] (0,1,0) y-faces
    //   to do so requires                        UFe on    [ 0:nx  , 0:ny  ] (1,1,0) xy-nodes
    //      which requires                        uold on   [ 0:nx  ,-1:ny  ] (1,0,0) x-faces
    //       and  requires                        vold on   [ -1:nx , 0:ny  ] (0,1,0) y-faces
    //
    // ybx is the y-face-centered box on which we update v  [ 0:nx-1, 0:ny  ] (0,1,0) y-faces
    //   to do so requires                        VFe on    [ 0:nx-1,-1:ny  ] (0,0,0) cc
    //   to do so requires                        VFx on    [ 0:nx  , 0:ny  ] (1,1,0) xy-nodes
    //
    // The x-flux of u and the y-flux of v at a rho point differ only in their metric
    // factor, as do the y-flux of u and the x-flux of v at a psi point, so each pair is
    // computed together from one strain term. At rho points the pair is computed on the
    // box of UFx, which holds the box of VFe except for the rows below and above it; at psi
    // points on the box of UFe, which holds the box of VFx except for the column to its
    // right. Those rows and that column get VFe or VFx alone. Where only one flux of a pair
    // is read, the other is computed but not read.
    const Box rbx  = growLo(xbx,0,1);
    const Box rbxv = growLo(ybx,1,1);
    const Box pbx  = growHi(xbx,1,1);
    const Box pbxv = growHi(ybx,0,1);
    AMREX_ASSERT(rbxv.smallEnd(0) >= rbx.smallEnd(0) && rbxv.bigEnd(0) <= rbx.bigEnd(0));
    AMREX_ASSERT(pbxv.smallEnd(1) >= pbx.smallEnd(1) && pbxv.bigEnd(1) <= pbx.bigEnd(1));
    AMREX_ASSERT(pbxv.smallEnd(0) == pbx.smallEnd(0));
    Box rbxv_lo = rbxv; rbxv_lo.setBig  (1, rbx.smallEnd(1)-1);
    Box rbxv_hi = rbxv; rbxv_hi.setSmall(1, rbx.bigEnd(1)+1);
    Box pbxv_hi = pbxv; pbxv_hi.setSmall(0, pbx.bigEnd(0)+1);

    ParallelFor(rbx, rbxv_lo, rbxv_hi,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const Real cff = strain_r(i,j,k,Hz,metr,pm,pn,uold,vold,nrhs);
        Real on_r = metr(i,j,0,MetricR::on_r);
        Real om_r = metr(i,j,0,MetricR::om_r);
        UFx(i,j,k) = on_r * on_r * visc2_r(i,j,0) * cff;
        VFe(i,j,k) = om_r * om_r * visc2_r(i,j,0) * cff;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const Real cff = strain_r(i,j,k,Hz,metr,pm,pn,uold,vold,nrhs);
        Real om_r = metr(i,j,0,MetricR::om_r);
        VFe(i,j,k) = om_r * om_r * visc2_r(i,j,0) * cff;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const Real cff = strain_r(i,j,k,Hz,metr,pm,pn,uold,vold,nrhs);
        Real om_r = metr(i,j,0,MetricR::om_r);
        VFe(i,j,k) = om_r * om_r * visc2_r(i,j,0) * cff;
    });

    ParallelFor(pbx, pbxv_hi,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const Real cff = strain_p(i,j,k,Hz,metp,mskp,pm,pn,uold,vold,nrhs);
        const Real om_p = metp(i,j,0,MetricP::om_p);
        const Real on_p = metp(i,j,0,MetricP::on_p);
        UFe(i,j,k) = om_p*om_p*visc2_p(i,j,0)*cff;
        VFx(i,j,k) = on_p*on_p*visc2_p(i,j,0)*cff;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const Real cff = strain_p(i,j,k,Hz,metp,mskp,pm,pn,uold,vold,nrhs);
        const Real on_p = metp(i,j,0,MetricP::on_p);
        VFx(i,j,k) = on_p*on_p*visc2_p(i,j,0)*cff;
    });

    // Update u and v, and add the depth integral of the viscous term to rufrc and rvfrc,
    // level by level from the bottom up
    ParallelFor(makeSlab(xbx,2,0), N+1,
    [=] AMREX_GPU_DEVICE (int i, int j, int, int k)
    {
        const Real cff=dt_lev*Real(0.25)*(pm(i-1,j,0)+pm(i,j,0))*(pn(i-1,j,0)+pn(i,j,0));
        const Real cff1=Real(0.5)*(pn(i-1,j,0)+pn(i,j,0))*(UFx(i,j  ,k)-UFx(i-1,j,k));
        const Real cff2=Real(0.5)*(pm(i-1,j,0)+pm(i,j,0))*(UFe(i,j+1,k)-UFe(i  ,j,k));
        const Real cff3=cff*(cff1+cff2);
        u(i,j,k,nnew)=u(i,j,k,nnew)+cff3;
        rufrc(i,j,0) += cff1+cff2;
    },
    makeSlab(ybx,2,0), N+1,
    [=] AMREX_GPU_DEVICE (int i, int j, int, int k)
    {
        const Real cff=dt_lev*Real(0.25)*(pm(i,j,0)+pm(i,j-1,0))*(pn(i,j,0)+pn(i,j-1,0));
        const Real cff1=Real(0.5)*(pn(i,j-1,0)+pn(i,j,0))*(VFx(i+1,j,k)-VFx(i,j  ,k));
        const Real cff2=Real(0.5)*(pm(i,j-1,0)+pm(i,j,0))*(VFe(i  ,j,k)-VFe(i,j-1,k));
        const Real cff3=cff*(cff1-cff2);

        v(i,j,k,nnew)=v(i,j,k,nnew)+cff3;
        rvfrc(i,j,0) += cff1-cff2;
    });
}
