#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] bx         box to operate on
 * @param[inout] state      scalar data
 * @param[  out] state_rhs  scalar data RHS
 * @param[in   ] diff2      diffusivity
 * @param[in   ] Hz         vertical cell height
 * @param[in   ] z_r        rho-point z coordinate
 * @param[in   ] pm         1/dx
 * @param[in   ] pn         1/dy
 * @param[in   ] msku       land-sea mask on u-points
 * @param[in   ] mskv       land-sea mask on v-points
 * @param[in   ] dt_lev     time step at level
 * @param[in   ] ncomp      number of components to do this calculation on
 * @param[in   ] N          number of vertical levels
 */
void
REMORA::t3dmix2   (const Box& bx,
                   const Array4<Real      >& state,
                   const Array4<Real      >& state_rhs,
                   const Array4<Real const>& diff2,
                   const Array4<Real const>& Hz,
                   const Array4<Real const>& z_r,
                   const Array4<Real const>& pm,
                   const Array4<Real const>& pn,
                   const Array4<Real const>& msku,
                   const Array4<Real const>& mskv,
                   const Real dt_lev, const int ncomp,
                   const int N) {
    if (solverChoice.harmonic_mixing_type == HarmonicMixingType::s) {
        t3dmix2_s(bx, state, state_rhs, diff2, Hz, pm, pn, msku, mskv, dt_lev, ncomp);
    } else if (solverChoice.harmonic_mixing_type == HarmonicMixingType::geopotential) {
        t3dmix2_geo(bx, state, state_rhs, diff2, Hz, z_r, pm, pn, msku, mskv, dt_lev, ncomp, N);
    }
}


/**
 * @param[in   ] bx         box to operate on
 * @param[inout] state      scalar data
 * @param[  out] state_rhs  scalar data RHS
 * @param[in   ] diff2      diffusivity
 * @param[in   ] Hz         vertical cell height
 * @param[in   ] pm         1/dx
 * @param[in   ] pn         1/dy
 * @param[in   ] msku       land-sea mask on u-points
 * @param[in   ] mskv       land-sea mask on v-points
 * @param[in   ] dt_lev     time step at level
 * @param[in   ] ncomp      number of components to do this calculation on
 */
void
REMORA::t3dmix2_s (const Box& bx,
                   const Array4<Real      >& state,
                   const Array4<Real      >& state_rhs,
                   const Array4<Real const>& diff2,
                   const Array4<Real const>& Hz,
                   const Array4<Real const>& pm,
                   const Array4<Real const>& pn,
                   const Array4<Real const>& msku,
                   const Array4<Real const>& mskv,
                   const Real dt_lev, const int ncomp)
{
    BL_PROFILE("REMORA::t3dmix2_s()");
    //-----------------------------------------------------------------------
    //  Add in harmonic diffusivity s terms.
    //-----------------------------------------------------------------------

    Box xbx(bx); xbx.surroundingNodes(0);
    Box ybx(bx); ybx.surroundingNodes(1);

    FArrayBox fab_FX(xbx,ncomp,The_Async_Arena());
    FArrayBox fab_FE(ybx,ncomp,The_Async_Arena());

    auto FX=fab_FX.array();
    auto FE=fab_FE.array();

    ParallelFor(xbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        const Real pmon_u = (pm(i-1,j,0)+pm(i,j,0))/(pn(i-1,j,0)+pn(i,j,0));

        const Real cff = Real(0.25) * (diff2(i,j,0,n) + diff2(i-1,j,0,n)) * pmon_u;
        FX(i,j,k,n) = cff * (Hz(i,j,k) + Hz(i-1,j,k)) * (state_rhs(i,j,k,n)-state_rhs(i-1,j,k,n));
        FX(i,j,k,n) *= msku(i,j,0);
    });

    ParallelFor(ybx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        const Real pnom_v = (pn(i,j-1,0)+pn(i,j,0))/(pm(i,j-1,0)+pm(i,j,0));

        const Real cff = Real(0.25)*(diff2(i,j,0,n)+diff2(i,j-1,0,n)) * pnom_v;
        FE(i,j,k,n) = cff * (Hz(i,j,k) + Hz(i,j-1,k)) * (state_rhs(i,j,k,n) - state_rhs(i,j-1,k,n));
        FE(i,j,k,n) *= mskv(i,j,0);
    });

    /*
     Time-step harmonic, S-surfaces diffusion term.
    */
    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        const Real cff = dt_lev*pm(i,j,0)*pn(i,j,0);

        state(i,j,k,n) += cff * ( (FX(i+1,j  ,k,n)-FX(i,j,k,n))
                                 +(FE(i  ,j+1,k,n)-FE(i,j,k,n)) );
    });
}

/**
 * @param[in   ] bx         box to operate on
 * @param[inout] state      scalar data
 * @param[  out] state_rhs  scalar data RHS
 * @param[in   ] diff2      diffusivity
 * @param[in   ] Hz         vertical cell height
 * @param[in   ] z_r        rho-point z coordinates
 * @param[in   ] pm         1/dx
 * @param[in   ] pn         1/dy
 * @param[in   ] msku       land-sea mask on u-points
 * @param[in   ] mskv       land-sea mask on v-points
 * @param[in   ] dt_lev     time step at level
 * @param[in   ] ncomp      number of components to do this calculation on
 * @param[in   ] N          number of vertical levels
 */
void
REMORA::t3dmix2_geo(const Box& bx,
                   const Array4<Real      >& state,
                   const Array4<Real      >& state_rhs,
                   const Array4<Real const>& diff2,
                   const Array4<Real const>& Hz,
                   const Array4<Real const>& z_r,
                   const Array4<Real const>& pm,
                   const Array4<Real const>& pn,
                   const Array4<Real const>& msku,
                   const Array4<Real const>& mskv,
                   const Real dt_lev, const int ncomp, const int N)
{
    BL_PROFILE("REMORA::t3dmix2_geo()");
    //-----------------------------------------------------------------------
    //  Add in harmonic diffusivity s terms.
    //-----------------------------------------------------------------------

    Box xbx(bx); xbx.surroundingNodes(0);
    Box ybx(bx); ybx.surroundingNodes(1);
    Box zbx(bx); zbx.surroundingNodes(2);

    FArrayBox fab_dZdx(xbx,ncomp,The_Async_Arena());
    FArrayBox fab_dTdx(xbx,ncomp,The_Async_Arena());
    FArrayBox fab_dZde(ybx,ncomp,The_Async_Arena());
    FArrayBox fab_dTde(ybx,ncomp,The_Async_Arena());
    FArrayBox fab_dTdz(grow(zbx,IntVect(1,1,0)),ncomp,The_Async_Arena());
    FArrayBox fab_FS(grow(zbx,IntVect(1,1,0)),ncomp,The_Async_Arena());
    FArrayBox fab_FX(xbx,ncomp,The_Async_Arena());
    FArrayBox fab_FE(ybx,ncomp,The_Async_Arena());

    auto dZdx = fab_dZdx.array();
    auto dTdx = fab_dTdx.array();
    auto dZde = fab_dZde.array();
    auto dTde = fab_dTde.array();
    auto dTdz = fab_dTdz.array();
    auto FS   = fab_FS.array();
    auto FX   = fab_FX.array();
    auto FE   = fab_FE.array();

    ParallelFor(xbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real cff = half * (pm(i,j,0) + pm(i-1,j,0)) * msku(i,j,0);
        dZdx(i,j,k,n) = cff * (z_r(i,j,k) - z_r(i-1,j,k));
        dTdx(i,j,k,n)=cff*(state_rhs(i  ,j,k,n)-state_rhs(i-1,j,k,n));
    });
    ParallelFor(ybx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real cff = half * (pn(i,j,0) + pn(i,j-1,0)) * mskv(i,j,0);
        dZde(i,j,k,n) = cff * (z_r(i,j,k) - z_r(i,j-1,k));
        dTde(i,j,k,n)=cff*(state_rhs(i,j,k,n)-state_rhs(i,j-1,k,n));
    });
    ParallelFor(makeSlab(grow(bx,IntVect(1,1,0)),2,0), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int , int n)
    {
        dTdz(i,j,0,n) = zero;
        dTdz(i,j,N+1,n) = zero;
        FS(i,j,0,n) = zero;
        FS(i,j,N+1,n) = zero;
    });
    ParallelFor(grow(zbx,IntVect(1,1,-1)), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real cff = one / (z_r(i,j,k)-z_r(i,j,k-1));
        dTdz(i,j,k,n) = cff * (state_rhs(i,j,k,n) - state_rhs(i,j,k-1,n));
    });

    //  Compute components of the rotated tracer flux (T m3/s) along
    //  geopotential surfaces.
    ParallelFor(xbx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real on_u = two / (pn(i,j,0) + pn(i-1,j,0));
        Real cff = Real(0.25) * (diff2(i,j,0,n)+diff2(i-1,j,0,n)) * on_u;
        FX(i,j,k,n) = cff *
                       (Hz(i,j,k)+Hz(i-1,j,k))*
                       (dTdx(i,j,k,n)-
                        half*(std::min(dZdx(i,j,k,n),zero)*
                                   (dTdz(i-1,j,k  ,n)+
                                    dTdz(i  ,j,k+1,n))+
                                std::max(dZdx(i,j,k,n),zero)*
                                   (dTdz(i-1,j,k+1,n)+
                                    dTdz(i  ,j,k  ,n))));
    });
    ParallelFor(ybx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real om_v = two / (pm(i,j,0) + pm(i,j-1,0));
        Real cff = Real(0.25) * (diff2(i,j,0,n)+diff2(i,j-1,0,n)) * om_v;
        FE(i,j,k,n) = cff *
                       (Hz(i,j,k)+Hz(i,j-1,k))*
                       (dTde(i,j,k,n)-
                        half*(std::min(dZde(i,j,k,n),zero)*
                                   (dTdz(i,j-1,k  ,n)+
                                    dTdz(i,j  ,k+1,n))+
                                std::max(dZde(i,j,k,n),zero)*
                                   (dTdz(i,j-1,k+1,n)+
                                    dTdz(i,j  ,k  ,n))));
    });
    ParallelFor(grow(zbx,IntVect(0,0,-1)), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real cff = half * diff2(i,j,0,n);
        Real cff1=std::min(dZdx(i  ,j,k-1,n),zero);
        Real cff2=std::min(dZdx(i+1,j,k  ,n),zero);
        Real cff3=std::max(dZdx(i  ,j,k  ,n),zero);
        Real cff4=std::max(dZdx(i+1,j,k-1,n),zero);
        FS(i,j,k,n) = cff *
                        (cff1*(cff1*dTdz(i,j,k,n)-dTdx(i  ,j,k-1,n))+
                         cff2*(cff2*dTdz(i,j,k,n)-dTdx(i+1,j,k  ,n))+
                         cff3*(cff3*dTdz(i,j,k,n)-dTdx(i  ,j,k  ,n))+
                         cff4*(cff4*dTdz(i,j,k,n)-dTdx(i+1,j,k-1,n)));
        cff1=std::min(dZde(i,j  ,k-1,n),zero);
        cff2=std::min(dZde(i,j+1,k  ,n),zero);
        cff3=std::max(dZde(i,j  ,k  ,n),zero);
        cff4=std::max(dZde(i,j+1,k-1,n),zero);
        FS(i,j,k,n)=FS(i,j,k,n)+
                            cff*
                            (cff1*(cff1*dTdz(i,j,k,n)-dTde(i,j  ,k-1,n))+
                             cff2*(cff2*dTdz(i,j,k,n)-dTde(i,j+1,k  ,n))+
                             cff3*(cff3*dTdz(i,j,k,n)-dTde(i,j  ,k  ,n))+
                             cff4*(cff4*dTdz(i,j,k,n)-dTde(i,j+1,k-1,n)));
    });

    // Time-step harmonic, geopotential diffusion term (m Tunits).

    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
    {
        Real cff=dt_lev*pm(i,j,0)*pn(i,j,0);
        Real cff1=cff*(FX(i+1,j  ,k,n)-FX(i,j,k,n));
        Real cff2=cff*(FE(i  ,j+1,k,n)-FE(i,j,k,n));
        Real cff3=dt_lev*(FS(i,j,k+1,n)-FS(i,j,k,n));
        Real cff4=cff1+cff2+cff3;
        state(i,j,k,n)=state(i,j,k,n)+cff4;
    });
}
