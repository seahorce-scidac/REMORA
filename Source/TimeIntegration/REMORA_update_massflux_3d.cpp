#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] lev        level of refinement (coarsest level is 0)
 * @param[in   ] bx             box on which to update
 * @param[in   ] ioff           offset in x-direction
 * @param[in   ] joff           offset in y-direction
 * @param[inout] phi            u or v
 * @param[inout] phibar         ubar or vbar
 * @param[inout] Hphi           phi-volume flux
 * @param[in   ] Hz             vertical cell height
 * @param[in   ] om_v_or_on_u   on_u for u, om_v for v
 * @param[in   ] Dphi_avg1      DU_avg1 or DV_avg1
 * @param[in   ] Dphi_avg2      DU_avg2 or DV_avg2
 * @param        DC             temporary
 * @param[in   ] msk            land-sea mask
 * @param[in   ] nnew           component of velocity
 */
void
REMORA::update_massflux_3d (int lev, const Box& bx,
                           const int ioff, const int joff,
                           const Array4<Real      >& phi,
                           const Array4<Real      >& phibar,
                           const Array4<Real      >& Hphi,
                           const Array4<Real const>& Hz,
                           const Array4<Real const>& om_v_or_on_u,
                           const Array4<Real const>& Dphi_avg1,
                           const Array4<Real const>& Dphi_avg2,
                           const Array4<Real      >& DC,
                           const Array4<Real const>& msk,
                           const int nnew)
{
    BL_PROFILE("REMORA::update_massflux_3d()");
    const Box& domain = geom[lev].Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    GeometryData const& geomdata = geom[0].data();
    bool is_periodic_in_x = geomdata.isPeriodic(0);
    bool is_periodic_in_y = geomdata.isPeriodic(1);

    auto N = Geom(0).Domain().size()[2]-1; // Number of vertical "levs" aka, NZ
    AMREX_ASSERT(bx.smallEnd(2) == 0 && bx.bigEnd(2) == N);

    auto bxD=bx;
    bxD.makeSlab(2,0);

    // Tile-private column sums, each built up level by level from the bottom: the column's
    // face area (then its inverse), the transport sum of DC*phi (then its correction), and
    // the sum of Hphi (then its correction)
    FArrayBox fab_sum(bxD,3,amrex::The_Async_Arena());
    auto Dsum = fab_sum.array(0);
    auto CF   = fab_sum.array(1);
    auto FC   = fab_sum.array(2);

    //Copied depth of water column calculation from DepthStretchTransform
    // Face areas of the u- or v-cells, DC(i,j,0:N) = on_u (or om_v) times the face thickness

    //This takes advantage of Hz being an extra grow cell size
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        DC(i,j,k) = Real(0.5) * om_v_or_on_u(i,j,0) * (Hz(i,j,k)+Hz(i-ioff,j-joff,k));
    });

    ParallelFor(bxD, N+1, [=] AMREX_GPU_DEVICE (int i, int j, int, int k)
    {
        Dsum(i,j,0) = ((k == 0) ? zero : Dsum(i,j,0)) + DC(i,j,k);
        CF(i,j,0)   = ((k == 0) ? zero : CF(i,j,0)) + DC(i,j,k) * phi(i,j,k,nnew);

        if (k == N) {
            Dsum(i,j,0) = one / Dsum(i,j,0);
            CF(i,j,0)  = Dsum(i,j,0) * (CF(i,j,0) - Dphi_avg1(i,j,0));

            phibar(i,j,0,0) = Dsum(i,j,0) * Dphi_avg1(i,j,0);
            phibar(i,j,0,1) = phibar(i,j,0,0);
        }
    });

    // Correct phi on the non-periodic domain edges, then average the volume flux with the
    // one from phi and sum it over the column
    ParallelFor(bxD, N+1, [=] AMREX_GPU_DEVICE (int i, int j, int, int k)
    {
        const bool edge = (i == dlo.x-joff && !is_periodic_in_x) ||
                          (i == dhi.x+1    && !is_periodic_in_x) ||
                          (j == dlo.y-ioff && !is_periodic_in_y) ||
                          (j == dhi.y+1    && !is_periodic_in_y);
        if (edge) {
            phi(i,j,k,nnew) -= CF(i,j,0);
            phi(i,j,k,nnew) *= msk(i,j,0);
        }

        Hphi(i,j,k) = Real(0.5) * (Hphi(i,j,k)+phi(i,j,k,nnew)*DC(i,j,k));
        FC(i,j,0) = ((k == 0) ? zero : FC(i,j,0)) + Hphi(i,j,k);

        if (k == N) {
            FC(i,j,0) = Dsum(i,j,0)*(FC(i,j,0)-Dphi_avg2(i,j,0));
        }
    });

    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Hphi(i,j,k) -= DC(i,j,k)*FC(i,j,0);
    });
}
