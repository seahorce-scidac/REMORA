#include <REMORA_BottomExchange.H>
#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] lev             level to operate on
 * @param[in   ] tbx             tile box
 * @param[in   ] tempold         scalar at last time
 * @param[in   ] tempcache       current time step's scalar value, read on tbx only
 * @param[in   ] Hz              vertical cell height
 * @param[in   ] Huon            u-volume flux
 * @param[in   ] Hvom            v-volume flux
 * @param[in   ] W               vertical velocity, on the columns of tbx
 * @param        FC              temporary
 * @param[  out] tempstore       scratch space for calculations on scalars
 * @param[in   ] pm              1/dx
 * @param[in   ] pn              1/dy
 * @param[in   ] msku            land-sea mask on u-points
 * @param[in   ] mskv            land-sea mask on v-points
 * @param[in   ] river_pos       river positions
 * @param[in   ] river_source    river source data to add, if using
 * @param[in   ] bottom_incoming incoming bottom tracer transport
 * @param[in   ] iic             which time step we're on
 * @param[in   ] ntfirst         what is the first time step?
 * @param[in   ] nrhs            index of RHS component
 * @param[in   ] N               number of vertical levels
 * @param[in   ] dt_lev          time step at this level
 */

void
REMORA::prestep_t_advection (int lev, const Box& tbx,
                            const Array4<Real      >& tempold,
                            const Array4<Real      >& tempcache,
                            const Array4<Real      >& Hz,
                            const Array4<Real      >& Huon,
                            const Array4<Real      >& Hvom,
                            const Array4<Real const>& W  ,
                            const Array4<Real      >& FC ,
                            const Array4<Real      >& tempstore,
                            const Array4<Real const>& pm,
                            const Array4<Real const>& pn,
                            const Array4<Real const>& msku,
                            const Array4<Real const>& mskv,
                            const Array4<int  const>& river_pos,
                            const Array4<Real const>& river_source,
                            const Array4<Real const>& bottom_incoming,
                            int iic, int ntfirst, int nrhs, int N,
                            Real dt_lev)
{
    BL_PROFILE("REMORA::prestep_t_advection()");
    const Box& domain = geom[lev].Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    GeometryData const& geomdata = geom[0].data();
    bool is_periodic_in_x = geomdata.isPeriodic(0);
    bool is_periodic_in_y = geomdata.isPeriodic(1);

    const auto hadv = solverChoice.tracer_Hadv_scheme;
    if (hadv != AdvectionScheme::upstream3 && hadv != AdvectionScheme::centered4) {
        Error("Not a valid horizontal advection scheme");
    }
    const int upstream3 = (hadv == AdvectionScheme::upstream3) ? 1 : 0;

    // The fluxes are needed on the faces of the tile, ubx and vbx; the gradients they
    // are built from one face further out, on utbxp1 and vtbxp1.
    Box ubx = surroundingNodes(tbx,0);
    Box vbx = surroundingNodes(tbx,1);
    Box utbxp1 = surroundingNodes(grow(tbx,0,1),0);
    Box vtbxp1 = surroundingNodes(grow(tbx,1,1),1);

    // Tile-private scratch: the gradients FX0 and FE0, and the fluxes FX and FE. Every
    // cell that is read is written first. The river pass writes on the tile grown by one.
    int ncomp = 0;
    const int FX0_comp = ncomp++;
    const int FE0_comp = ncomp++;
    const int FX_comp  = ncomp++;
    const int FE_comp  = ncomp++;
    FArrayBox fab(grow(tbx,IntVect(NGROW,NGROW,0)),ncomp,amrex::The_Async_Arena());
    auto FX0 = fab.array(FX0_comp);
    auto FE0 = fab.array(FE0_comp);
    auto FX  = fab.array(FX_comp);
    auto FE  = fab.array(FE_comp);

    // Gradients across the faces. Outside a non-periodic domain edge the gradient is the
    // one on the edge face, so the face index is clamped to the domain's faces there.
    const int ilo = is_periodic_in_x ? std::numeric_limits<int>::lowest() : dlo.x;
    const int ihi = is_periodic_in_x ? std::numeric_limits<int>::max()    : dhi.x+1;
    const int jlo = is_periodic_in_y ? std::numeric_limits<int>::lowest() : dlo.y;
    const int jhi = is_periodic_in_y ? std::numeric_limits<int>::max()    : dhi.y+1;
    ParallelFor(utbxp1, vtbxp1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const int ii = amrex::min(amrex::max(i,ilo),ihi);
        //should be t index 3
        FX0(i,j,k)=(tempold(ii,j,k,nrhs)-tempold(ii-1,j,k,nrhs)) * msku(ii,j,0);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const int jj = amrex::min(amrex::max(j,jlo),jhi);
        //should be t index 3
        FE0(i,j,k)=(tempold(i,jj,k,nrhs)-tempold(i,jj-1,k,nrhs)) * mskv(i,jj,0);
    });

    // Fluxes on the faces of the tile, from the curvature (upstream3) or the averaged
    // gradient (centered4) at the cells on either side of the face
    Real cffa=one/Real(6.0);
    Real cffb=one/Real(3.0);
    AnyCTO(TypeList<CompileTimeOptions<0,1>>{}, {upstream3},
    [&] (auto const& fx, auto const& fy) { ParallelFor(ubx, vbx, fx, fy); },
    [=] AMREX_GPU_DEVICE (int i, int j, int k, auto up3)
    {
        if (up3.value) {
            const Real curv_i   = -FX0(i  ,j,k)+FX0(i+1,j,k);
            const Real curv_im1 = -FX0(i-1,j,k)+FX0(i  ,j,k);

            Real max_Huon = std::max(Huon(i,j,k),Real(0.0));
            Real min_Huon = std::min(Huon(i,j,k),Real(0.0));

            FX(i,j,k)=Huon(i,j,k)*Real(0.5)*(tempold(i,j,k)+tempold(i-1,j,k))-
                cffa*(curv_i*min_Huon+ curv_im1*max_Huon);
        } else {
            const Real grad_i   = Real(0.5)*(FX0(i  ,j,k)+FX0(i+1,j,k));
            const Real grad_im1 = Real(0.5)*(FX0(i-1  ,j,k)+FX0(i  ,j,k));
            FX(i,j,k)=Huon(i,j,k)*Real(0.5)*(tempold(i,j,k)+tempold(i-1,j,k)-
                                       cffb*(grad_i-grad_im1));
        }
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k, auto up3)
    {
        if (up3.value) {
            const Real curv_j   = -FE0(i,j  ,k)+FE0(i,j+1,k);
            const Real curv_jm1 = -FE0(i,j-1,k)+FE0(i,j  ,k);

            Real max_Hvom = std::max(Hvom(i,j,k),Real(0.0));
            Real min_Hvom = std::min(Hvom(i,j,k),Real(0.0));

            FE(i,j,k)=Hvom(i,j,k)*Real(0.5)*(tempold(i,j,k)+tempold(i,j-1,k))-
                cffa*(curv_j*min_Hvom+ curv_jm1*max_Hvom);
        } else {
            const Real grad_j   = Real(0.5)*(FE0(i,j  ,k)+FE0(i,j+1,k));
            const Real grad_jm1 = Real(0.5)*(FE0(i,j-1,k)+FE0(i,j  ,k));
            FE(i,j,k)=Hvom(i,j,k)*Real(0.5)*(tempold(i,j,k)+tempold(i,j-1,k)-
                                       cffb*(grad_j- grad_jm1));
        }
    });

    bool do_rivers_cons = (river_source.size() > 0);
    if (solverChoice.do_rivers) {
        int* river_direction_d = river_direction.data();
        ParallelFor(grow(tbx,IntVect(NGROW-1,NGROW-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            int iriver = river_pos(i,j,0);
            if (iriver >= 0) {
                if (river_direction_d[iriver] == 0) {
                    FX(i,j,k) = (!do_rivers_cons) ? Real(0.0) : Huon(i,j,k) * river_source(iriver,0,k);
                } else {
                    FE(i,j,k) = (!do_rivers_cons) ? Real(0.0) : Hvom(i,j,k) * river_source(iriver,0,k);
                }
            }
        });
    }

    Real c1=Real(0.5);
    Real c2=Real(7.0)/Real(12.0);
    Real c3=one/Real(12.0);

    //
    // Vertical advective flux at the w-levels of the tile's columns: fourth-order in the
    // interior, one-sided next to the bottom and the surface, and prescribed at the bottom.
    //
    ParallelFor(surroundingNodes(tbx,2), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if (k == 0) {
            FC(i,j,0) = bottom_incoming ? BottomExchange::tracer_flux(
                W(i,j,0), bottom_incoming(i,j,0), tempold(i,j,0,nrhs)) : Real(0.0);
        } else if (k == N+1) {
            FC(i,j,k) = Real(0.0);
        } else if (k >= 2 && k <= N-1) {
            FC(i,j,k)=( c2*(tempold(i  ,j,k-1,nrhs)+ tempold(i,j,k  ,nrhs))
                             -c3*(tempold(i  ,j,k-2 ,nrhs)+ tempold(i,j,k+1,nrhs)) )*
                                ( W(i,j,k));
        } else if (k == 1) {
            FC(i,j,  1) = ( c2*tempold(i,j,  1,nrhs)+ c1*tempold(i,j,0,nrhs)-c3*tempold(i,j,2,nrhs) )
                        * W(i,j,1);
        } else {
            FC(i,j,N) = ( c2*tempold(i,j,N-1,nrhs)+ c1*tempold(i,j,N,nrhs)-c3*tempold(i,j,N-2,nrhs) )
                      * W(i,j,N);
        }
    });

    //Intermediate tracer at 3
    //
    //  Time-step horizontal advection (m Tunits), then vertical advection of tracers
    //  (Tunits), imposing the artificial continuity equation through DC.
    //
    Real cff1 = zero, cff2 = zero, cff;

    Real GammaT = one/Real(6.0);

    if (iic==ntfirst)
    {
        cff=Real(0.5)*dt_lev;
        cff1=one;
        cff2=zero;
    } else {
        cff=(one-GammaT)*dt_lev;
        cff1=Real(0.5)+GammaT;
        cff2=Real(0.5)-GammaT;
    }

    ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const Real ts = Hz(i,j,k)*( cff1 * tempold(i,j,k)+
                                    cff2 * tempcache(i,j,k) )-
                                    cff  * pm(i,j,0)*pn(i,j,0) * (FX(i+1,j,k)-FX(i,j,k)+
                                                                  FE(i,j+1,k)-FE(i,j,k));

        const Real DC=one/(Hz(i,j,k)-
                        cff*pm(i,j,0)*pn(i,j,0)*
                        (Huon(i+1,j,k)-Huon(i,j,k)+
                         Hvom(i,j+1,k)-Hvom(i,j,k)+
                         (W(i,j,k+1)-W(i,j,k))));

        Real c_p = cff*pm(i,j,0)*pn(i,j,0);

        Real c4 = FC(i,j,k+1)-FC(i,j,k);

        tempstore(i,j,k) = DC*(ts-c_p*c4);
    });
}
