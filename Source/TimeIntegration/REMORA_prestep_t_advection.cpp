#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] lev             level to operate on
 * @param[in   ] tbx             tile box
 * @param[in   ] gbx             grown tile box
 * @param[in   ] tempold         scalar at last time
 * @param[in   ] tempcache       cached current time step's scalar value
 * @param[in   ] Hz              vertical cell height
 * @param[in   ] Huon            u-volume flux
 * @param[in   ] Hvom            v-volume flux
 * @param[in   ] Akv             vertical viscosity coefficient
 * @param[inout] W               vertical velocity
 * @param        DC              temporary
 * @param        FC              temporary
 * @param[  out] tempstore       scratch space for calculations on scalars
 * @param[in   ] z_w             z coordinates at w points
 * @param[in   ] h               bathymetry
 * @param[in   ] pm              1/dx
 * @param[in   ] pn              1/dy
 * @param[in   ] msku            land-sea mask on u-points
 * @param[in   ] mskv            land-sea mask on v-points
 * @param[in   ] river_pos       river positions
 * @param[in   ] river_source    river source data to add, if using
 * @param[in   ] iic             which time step we're on
 * @param[in   ] ntfirst         what is the first time step?
 * @param[in   ] nrhs            index of RHS component
 * @param[in   ] N               number of vertical levels
 * @param[in   ] dt_lev          time step at this level
 */

void
REMORA::prestep_t_advection (int lev, const Box& tbx, const Box& gbx,
                            const Array4<Real      >& tempold,
                            const Array4<Real      >& tempcache,
                            const Array4<Real      >& Hz,
                            const Array4<Real      >& Huon,
                            const Array4<Real      >& Hvom,
                            const Array4<Real      >& W  ,
                            const Array4<Real      >& DC,
                            const Array4<Real      >& FC ,
                            const Array4<Real      >& tempstore,
                            const Array4<Real const>& z_w,
                            const Array4<Real const>& h,
                            const Array4<Real const>& pm,
                            const Array4<Real const>& pn,
                            const Array4<Real const>& msku,
                            const Array4<Real const>& mskv,
                            const Array4<int  const>& river_pos,
                            const Array4<Real const>& river_source,
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

    //copy the tilebox
    Box gbx1 = tbx;
    Box gbx2 = tbx;
    Box tbxp1 = tbx;
    Box tbxp2 = tbx;

    tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
    tbxp2.grow(IntVect(NGROW,NGROW,0));

    int ncomp = 0;
    int FX_comp = ncomp++;
    int FE_comp = ncomp++;
    int curv_comp = ncomp++;
    int grad_comp = ncomp++;

    FArrayBox fab(tbxp2,ncomp,amrex::The_Async_Arena());
    fab.template setVal<RunOn::Device>(0.0_rt);

    auto FX=fab.array(FX_comp);
    auto FE=fab.array(FE_comp);
    auto curv=fab.array(curv_comp);
    auto grad=fab.array(grad_comp);

    Box ubx = surroundingNodes(tbx,0);
    Box vbx = surroundingNodes(tbx,1);

    Box utbxp1 = surroundingNodes(tbxp1,0);
    Box vtbxp1 = surroundingNodes(tbxp1,1);

    Box gbx3uneven_init(IntVect(AMREX_D_DECL(tbx.smallEnd(0)-3,tbx.smallEnd(1)-3,tbx.smallEnd(2))),
                   IntVect(AMREX_D_DECL(tbx.bigEnd(0)+2,tbx.bigEnd(1)+2,tbx.bigEnd(2))));
    BoxArray ba_gbx3uneven = intersect(BoxArray(gbx3uneven_init), gbx);
    AMREX_ASSERT((ba_gbx3uneven.size() == 1));
    Box gbx3uneven = ba_gbx3uneven[0];

    gbx2.grow(IntVect(NGROW,NGROW,0));
    BoxArray ba_gbx2 = intersect(BoxArray(gbx2), gbx);
    AMREX_ASSERT((ba_gbx2.size() == 1));
    gbx2 = ba_gbx2[0];

    gbx1.grow(IntVect(NGROW-1,NGROW-1,0));
    BoxArray ba_gbx1 = intersect(BoxArray(gbx1), gbx);
    AMREX_ASSERT((ba_gbx1.size() == 1));
    gbx1 = ba_gbx1[0];

    //------------------------------------------------------------------------
    //  Vertically integrate horizontal mass flux divergence.
    //------------------------------------------------------------------------
    //
    //Should really use gbx3uneven
    Box gbx3unevenD = gbx3uneven;
    gbx3unevenD.makeSlab(2,0);
    Box gbx1D = gbx1;
    gbx1D.makeSlab(2,0);

    // We used to set W(i,j,0) = 0.0, but this should have already been set before passing into the function
    ParallelFor(gbx1D, N+1,
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
    ParallelFor(gbx1D, [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        W(i,j,N+1)=W(i,j,N+1)/(z_w(i,j,N+1)+h(i,j,0,0)); // wrk_i
    });
    ParallelFor(gbx1D, N,
    [=] AMREX_GPU_DEVICE (int i, int j, int , int kk)
    {
        int k = kk + 1;
        W(i,j,k) = W(i,j,k)- W(i,j,N+1)*(z_w(i,j,k)+h(i,j,0,0));
    });
    ParallelFor(gbx1D, [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        W(i,j,N+1) = 0.0_rt;
    });

    //Use FC and DC as intermediate arrays for FX and FE
    //First pass do centered 2d terms
    if (solverChoice.flat_bathymetry) {
        ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            FX(i,j,k)=Box(tempold).contains(i-1,j,k) ? Huon(i,j,k)*
                        0.5_rt*(tempold(i-1,j,k)+tempold(i  ,j,k)) : 1e34_rt;
        });
        ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            FE(i,j,k)=Box(tempold).contains(i,j-1,k) ? Hvom(i,j,k)*
                        0.5_rt*(tempold(i,j-1,k)+tempold(i,j,k)) : 1e34_rt;
        });

    } else {

        ParallelFor(utbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            //should be t index 3
            FX(i,j,k)=(tempold(i,j,k,nrhs)-tempold(i-1,j,k,nrhs)) * msku(i,j,0);
        });

        Box utbxp1_slab_lo = makeSlab(utbxp1,0,dlo.x-1) & utbxp1;
        Box utbxp1_slab_hi = makeSlab(utbxp1,0,dhi.x+1) & utbxp1;
        if (utbxp1_slab_lo.ok() && !is_periodic_in_x) {
            ParallelFor(utbxp1_slab_lo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                FX(i,j,k) = FX(i+1,j,k);
            });
        }
        if (utbxp1_slab_hi.ok() && !is_periodic_in_x) {
            ParallelFor(utbxp1_slab_hi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                FX(i+1,j,k) = FX(i,j,k);
            });
        }

        ParallelFor(vtbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            //should be t index 3
            FE(i,j,k)=(tempold(i,j,k,nrhs)-tempold(i,j-1,k,nrhs)) * mskv(i,j,0);
        });

        Box vtbxp1_slab_lo = makeSlab(vtbxp1,1,dlo.y-1) & vtbxp1;
        Box vtbxp1_slab_hi = makeSlab(vtbxp1,1,dhi.y+1) & vtbxp1;
        if (vtbxp1_slab_lo.ok() && !is_periodic_in_y) {
            ParallelFor(vtbxp1_slab_lo, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                FE(i,j,k) = FE(i,j+1,k);
            });
        }
        if (vtbxp1_slab_hi.ok() && !is_periodic_in_y) {
            ParallelFor(vtbxp1_slab_hi, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                FE(i,j+1,k) = FE(i,j,k);
            });
        }

        Real cffa=1.0_rt/6.0_rt;
        Real cffb=1.0_rt/3.0_rt;
        if (solverChoice.tracer_Hadv_scheme == AdvectionScheme::upstream3)
        {
            ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                //Upstream3
                curv(i,j,k)=-FX(i,j,k)+FX(i+1,j,k);
            });

            ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real max_Huon = std::max(Huon(i,j,k),0.0_rt);
                Real min_Huon = std::min(Huon(i,j,k),0.0_rt);

                FX(i,j,k)=Huon(i,j,k)*0.5_rt*(tempold(i,j,k)+tempold(i-1,j,k))-
                    cffa*(curv(i,j,k)*min_Huon+ curv(i-1,j,k)*max_Huon);
            });

        } else if (solverChoice.tracer_Hadv_scheme == AdvectionScheme::centered4) {

            ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                //Centered4
                grad(i,j,k)=0.5_rt*(FX(i,j,k)+FX(i+1,j,k));
            });

            ParallelFor(ubx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                FX(i,j,k)=Huon(i,j,k)*0.5_rt*(tempold(i,j,k)+tempold(i-1,j,k)-
                                           cffb*(grad(i,j,k)-grad(i-1,j,k)));
            });
        } else {
           Error("Not a valid horizontal advection scheme");
        }

        if (solverChoice.tracer_Hadv_scheme == AdvectionScheme::upstream3)
        {
            ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                curv(i,j,k)=-FE(i,j,k)+FE(i,j+1,k);
            });

            ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real max_Hvom = std::max(Hvom(i,j,k),0.0_rt);
                Real min_Hvom = std::min(Hvom(i,j,k),0.0_rt);

                FE(i,j,k)=Hvom(i,j,k)*0.5_rt*(tempold(i,j,k)+tempold(i,j-1,k))-
                    cffa*(curv(i,j,k)*min_Hvom+ curv(i,j-1,k)*max_Hvom);
            });

        } else if (solverChoice.tracer_Hadv_scheme == AdvectionScheme::centered4) {

            ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                grad(i,j,k)=0.5_rt*(FE(i,j,k)+FE(i,j+1,k));
            });

            ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                FE(i,j,k)=Hvom(i,j,k)*0.5_rt*(tempold(i,j,k)+tempold(i,j-1,k)-
                                           cffb*(grad(i,j,k)- grad(i,j-1,k)));
            });
        } else {
            Error("Not a valid horizontal advection scheme");
        }
    } // not flat

    bool do_rivers_cons = (river_source.size() > 0);
    if (solverChoice.do_rivers) {
        int* river_direction_d = river_direction.data();
        ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            int iriver = river_pos(i,j,0);
            if (iriver >= 0) {
                if (river_direction_d[iriver] == 0) {
                    FX(i,j,k) = (!do_rivers_cons) ? 0.0_rt : Huon(i,j,k) * river_source(iriver,0,k);
                } else {
                    FE(i,j,k) = (!do_rivers_cons) ? 0.0_rt : Hvom(i,j,k) * river_source(iriver,0,k);
                }
            }
        });
    }

    //Intermediate tracer at 3
    //
    //  Time-step horizontal advection (m Tunits).
    //
    Real cff1 = 0.0_rt, cff2 = 0.0_rt, cff;

    Real GammaT = 1.0_rt/6.0_rt;

    if (iic==ntfirst)
    {
        cff=0.5_rt*dt_lev;
        cff1=1.0_rt;
        cff2=0.0_rt;
    } else {
        cff=(1.0_rt-GammaT)*dt_lev;
        cff1=0.5_rt+GammaT;
        cff2=0.5_rt-GammaT;
    }

    ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        tempstore(i,j,k)=Hz(i,j,k)*( cff1 * tempold(i,j,k)+
                                     cff2 * tempcache(i,j,k) )-
                                     cff  * pm(i,j,0)*pn(i,j,0) * (FX(i+1,j,k)-FX(i,j,k)+
                                                                   FE(i,j+1,k)-FE(i,j,k));
    });

    Real c1=0.5_rt;
    Real c2=7.0_rt/12.0_rt;
    Real c3=1.0_rt/12.0_rt;

    //
    // Time-step vertical advection of tracers (Tunits). Impose artificial
    // continuity equation.
    //
    ParallelFor(growHi(growLo(convert(tbx,IntVect(0,0,1)),2,-2),2,-1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //-----------------------------------------------------------------------
        //  Add in vertical advection.
        //-----------------------------------------------------------------------

        FC(i,j,k)=( c2*(tempold(i  ,j,k-1,nrhs)+ tempold(i,j,k  ,nrhs))
                         -c3*(tempold(i  ,j,k-2 ,nrhs)+ tempold(i,j,k+1,nrhs)) )*
                            ( W(i,j,k));
    });
    ParallelFor(makeSlab(tbx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int )
    {
        FC(i,j,N+1)=0.0_rt;
        FC(i,j,N) = ( c2*tempold(i,j,N-1,nrhs)+ c1*tempold(i,j,N,nrhs)-c3*tempold(i,j,N-2,nrhs) )
                  * W(i,j,N);
        FC(i,j,  1) = ( c2*tempold(i,j,  1,nrhs)+ c1*tempold(i,j,0,nrhs)-c3*tempold(i,j,2,nrhs) )
                    * W(i,j,1);
        FC(i,j,  0) = 0.0_rt;
    });

    ParallelFor(tbxp1, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        DC(i,j,k)=1.0_rt/(Hz(i,j,k)-
                        cff*pm(i,j,0)*pn(i,j,0)*
                        (Huon(i+1,j,k)-Huon(i,j,k)+
                         Hvom(i,j+1,k)-Hvom(i,j,k)+
                         (W(i,j,k+1)-W(i,j,k))));
    });

    ParallelFor(tbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real c_p = cff*pm(i,j,0)*pn(i,j,0);

        Real c4 = FC(i,j,k+1)-FC(i,j,k);

        tempstore(i,j,k) = DC(i,j,k)*(tempstore(i,j,k)-c_p*c4);
    });
}
