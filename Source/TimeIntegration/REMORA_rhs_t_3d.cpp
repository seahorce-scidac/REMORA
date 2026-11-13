#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] lev           level to operate on
 * @param[in   ] bx            tilebox
 * @param[inout] t             tracer data
 * @param[in   ] sstore        scratch space for tracer calculations
 * @param[in   ] Huon          u-volume flux
 * @param[in   ] Hvom          v-volume flux
 * @param[in   ] Hz            vertical cell height
 * @param[in   ] pn            1/dx
 * @param[in   ] pm            1/dy
 * @param[in   ] W             vertical velocity
 * @param        FC            temporary
 * @param[in   ] mskr          land-sea mask on rho-points
 * @param[in   ] msku          land-sea mask on u-points
 * @param[in   ] mskv          land-sea mask on v-points
 * @param[in   ] river_pos     river positions
 * @param[in   ] river_source  river source data to add, if using
 * @param[in   ] nrhs          index of RHS component
 * @param[in   ] nnew          index of current time step
 * @param[in   ] N             number of vertical levels
 * @param[in   ] dt_lev        time step at this level
 */
void
REMORA::rhs_t_3d (int lev, const Box& bx,
                 const Array4<Real      >& t,
                 const Array4<Real const>& sstore,
                 const Array4<Real const>& Huon,
                 const Array4<Real const>& Hvom,
                 const Array4<Real const>& Hz,
                 const Array4<Real const>& pn,
                 const Array4<Real const>& pm,
                 const Array4<Real const>& W ,
                 const Array4<Real      >& FC,
                 const Array4<Real const>& mskr,
                 const Array4<Real const>& msku,
                 const Array4<Real const>& mskv,
                 const Array4<int  const>& river_pos,
                 const Array4<Real const>& river_source,
                 int nrhs, int nnew, int N, Real dt_lev,
                 const Array4<Real>& fx_reg,
                 const Array4<Real>& fy_reg)
{
    BL_PROFILE("REMORA::rhs_t_3d()");
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
    Box ubx = surroundingNodes(bx, 0);
    Box vbx = surroundingNodes(bx, 1);
    Box utbxp1 = surroundingNodes(grow(bx,0,1), 0);
    Box vtbxp1 = surroundingNodes(grow(bx,1,1), 1);

    //
    // Tile-private scratch: the gradients FX0 and FE0, and the fluxes FX and FE. Every
    // cell that is read is written first. The river pass writes on the tile grown by one.
    //
    int ncomp = 0;
    const int FX0_comp = ncomp++;
    const int FE0_comp = ncomp++;
    const int FX_comp  = ncomp++;
    const int FE_comp  = ncomp++;
    FArrayBox fab(grow(bx,IntVect(NGROW,NGROW,0)),ncomp,amrex::The_Async_Arena());
    auto FX0 = fab.array(FX0_comp);
    auto FE0 = fab.array(FE0_comp);
    auto FX  = fab.array(FX_comp);
    auto FE  = fab.array(FE_comp);

    BL_PROFILE_VAR("REMORA::rhs_t_3d()::hadv",phadv);
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
        FX0(i,j,k)=(sstore(ii,j,k,nrhs)-sstore(ii-1,j,k,nrhs)) * msku(ii,j,0);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        const int jj = amrex::min(amrex::max(j,jlo),jhi);
        FE0(i,j,k)=(sstore(i,jj,k,nrhs)-sstore(i,jj-1,k,nrhs)) * mskv(i,jj,0);
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
            FX(i,j,k)=Huon(i,j,k)*Real(0.5)*(sstore(i,j,k)+sstore(i-1,j,k))-
                cffa*(curv_i*min_Huon+ curv_im1*max_Huon);
        } else {
            const Real grad_i   = Real(0.5)*(FX0(i  ,j,k)+FX0(i+1,j,k));
            const Real grad_im1 = Real(0.5)*(FX0(i-1,j,k)+FX0(i  ,j,k));
            FX(i,j,k)=Huon(i,j,k)*Real(0.5)*(sstore(i,j,k)+sstore(i-1,j,k)-
                                       cffb*(grad_i- grad_im1));
        }
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k, auto up3)
    {
        if (up3.value) {
            const Real curv_j   = -FE0(i,j  ,k)+FE0(i,j+1,k);
            const Real curv_jm1 = -FE0(i,j-1,k)+FE0(i,j  ,k);

            Real max_Hvom = std::max(Hvom(i,j,k),Real(0.0));
            Real min_Hvom = std::min(Hvom(i,j,k),Real(0.0));

            FE(i,j,k)=Hvom(i,j,k)*Real(0.5)*(sstore(i,j,k)+sstore(i,j-1,k))-
                cffa*(curv_j*min_Hvom+ curv_jm1*max_Hvom);
        } else {
            const Real grad_j   = Real(0.5)*(FE0(i,j  ,k)+FE0(i,j+1,k));
            const Real grad_jm1 = Real(0.5)*(FE0(i,j-1,k)+FE0(i,j  ,k));
            FE(i,j,k)=Hvom(i,j,k)*Real(0.5)*(sstore(i,j,k)+sstore(i,j-1,k)-
                                       cffb*(grad_j- grad_jm1));
        }
    });

    bool do_rivers_cons = (river_source.size() > 0);
    if (solverChoice.do_rivers) {
        int* river_direction_d = river_direction.data();
        ParallelFor(grow(bx,IntVect(NGROW-1,NGROW-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            int iriver = river_pos(i,j,0);
            if (iriver >= 0) {
                if (river_direction_d[iriver] == 0) {
                    if (do_rivers_cons) {
                        FX(i,j,k) = Huon(i,j,k) * river_source(iriver,0,k);
                    } else if ((mskr(i,j,0)==0) && (mskr(i-1,j,0)==1)) {
                        FX(i,j,k) = Huon(i,j,k) * sstore(i-1,j,k);
                    } else if ((mskr(i,j,0)==1) && (mskr(i-1,j,0)==0)) {
                        FX(i,j,k) = Huon(i,j,k) * sstore(i,j,k);
                    }
                } else if (river_direction_d[iriver] == 1) {
                    if (do_rivers_cons) {
                        FE(i,j,k) = Hvom(i,j,k) * river_source(iriver,0,k);
                    } else if ((mskr(i,j,0)==0) && (mskr(i,j-1,0)==1)) {
                        FE(i,j,k) = Hvom(i,j,k) * sstore(i,j-1,k);
                    } else if ((mskr(i,j,0)==1) && (mskr(i,j-1,0)==0)) {
                        FE(i,j,k) = Hvom(i,j,k) * sstore(i,j,k);
                    }
                }
            }
        });
    }

    // Hand the same fluxes to the flux register, scaled so its internal dt/dx reproduces
    // the dt*pm*pn above. REMORA's divergence uses the receiving cell's metric, not a face
    // one, so the face weight is an average: exact where pm and pn are uniform.
    if (fx_reg) {
        const Real dx0 = geom[lev].CellSize(0);
        const Real dx1 = geom[lev].CellSize(1);
        ParallelFor(surroundingNodes(bx,0), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real w = Real(0.5) * (pm(i-1,j,0)*pn(i-1,j,0) + pm(i,j,0)*pn(i,j,0));
            fx_reg(i,j,k) = FX(i,j,k) * w * dx0;
        });
        ParallelFor(surroundingNodes(bx,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real w = Real(0.5) * (pm(i,j-1,0)*pn(i,j-1,0) + pm(i,j,0)*pn(i,j,0));
            fy_reg(i,j,k) = FE(i,j,k) * w * dx1;
        });
    }

    BL_PROFILE_VAR_STOP(phadv);
    //-----------------------------------------------------------------------
    //  Time-step vertical advection term.
    //-----------------------------------------------------------------------
    //Check which type of differences:
    //
    //  Fourth-order, central differences vertical advective flux
    //  (Tunits m3/s).
    //

    BL_PROFILE_VAR("REMORA::rhs_t_3d()::vadv",pvadv);
    ParallelFor(surroundingNodes(bx,2), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //-----------------------------------------------------------------------
        //  Add in vertical advection.
        //-----------------------------------------------------------------------

        Real cff1=Real(0.5);
        Real cff2=Real(7.0)/Real(12.0);
        Real cff3=one/Real(12.0);

        if (k>=2 && k<=N-1)
        {
            FC(i,j,k)=( cff2*(sstore(i  ,j,k-1)+ sstore(i,j,k))
                        -cff3*(sstore(i  ,j,k-2)+ sstore(i,j,k+1)) ) * ( W(i,j,k));
        } else if (k==N+1) {
            FC(i,j,N+1)=Real(0.0);
        } else if (k==N) {
            FC(i,j,N)=( cff2*sstore(i  ,j,N-1)+ cff1*sstore(i,j,N  )
                         -cff3*sstore(i  ,j,N-2) ) * ( W(i  ,j,N));
        } else if (k==1) {
            FC(i,j,1)=( cff2*sstore(i  ,j,1)+ cff1*sstore(i,j,0)
                       -cff3*sstore(i  ,j,2) ) * ( W(i  ,j,1));
        } else if (k==0) {
            FC(i,j,0) = Real(0.0);
        }
    });

    // Horizontal and then vertical advection of the tracer, in one pass
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        //
        //  Add in horizontal advection.
        //
        Real cff = dt_lev*pm(i,j,0)*pn(i,j,0);
        Real cff1=cff*(FX(i+1,j,k)-FX(i,j,k));
        Real cff2=cff*(FE(i,j+1,k)-FE(i,j,k));
        Real cff3=cff1+cff2;

        const Real t_h = t(i,j,k,nnew) - cff3;

        Real cff4=FC(i,j,k+1)-FC(i,j,k);

        t(i,j,k,nnew) = (t_h-cff*cff4) / Hz(i,j,k);
    });
    BL_PROFILE_VAR_STOP(pvadv);
}
