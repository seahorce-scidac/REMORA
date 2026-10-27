#include <REMORA.H>
#include <AMReX_OpenMP.H>

using namespace amrex;

/** Corresponds to step3d_uv.F and step3d_t.F in ROMS
 *
 * @param[in   ] lev        level of refinement (coarsest level is 0)
 * @param[inout] mf_cons    scalar variables (temp, salt, scalar, etc)
 * @param[inout] mf_u       u-velocity
 * @param[inout] mf_v       v-velocity
 * @param[inout] mf_sstore  intermediate-step scalar variables
 * @param[inout] mf_ru      RHS of total u-velocity
 * @param[inout] mf_rv      RHS of total v-velocity
 * @param[in   ] mf_DU_avg1 time-averaged u-flux for 2D equations
 * @param[in   ] mf_DU_avg2 time-averaged u-flux for 3D equation coupling
 * @param[in   ] mf_DV_avg1 time-averaged v-flux for 2D equations
 * @param[in   ] mf_DV_avg2 time-averaged v-flux for 3D equation coupling
 * @param[inout] mf_ubar    vertically integrated u-momentum
 * @param[inout] mf_vbar    vertically integrated v-momentum
 * @param[inout] mf_Akv     vertical viscosity coefficient
 * @param[inout] mf_Akt     vertical diffusivity coefficient
 * @param[inout] mf_Hz      vertical height of cells
 * @param[inout] mf_Huon    u-volume flux
 * @param[inout] mf_Huon    v-volume flux
 * @param[inout] mf_z_w     vertical coordinates on w points
 * @param[in   ] mf_h       bathymetry
 * @param[in   ] mf_pm      1 / dx
 * @param[in   ] mf_pn      1 / dy
 * @param[in   ] mf_msku    land-sea mask at u-points
 * @param[in   ] mf_mskv    land-sea mask at v-points
 * @param[in   ] N          number of vertical levels
 * @param[in   ] dt_lev     time step at this refinement level
 */
void
REMORA::advance_3d (int lev, MultiFab& mf_cons,
                   MultiFab& mf_u        , MultiFab& mf_v ,
                   MultiFab* mf_sstore,
                   MultiFab* mf_ru       , MultiFab* mf_rv,
                   std::unique_ptr<MultiFab>& mf_DU_avg1,
                   std::unique_ptr<MultiFab>& mf_DU_avg2,
                   std::unique_ptr<MultiFab>& mf_DV_avg1,
                   std::unique_ptr<MultiFab>& mf_DV_avg2,
                   std::unique_ptr<MultiFab>& mf_ubar,
                   std::unique_ptr<MultiFab>& mf_vbar,
                   std::unique_ptr<MultiFab>& mf_Akv,
                   std::unique_ptr<MultiFab>& mf_Akt,
                   std::unique_ptr<MultiFab>& mf_Hz,
                   std::unique_ptr<MultiFab>& mf_Huon,
                   std::unique_ptr<MultiFab>& mf_Hvom,
                   std::unique_ptr<MultiFab>& mf_z_w,
                   MultiFab const* mf_h,
                   MultiFab const* mf_pm,
                   MultiFab const* mf_pn,
                   MultiFab const* mf_mskr,
                   MultiFab const* mf_msku,
                   MultiFab const* mf_mskv,
                   const int N, Real dt_lev)
{
    BL_PROFILE("REMORA::advance_3d()");
    const int nrhs  = 0;
    int nnew  = 0;

    int iic = istep[lev];
    int ntfirst = 0;

    // ROMS accumulates the finer level's perimeter tracer flux over its substeps and
    // subtracts the mismatch against the coarse flux from the coarse ring (correct_tracer,
    // nesting.F). YAFluxRegister is the same construction.
    const bool refluxing = do_reflux && do_substep && finest_level > 0 &&
                           (solverChoice.coupling_type == CouplingType::two_way);

    // Because zeta may have changed
    stretch_transform(lev);

    const BoxArray&            ba = mf_cons.boxArray();
    const DistributionMapping& dm = mf_cons.DistributionMap();

    //Only used locally, probably should be rearranged into FArrayBox declaration
    MultiFab mf_AK (convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW,NGROW,0));       //2d missing j coordinate
    MultiFab mf_DC (ba,dm,1,IntVect(NGROW,NGROW,NGROW-1)); //2d missing j coordinate
    MultiFab mf_Hzk(ba,dm,1,IntVect(NGROW,NGROW,NGROW-1)); //2d missing j coordinate

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real      > const& u = mf_u.array(mfi);
        Array4<Real      > const& v = mf_v.array(mfi);

        Array4<Real      > const& ru = mf_ru->array(mfi);
        Array4<Real      > const& rv = mf_rv->array(mfi);

        Array4<Real      > const& AK = mf_AK.array(mfi);
        Array4<Real      > const& DC = mf_DC.array(mfi);

        Array4<Real      > const& Hzk = mf_Hzk.array(mfi);
        Array4<Real      > const& Akv = mf_Akv->array(mfi);

        Array4<Real const> const& Hz  = mf_Hz->const_array(mfi);

        Array4<Real const> const& DU_avg1  = mf_DU_avg1->const_array(mfi);
        Array4<Real const> const& DV_avg1  = mf_DV_avg1->const_array(mfi);

        Array4<Real const> const& pm  = mf_pm->const_array(mfi);
        Array4<Real const> const& pn  = mf_pn->const_array(mfi);
        Array4<Real const> const& on_u = vec_metric_u[lev]->const_array(mfi,MetricU::on_u);
        Array4<Real const> const& om_v = vec_metric_v[lev]->const_array(mfi,MetricV::om_v);

        Array4<Real const> const& msku  = mf_msku->const_array(mfi);
        Array4<Real const> const& mskv  = mf_mskv->const_array(mfi);

        Box bx = mfi.tilebox();
        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        Box gbx21 = mfi.growntilebox(IntVect(NGROW,NGROW,NGROW-1));

        Box xbx = mfi.nodaltilebox(0);
        Box ybx = mfi.nodaltilebox(1);

        Box gbx2D = gbx2;
        gbx2D.makeSlab(2,0);

        Box tbxp1 = bx;
        Box tbxp11 = bx;
        Box tbxp2 = bx;
        tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
        tbxp2.grow(IntVect(NGROW,NGROW,0));
        tbxp11.grow(IntVect(NGROW-1,NGROW-1,NGROW-1));

        FArrayBox fab_FC(convert(gbx2,IntVect(0,0,1)),1,amrex::The_Async_Arena());
        FArrayBox fab_BC(convert(gbx2,IntVect(0,0,1)),1,amrex::The_Async_Arena());
        FArrayBox fab_CF(gbx21,1,amrex::The_Async_Arena());
        FArrayBox fab_W(tbxp2,1,amrex::The_Async_Arena());

        auto FC = fab_FC.array();
        auto BC = fab_BC.array();
        auto CF = fab_CF.array();

        Real cff;
        if (iic==ntfirst) {
            cff=Real(0.25)*dt_lev;
        } else if (iic==ntfirst+1) {
            cff=Real(0.25)*dt_lev*Real(3.0)/Real(2.0);
        } else {
            cff=Real(0.25)*dt_lev*Real(23.0)/Real(12.0);
        }

        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            u(i,j,k) += cff * (pm(i,j,0)+pm(i-1,j,0)) * (pn(i,j,0)+pn(i-1,j,0)) * ru(i,j,k,nrhs);
            u(i,j,k) *= two / (Hz(i-1,j,k) + Hz(i,j,k));
        });

        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            v(i,j,k) += cff * (pm(i,j,0)+pm(i,j-1,0)) * (pn(i,j,0)+pn(i,j-1,0)) * rv(i,j,k,nrhs);
            v(i,j,k) *= two / (Hz(i,j-1,k) + Hz(i,j,k));
        });

        // NOTE: DC is only used as scratch in vert_visc_3d -- no need to pass or return a value
        // NOTE: may not actually need to set these to zero

        // Reset to zero on the box on which they'll be used
        mf_DC[mfi].template setVal<RunOn::Device>(zero,xbx);
        fab_CF.template     setVal<RunOn::Device>(zero,xbx);

        vert_visc_3d(xbx,1,0,u,Hz,Hzk,AK,Akv,BC,DC,FC,CF,nnew,N,dt_lev);

        // Reset to zero on the box on which they'll be used
        mf_DC[mfi].template setVal<RunOn::Device>(zero,ybx);
        fab_CF.template     setVal<RunOn::Device>(zero,ybx);

        vert_visc_3d(ybx,0,1,v,Hz,Hzk,AK,Akv,BC,DC,FC,CF,nnew,N,dt_lev);

        // Reset to zero on the box on which they'll be used
        mf_DC[mfi].template setVal<RunOn::Device>(zero,xbx);
        fab_CF.template     setVal<RunOn::Device>(zero,xbx);

        vert_mean_3d(xbx,1,0,u,Hz,DU_avg1,DC,CF,on_u,msku,nnew,N);

        // Reset to zero on the box on which they'll be used
        mf_DC[mfi].template setVal<RunOn::Device>(zero,ybx);
        fab_CF.template     setVal<RunOn::Device>(zero,ybx);

        vert_mean_3d(ybx,0,1,v,Hz,DV_avg1,DC,CF,om_v,mskv,nnew,N);
    }

    // Apply physical boundary conditions to u and v
    (*physbcs[lev])(mf_u,*mf_msku,0,1,mf_u.nGrowVect(),t_old[lev],xvel_bc(),0,*xvel_old[lev]);
    (*physbcs[lev])(mf_v,*mf_mskv,0,1,mf_v.nGrowVect(),t_old[lev],yvel_bc(),0,*yvel_old[lev]);

#ifdef REMORA_USE_NETCDF
        // Fill the data which is stored in the boundary data read from netcdf files
        if (solverChoice.boundary_from_netcdf)
        {
            fill_from_bdyfiles(lev, mf_u,*mf_msku,t_old[lev],xvel_bc(),BdyVars::u,0,0,*xvel_old[lev],dt_lev);
            fill_from_bdyfiles(lev, mf_v,*mf_mskv,t_old[lev],yvel_bc(),BdyVars::v,0,0,*yvel_old[lev],dt_lev);
        }

    if (solverChoice.do_rivers) {
        river_source_transport->update_interpolated_to_time(model_time(t_old[lev]));
    }
#endif
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real      > const& u = mf_u.array(mfi);
        Array4<Real      > const& v = mf_v.array(mfi);

        Array4<Real      > const& DC = mf_DC.array(mfi);

        Array4<Real const> const& Hz  = mf_Hz->const_array(mfi);

        Array4<Real const> const& DU_avg1  = mf_DU_avg1->const_array(mfi);
        Array4<Real const> const& DV_avg1  = mf_DV_avg1->const_array(mfi);

        Array4<Real const> const& DU_avg2  = mf_DU_avg2->const_array(mfi);
        Array4<Real const> const& DV_avg2  = mf_DV_avg2->const_array(mfi);

        Array4<Real> const& ubar = mf_ubar->array(mfi);
        Array4<Real> const& vbar = mf_vbar->array(mfi);

        Array4<Real> const& Huon = mf_Huon->array(mfi);
        Array4<Real> const& Hvom = mf_Hvom->array(mfi);

        Array4<Real const> const& on_u = vec_metric_u[lev]->const_array(mfi,MetricU::on_u);
        Array4<Real const> const& om_v = vec_metric_v[lev]->const_array(mfi,MetricV::om_v);
        Array4<Real const> const& msku = mf_msku->const_array(mfi);
        Array4<Real const> const& mskv = mf_mskv->const_array(mfi);

        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));

        FArrayBox fab_FC(gbx2,1,amrex::The_Async_Arena());
        auto FC = fab_FC.array();

#ifdef REMORA_USE_NETCDF
        if (solverChoice.do_rivers) {
            Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
            Array4<Real const> const& z_w = mf_z_w->const_array(mfi);
            Array4<int  const> const& river_pos = vec_river_position[lev]->const_array(mfi);
            Array4<Real const> const& river_transport = river_source_transport->fab_interp->array();
            int* river_direction_d = river_direction.data();
            ParallelFor(gbx1, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                int iriver = river_pos(i,j,0);
                if (iriver >= 0) {
                    if (river_direction_d[iriver] == 0) {
                        Real cff = one / (on_u(i,j,0) * Real(0.5) * (z_w(i-1,j,k+1) - z_w(i-1,j,k) + z_w(i,j,k+1) - z_w(i,j,k)));
                        u(i,j,k) = cff * river_transport(iriver,0,k);
                    } else {
                        Real cff = one / (om_v(i,j,0) * Real(0.5) * (z_w(i,j-1,k+1) - z_w(i,j-1,k) + z_w(i,j,k+1) - z_w(i,j,k)));
                        v(i,j,k) = cff * river_transport(iriver,0,k);
                    }
                }
            });
        }
#endif

#if 0
        // Reset to zero on the box on which they'll be used
        mf_DC[mfi].template setVal<RunOn::Device>(zero,grow(xbx,IntVect(0,0,1)));
        fab_CF.template     setVal<RunOn::Device>(zero,grow(xbx,IntVect(0,0,1)));

        update_massflux_3d(xbx,1,0,u,ubar,Huon,Hz,pn,DU_avg1,DU_avg2,DC,FC,nnew);

        // Reset to zero on the box on which they'll be used
        mf_DC[mfi].template setVal<RunOn::Device>(zero,grow(ybx,IntVect(0,0,1)));
        fab_CF.template     setVal<RunOn::Device>(zero,grow(ybx,IntVect(0,0,1)));

        update_massflux_3d(ybx,0,1,v,vbar,Hvom,Hz,pm,DV_avg1,DV_avg2,DC,FC,nnew);

#else

        // Reset to zero on the box on which they'll be used
        fab_FC.template setVal<RunOn::Device>(zero,gbx2);
        mf_DC[mfi].template setVal<RunOn::Device>(zero,grow(gbx2,IntVect(0,0,1)));
        update_massflux_3d(lev,gbx2,1,0,u,ubar,Huon,Hz,on_u,DU_avg1,DU_avg2,DC,FC,msku,nnew);

        // Reset to zero on the box on which they'll be used
        fab_FC.template     setVal<RunOn::Device>(zero,gbx2);
        mf_DC[mfi].template setVal<RunOn::Device>(zero,grow(gbx2,IntVect(0,0,1)));
        update_massflux_3d(lev,gbx2,0,1,v,vbar,Hvom,Hz,om_v,DV_avg1,DV_avg2,DC,FC,mskv,nnew);

#endif
    }

    // WE BELIEVE THESE VALUES SHOULD ALREADY BE FILLED
    // mf_Huon->FillBoundary(geom[lev].periodicity());
    // mf_Hvom->FillBoundary(geom[lev].periodicity());

    // ************************************************************************
    // This should fill both temp and salt with temp/salt currently in cons_old
    // ************************************************************************

    MultiFab mf_W(convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW+1,NGROW+1,0));
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& Huon = mf_Huon->array(mfi);
        Array4<Real> const& Hvom = mf_Hvom->array(mfi);

        Array4<Real const> const& z_w = mf_z_w->const_array(mfi);
        Array4<Real const> const& h   = mf_h->const_array(mfi);

        Box bx = mfi.tilebox();
        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        Box gbx21 = mfi.growntilebox(IntVect(NGROW,NGROW,NGROW-1));

        Box tbxp1 = bx;
        Box tbxp11 = bx;
        Box tbxp2 = bx;
        tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
        tbxp2.grow(IntVect(NGROW,NGROW,0));
        tbxp11.grow(IntVect(NGROW-1,NGROW-1,NGROW-1));

        FArrayBox fab_FC(surroundingNodes(gbx2,2),1,amrex::The_Async_Arena());
        FArrayBox fab_BC(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_CF(gbx21,1,amrex::The_Async_Arena());

        auto W   = mf_W.array(mfi);

        //
        //------------------------------------------------------------------------
        //  Vertically integrate horizontal mass flux divergence.
        //------------------------------------------------------------------------
        //
        //Should really use gbx3uneven
        //TODO: go over these boxes and compare to other spots where we do the same thing
        Box gbx1D = gbx1;
        gbx1D.makeSlab(2,0);

        compute_W(gbx1D, W, Huon, Hvom, z_w, h, N);
    }

    const int nstp = (iic) % 2;
    nnew = 1-nstp;
    if (solverChoice.vert_mixing_type == VertMixingType::GLS) {
        gls_corrector(lev, vec_gls[lev].get(), vec_tke[lev].get(), mf_W, vec_Akv[lev].get(),
                  vec_Akt[lev].get(),vec_Akk[lev].get(), vec_Akp[lev].get(), vec_mskr[lev].get(),
                  vec_msku[lev].get(), vec_mskv[lev].get(),
                  nstp, nnew, N, dt_lev);
    }
    advance_biology(lev, *cons_old[lev], mf_cons, N, dt_lev);
    nnew = 0;

    // FineAdd accumulates a tile's fluxes into the coarse cells under it with atomic adds,
    // and a tile edge that is not a multiple of the refinement ratio shares a coarse cell
    // with the neighbouring tile, so with several threads the order of the adds, and the
    // rounding, would follow thread timing. With more than one thread the fluxes are kept
    // on this level's faces and added after the loop, tile by tile in the fixed order a
    // single thread uses.
    const bool defer_fine_add = refluxing && lev > 0 && OpenMP::get_max_threads() > 1;
    if (defer_fine_add) {
        if (int(vec_reflux_flux.size()) <= lev) { vec_reflux_flux.resize(lev+1); }
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            MultiFab& f = vec_reflux_flux[lev][d];
            const BoxArray fba = convert(mf_cons.boxArray(), IntVect::TheDimensionVector(d));
            if (!f.ok() || f.boxArray() != fba || f.DistributionMap() != mf_cons.DistributionMap()
                || f.nComp() != ncons) {
                f.define(fba, mf_cons.DistributionMap(), ncons, 0);
            }
        }
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& Hz  = mf_Hz->array(mfi);

        Array4<Real> const& Huon = mf_Huon->array(mfi);
        Array4<Real> const& Hvom = mf_Hvom->array(mfi);

        Array4<Real const> const& pm  = mf_pm->const_array(mfi);
        Array4<Real const> const& pn  = mf_pn->const_array(mfi);
        Array4<Real const> const& mskr  = mf_mskr->const_array(mfi);
        Array4<Real const> const& msku  = mf_msku->const_array(mfi);
        Array4<Real const> const& mskv  = mf_mskv->const_array(mfi);

        Box bx = mfi.tilebox();
        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        Box gbx21 = mfi.growntilebox(IntVect(NGROW,NGROW,NGROW-1));

        Box tbxp1 = bx;
        Box tbxp11 = bx;
        Box tbxp2 = bx;
        tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
        tbxp2.grow(IntVect(NGROW,NGROW,0));
        tbxp11.grow(IntVect(NGROW-1,NGROW-1,NGROW-1));

        FArrayBox fab_FC(surroundingNodes(gbx2,2),1,amrex::The_Async_Arena());
        FArrayBox fab_BC(gbx2,1,amrex::The_Async_Arena());
        FArrayBox fab_CF(gbx21,1,amrex::The_Async_Arena());

        auto FC  = fab_FC.array();
        auto W   = mf_W.array(mfi);

        // Advective tracer fluxes for the coarse-fine correction, all components at once.
        FArrayBox fab_fx_reg, fab_fy_reg, fab_fz_reg;
        if (refluxing) {
            fab_fx_reg.resize(surroundingNodes(bx,0), ncons, The_Async_Arena());
            fab_fy_reg.resize(surroundingNodes(bx,1), ncons, The_Async_Arena());
            fab_fz_reg.resize(surroundingNodes(bx,2), ncons, The_Async_Arena());
            fab_fx_reg.template setVal<RunOn::Device>(zero);
            fab_fy_reg.template setVal<RunOn::Device>(zero);
            // No refinement in the vertical, so nothing crosses a z interface.
            fab_fz_reg.template setVal<RunOn::Device>(zero);
        }
        //
        //-----------------------------------------------------------------------
        // rhs_t_3d
        //-----------------------------------------------------------------------
        //
        for (int i_comp=0; i_comp < ncons; i_comp++)
        {
#ifdef REMORA_USE_NETCDF
            FArrayBox* fab_river_source;
            if (solverChoice.do_rivers_cons[i_comp]) {
                fab_river_source = river_source_cons[i_comp]->fab_interp;
            }
            const Array4<const int>& river_pos = (solverChoice.do_rivers) ? vec_river_position[lev]->const_array(mfi) : Array4<const int>();
            const Array4<const Real>& river_source = (solverChoice.do_rivers_cons[i_comp]) ? fab_river_source->const_array() : Array4<const Real>();
#else
            const Array4<const int >& river_pos = Array4<const int>();
            const Array4<const Real>& river_source = Array4<const Real>();
#endif
            Array4<Real> const& sstore = mf_sstore->array(mfi, i_comp);
            rhs_t_3d(lev,bx, mf_cons.array(mfi,i_comp), sstore, Huon, Hvom,
                     Hz, pn, pm, W, FC, mskr, msku, mskv, river_pos, river_source, nrhs, nnew, N,dt_lev,
                     refluxing ? fab_fx_reg.array(i_comp) : Array4<Real>(),
                     refluxing ? fab_fy_reg.array(i_comp) : Array4<Real>());
        }

        if (refluxing) {
            std::array<FArrayBox const*, AMREX_SPACEDIM> flux
                { &fab_fx_reg, &fab_fy_reg, &fab_fz_reg };
            const Real* dx = geom[lev].CellSize();
            if (lev < finest_level) {
                getAdvFluxReg(lev+1)->CrseAdd(mfi, flux, dx, dt_lev, RunOn::Device);
            }
            if (lev > 0 && !defer_fine_add) {
                getAdvFluxReg(lev)->FineAdd(mfi, flux, dx, dt_lev, RunOn::Device);
            } else if (lev > 0) {
                // nodaltilebox gives each face to exactly one tile
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    const Box& fbx = mfi.nodaltilebox(d);
                    Array4<Real      > const& dst = vec_reflux_flux[lev][d].array(mfi);
                    Array4<Real const> const& src = flux[d]->const_array();
                    ParallelFor(fbx, ncons, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                    {
                        dst(i,j,k,n) = src(i,j,k,n);
                    });
                }
            }
            // The flux FABs are temporaries scoped to this tile, so don't move on until the
            // register has finished reading them. ERF does the same at its own CrseAdd.
            Gpu::streamSynchronize();
        }
    } // mfi

    if (defer_fine_add) {
        const Real* dx = geom[lev].CellSize();
        for ( MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            // FineAdd takes every face the flux FAB holds inside the coarsened tile, so a
            // tile that is not a multiple of the refinement ratio would also add its
            // neighbour's faces from the level-wide FAB. Hand it the tile's own faces, as
            // the loop above did.
            std::array<FArrayBox, AMREX_SPACEDIM> ftile;
            std::array<FArrayBox const*, AMREX_SPACEDIM> flux;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                ftile[d].resize(surroundingNodes(mfi.tilebox(), d), ncons, The_Async_Arena());
                ftile[d].template copy<RunOn::Device>(vec_reflux_flux[lev][d][mfi], 0, 0, ncons);
                flux[d] = &ftile[d];
            }
            getAdvFluxReg(lev)->FineAdd(mfi, flux, dx, dt_lev, RunOn::Device);
            Gpu::streamSynchronize();
        }
    }

    FillPatch(lev, t_old[lev], mf_cons, cons_new, BCVars::cons_bc, BdyVars::t,0,true,false,0,0,dt_lev,*cons_old[lev]);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& AK = mf_AK.array(mfi);
        Array4<Real> const& DC = mf_DC.array(mfi);

        Array4<Real> const& Hzk = mf_Hzk.array(mfi);
        Array4<Real const> const& Hz  = mf_Hz->const_array(mfi);

        Box bx = mfi.tilebox();

        // Copy the tilebox
        Box tbxp1 = bx;
        Box tbxp11 = bx;
        Box tbxp2 = bx;
        Box tbxp21 = bx;
        //make only gbx be grown to match multifabs
        tbxp21.grow(IntVect(NGROW,NGROW,NGROW-1));
        tbxp2.grow(IntVect(NGROW,NGROW,0));
        tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
        tbxp11.grow(IntVect(NGROW-1,NGROW-1,NGROW-1));

        FArrayBox fab_FC(tbxp2,1,amrex::The_Async_Arena());
        FArrayBox fab_BC(tbxp2,1,amrex::The_Async_Arena());
        FArrayBox fab_CF(tbxp21,1,amrex::The_Async_Arena());
        FArrayBox fab_W(tbxp2,1,amrex::The_Async_Arena());

        auto FC = fab_FC.array();
        auto BC = fab_BC.array();
        auto CF = fab_CF.array();

        for (int i_comp=0; i_comp < ncons; i_comp++) {
            vert_visc_3d(bx,0,0,mf_cons.array(mfi,i_comp),Hz,Hzk,
                    AK,mf_Akt->array(mfi,akt_comp(i_comp)),BC,DC,FC,CF,nnew,N,dt_lev);
        }
    } // MFiter
    FillPatch(lev, t_old[lev], *cons_new[lev], cons_new, BCVars::cons_bc, BdyVars::t,0,true,false,0,0,dt_lev,*cons_old[lev]);

#ifdef REMORA_USE_NETCDF
    // Nudge every tracer that was given a climatology toward it. temp and salt are just
    // components 0 and 1 of this loop, so their behavior is unchanged.
    for (int icomp = 0; icomp < ncons; ++icomp) {
        if (!solverChoice.do_cons_clim_nudg[icomp]) { continue; }

        cons_clim_data_from_file[icomp]->update_interpolated_to_time(model_time(t_old[lev]), lev, cons_new[lev], geom, ref_ratio);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(mf_cons, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            Box bx = mfi.growntilebox(IntVect(1,1,0));
            Array4<const Real> const& cons_nudg_coeff = vec_nudg_coeff[BdyVars::cons(icomp)][lev]->const_array(mfi);
            Array4<const Real> const& cons_clim = cons_clim_data_from_file[icomp]->get_interpolated_mf(lev)->const_array(mfi);
            Array4<      Real> const& cons = cons_new[lev]->array(mfi, icomp);
            Array4<const Real> const& Hz   = vec_Hz[lev]->const_array(mfi);

            // Tracers are nudged at rho points: no face metrics
            apply_clim_nudg(bx, 0, 0, cons, cons, cons_clim, cons_nudg_coeff, Hz,
                            Array4<const Real>(), Array4<const Real>(), dt_lev);
        }
    }
#endif
}
