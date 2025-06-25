#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] lev            level to operate on
 * @param[in   ] mf_uold        u-velocity from last time step
 * @param[in   ] mf_vold        v-velocity from last time step
 * @param[  out] mf_u           u-velocity at current time step
 * @param[  out] mf_v           v-velocity at current time step
 * @param[inout] mf_ru          u-velocity RHS at current time step
 * @param[inout] mf_rv          v-velocity RHS at current time step
 * @param[in   ] S_old          scalar variables at last time step
 * @param[inout] S_new          scalar variables at current time step
 * @param[in   ] mf_W           vertical velocity
 * @param        mf_DC          temporary variable container
 * @param[in   ] mf_z_r         z coordinates at rho points (cell centers)
 * @param[in   ] mf_z_w         z coordinates at w points
 * @param[in   ] mf_h           bathymetry
 * @param[in   ] mf_pm          1/dx
 * @param[in   ] mf_pn          1/dy
 * @param[in   ] mf_sustr       u-direction surface momentum flux
 * @param[in   ] mf_svstr       v-direction surface momentum flux
 * @param[in   ] mf_bustr       u-direction bottom stress
 * @param[in   ] mf_bvstr       v-direction bottom stress
 * @param[in   ] mf_msku        land-sea mask on u-points
 * @param[in   ] mf_mskv        land-sea mask on v-points
 * @param[in   ] iic            which time step we're on
 * @param[in   ] ntfirst        what is the first time step?
 * @param[in   ] nnew           index of time step to update
 * @param[in   ] nstp           index of last time step
 * @param[in   ] nrhs           index of RHS component
 * @param[in   ] N              number of vertical levels
 * @param[in   ] dt_lev         time step at this level
 */

void
REMORA::prestep (int lev,
                MultiFab& mf_uold, MultiFab& mf_vold,
                MultiFab& mf_u, MultiFab& mf_v,
                      MultiFab* mf_ru,
                      MultiFab* mf_rv,
                MultiFab& S_old, MultiFab& S_new,
                MultiFab& mf_W, MultiFab& mf_DC,
                const MultiFab* mf_z_r,
                const MultiFab* mf_z_w,
                const MultiFab* mf_h,
                const MultiFab* mf_pm,
                const MultiFab* mf_pn,
                const MultiFab* mf_sustr,
                const MultiFab* mf_svstr,
                const MultiFab* mf_bustr,
                const MultiFab* mf_bvstr,
                const MultiFab* mf_msku,
                const MultiFab* mf_mskv,
                const int iic, const int ntfirst,
                const int nnew, int nstp, int nrhs,
                int N, const Real dt_lev)
{
    const BoxArray&            ba = S_old.boxArray();
    const DistributionMapping& dm = S_old.DistributionMap();

    // Maybe not the best way to do this, but need to cache salt and temp since
    // they get rewritten by prestep_t
    MultiFab mf_saltcache(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_tempcache(ba,dm,1,IntVect(NGROW,NGROW,0));

    MultiFab mf_scalarcache(ba,dm,NCONS,IntVect(NGROW,NGROW,0));
    MultiFab::Copy(mf_scalarcache,S_new,0,0,NCONS,IntVect(NGROW,NGROW,0));

    MultiFab::Copy(mf_saltcache,S_new,Salt_comp,0,1,IntVect(NGROW,NGROW,0));
    MultiFab::Copy(mf_tempcache,S_new,Temp_comp,0,1,IntVect(NGROW,NGROW,0));

    for ( MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& DC = mf_DC.array(mfi);
        Array4<Real> const& Akv   = vec_Akv[lev]->array(mfi);
        Array4<Real> const& Akt   = vec_Akt[lev]->array(mfi);
        Array4<Real> const& Hz    = vec_Hz[lev]->array(mfi);
        Array4<Real> const& Huon  = vec_Huon[lev]->array(mfi);
        Array4<Real> const& Hvom  = vec_Hvom[lev]->array(mfi);
        Array4<Real const> const& uold  = mf_uold.const_array(mfi);
        Array4<Real const> const& vold = mf_vold.const_array(mfi);
        Array4<Real> const& u = (mf_u).array(mfi);
        Array4<Real> const& v = (mf_v).array(mfi);

        Array4<Real const> const& z_r   = mf_z_r->const_array(mfi);
        Array4<Real const> const& z_w   = mf_z_w->const_array(mfi);
        Array4<Real const> const& h     = mf_h->const_array(mfi);
        Array4<Real const> const& pm    = mf_pm->const_array(mfi);
        Array4<Real const> const& pn    = mf_pn->const_array(mfi);

        Array4<Real      > const& ru = mf_ru->array(mfi);
        Array4<Real      > const& rv = mf_rv->array(mfi);
        Array4<Real      > const& W = (mf_W).array(mfi);
        Array4<Real const> const& sustr = mf_sustr->const_array(mfi);
        Array4<Real const> const& svstr = mf_svstr->const_array(mfi);
        Array4<Real const> const& bustr = mf_bustr->const_array(mfi);
        Array4<Real const> const& bvstr = mf_bvstr->const_array(mfi);
        Array4<Real const> const& msku  = mf_msku->const_array(mfi);
        Array4<Real const> const& mskv  = mf_mskv->const_array(mfi);

        Real lambda = 1.0_rt;

        Box bx = mfi.tilebox();
        Box gbx = mfi.growntilebox();
        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));
        Box gbx2 = mfi.growntilebox(IntVect(NGROW,NGROW,0));
        //Box gbx11 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,NGROW-1));

        //copy the tilebox
        Box tbxp1 = bx;
        Box tbxp11 = bx;
        Box tbxp2 = bx;

        //TODO: adjust for tiling
        //Box gbx3uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-3,bx.smallEnd(1)-3,bx.smallEnd(2))),
        //               IntVect(AMREX_D_DECL(bx.bigEnd(0)+2,bx.bigEnd(1)+2,bx.bigEnd(2))));
        //Box gbx2uneven(IntVect(AMREX_D_DECL(bx.smallEnd(0)-2,bx.smallEnd(1)-2,bx.smallEnd(2))),
        //               IntVect(AMREX_D_DECL(bx.bigEnd(0)+1,bx.bigEnd(1)+1,bx.bigEnd(2))));
        //make only gbx be grown to match multifabs
        tbxp2.grow(IntVect(NGROW,NGROW,0));
        tbxp1.grow(IntVect(NGROW-1,NGROW-1,0));
        tbxp11.grow(IntVect(NGROW-1,NGROW-1,NGROW-1));

        Box bxD = bx;
        bxD.makeSlab(2,0);
        Box gbx1D = gbx1;
        gbx1D.makeSlab(2,0);
        Box gbx2D = gbx2;
        gbx2D.makeSlab(2,0);

        Box tbxp1D = tbxp1;
        tbxp1D.makeSlab(2,0);
        Box tbxp2D = tbxp2;
        tbxp2D.makeSlab(2,0);

        FArrayBox fab_FC(convert(tbxp2,IntVect(0,0,1)),1,amrex::The_Async_Arena()); //3D

        auto FC=fab_FC.array();


        for (int i_comp=0; i_comp < NCONS; i_comp++) {
            Array4<Real> const& sstore = (vec_sstore[lev])->array(mfi,i_comp);
#ifdef REMORA_USE_NETCDF
            FArrayBox* fab_river_source;
            if (solverChoice.do_rivers_cons[i_comp]) {
                river_source_cons[i_comp]->update_interpolated_to_time(t_old[lev]);
                fab_river_source = river_source_cons[i_comp]->fab_interp;
            }
            const Array4<const int>& river_pos = (solverChoice.do_rivers) ? vec_river_position[lev]->const_array(mfi) : Array4<const int>();
            const Array4<const Real>& river_source = (solverChoice.do_rivers_cons[i_comp]) ? fab_river_source->const_array() : Array4<const Real>();
#else
            const Array4<const int >& river_pos = Array4<const int>();
            const Array4<const Real>& river_source = Array4<const Real>();
#endif
            prestep_t_advection(bx, gbx, S_old.array(mfi,i_comp),
                                mf_scalarcache.array(mfi,i_comp), Hz, Huon, Hvom,
                                W, DC, FC, sstore, z_w, h, pm, pn, msku, mskv, river_pos,
                                river_source, iic, ntfirst,
                                nrhs, N, dt_lev);
        }

        // Only do diffusion for salt and temperature, not other tracer(s)
        for (int i_comp=0; i_comp < NCONS; i_comp++) {
            const Array4<Real const>& stflx = vec_stflx[lev]->const_array(mfi,i_comp);
            const Array4<Real const>& btflx = vec_btflx[lev]->const_array(mfi,i_comp);
            prestep_diffusion(bx,gbx,0,0,S_new.array(mfi,i_comp), S_old.array(mfi,i_comp), ru,
                              Hz, Akt, DC, FC, stflx, btflx, z_r, pm, pn, iic, iic, nnew, nstp,
                              nrhs, N, lambda, dt_lev);
        }

        //
        //-----------------------------------------------------------------------
        // prestep_uv_3d
        //-----------------------------------------------------------------------
        //
        //updates u,v,ru,rv (ru and rv have multiple components)

        prestep_diffusion(tbxp1, gbx, 1, 0, u, uold, ru, Hz, Akv, DC, FC,
                          sustr, bustr, z_r, pm, pn, iic, ntfirst, nnew, nstp, nrhs, N, lambda, dt_lev);

        prestep_diffusion(tbxp1, gbx, 0, 1, v, vold, rv, Hz, Akv, DC, FC,
                          svstr, bvstr, z_r, pm, pn, iic, ntfirst, nnew, nstp, nrhs, N, lambda, dt_lev);
    }
}
