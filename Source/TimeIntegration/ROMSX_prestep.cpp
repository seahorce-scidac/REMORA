#include <ROMSX.H>

using namespace amrex;

/**
 * prestep
 *
 * @param[in   ] lev
 * @param[in   ] mf_uold
 * @param[in   ] mf_vold
 * @param[inout] mf_u (looks like reset in update_vel <- prestep_uv_3d, so maybe just out
 * @param[inout] mf_v maybe just out
 * @param[inout] mf_ru
 * @param[inout] mf_rv
 * @param[in   ] S_old
 * @param[inout] S_new
 * @param[out  ] mf_W
 * @param[none ] mf_DC (temp)
 * @param[in   ] mf_z_r
 * @param[in   ] mf_z_w
 * @param[in   ] mf_h
 * @param[in   ] mf_sustr
 * @param[in   ] mf_svstr
 * @param[in   ] mf_bustr
 * @param[in   ] mf_bvstr
 * @param[in   ] iic
 * @param[in   ] ntfirst
 * @param[in   ] nnew
 * @param[in   ] nstp
 * @param[in   ] nrhs
 * @param[in   ] N
 * @param[in   ] dt_lev
 */

void
ROMSX::prestep (int lev,
                MultiFab& mf_uold, MultiFab& mf_vold,
                MultiFab& mf_u, MultiFab& mf_v,
                std::unique_ptr<MultiFab>& mf_ru,
                std::unique_ptr<MultiFab>& mf_rv,
                MultiFab& S_old, MultiFab& S_new,
                MultiFab& mf_W, MultiFab& mf_DC,
                std::unique_ptr<MultiFab>& mf_z_r,
                std::unique_ptr<MultiFab>& mf_z_w,
                std::unique_ptr<MultiFab>& mf_h,
                std::unique_ptr<MultiFab>& mf_sustr,
                std::unique_ptr<MultiFab>& mf_svstr,
                std::unique_ptr<MultiFab>& mf_bustr,
                std::unique_ptr<MultiFab>& mf_bvstr,
                const int iic, const int ntfirst,
                const int nnew, int nstp, int nrhs,
                int N, const Real dt_lev)
{
    const auto dxi              = Geom(lev).InvCellSizeArray();

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
        Array4<Real> const& Akv = (vec_Akv[lev])->array(mfi);
        Array4<Real> const& Akt = (vec_Akt[lev])->array(mfi);
        Array4<Real> const& Hz  = (vec_Hz[lev])->array(mfi);
        Array4<Real> const& Huon  = (vec_Huon[lev])->array(mfi);
        Array4<Real> const& Hvom  = (vec_Hvom[lev])->array(mfi);
        Array4<Real> const& z_r = (mf_z_r)->array(mfi);
        Array4<Real> const& z_w= (mf_z_w)->array(mfi);
        Array4<Real> const& h= (mf_h)->array(mfi);
        Array4<Real> const& uold = (mf_uold).array(mfi);
        Array4<Real> const& vold = (mf_vold).array(mfi);
        Array4<Real> const& u = (mf_u).array(mfi);
        Array4<Real> const& v = (mf_v).array(mfi);

        Array4<Real> const& ru = (mf_ru)->array(mfi);
        Array4<Real> const& rv = (mf_rv)->array(mfi);
        Array4<Real> const& W = (mf_W).array(mfi);
        Array4<Real> const& sustr = (mf_sustr)->array(mfi);
        Array4<Real> const& svstr = (mf_svstr)->array(mfi);
        Array4<Real> const& bustr = (mf_bustr)->array(mfi);
        Array4<Real> const& bvstr = (mf_bvstr)->array(mfi);

        Real lambda = 1.0;

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

        FArrayBox fab_FC(tbxp2,1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_pm(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pn(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_on_u(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_om_v(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_om_u(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_on_v(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_om_r(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_on_r(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_om_p(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_on_p(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pmon_u(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pnom_u(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pmon_v(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pnom_v(tbxp2D,1,amrex::The_Async_Arena());
        //FArrayBox fab_oHz(gbx11,1,amrex::The_Async_Arena());

        auto FC=fab_FC.array();
        auto pm=fab_pm.array();
        auto pn=fab_pn.array();


        ParallelFor(tbxp2D,
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            //Note: are the comment definitons right? Don't seem to match metrics.f90
            pm(i,j,0) = dxi[0];
            pn(i,j,0) = dxi[1];
        });

        //From ini_fields and .in file
        //fab_Akt.setVal(1e-6);
        FArrayBox fab_stflux(tbxp2,1,amrex::The_Async_Arena());
        auto stflux= fab_stflux.array();
        FArrayBox fab_btflux(tbxp2,1,amrex::The_Async_Arena());
        auto btflux= fab_btflux.array();

        //From ini_fields and .in file
        //fab_stflux.setVal(0.0);
        //also set btflux=0 (as in ana_btflux.H)
        ParallelFor(tbxp2,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            stflux(i,j,k)=0.0;
            btflux(i,j,k)=0.0;
        });

        for (int i_comp=0; i_comp < NCONS; i_comp++) {
            Array4<Real> const& sstore = (vec_sstore[lev])->array(mfi,i_comp);
            prestep_t_advection(bx, gbx, S_old.array(mfi,i_comp),
                                mf_scalarcache.array(mfi,i_comp), Hz, Huon, Hvom,
                                pm, pn, W, DC, FC, sstore, z_w, h, iic, ntfirst,
                                nrhs, N, dt_lev);
        }

        // Only do diffusion for salt and temperature, not other tracer(s)
        for (int i_comp=0; i_comp < NCONS; i_comp++) {
            prestep_diffusion(bx,gbx,0,0,S_new.array(mfi,i_comp), S_old.array(mfi,i_comp), ru,
                              Hz, Akt, DC, FC, stflux, btflux, z_r, pm, pn, iic, iic, nnew, nstp,
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
