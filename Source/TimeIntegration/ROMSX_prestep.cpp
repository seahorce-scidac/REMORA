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
    const auto prob_lo          = Geom(lev).ProbLoArray();
    const auto dxi              = Geom(lev).InvCellSizeArray();
    const auto dx               = Geom(lev).CellSizeArray();
    const int Mm = Geom(lev).Domain().size()[1];

    const BoxArray&            ba = S_old.boxArray();
    const DistributionMapping& dm = S_old.DistributionMap();

    // Maybe not the best way to do this, but need to cache salt and temp since
    // they get rewritten by prestep_t
    MultiFab mf_saltcache(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_tempcache(ba,dm,1,IntVect(NGROW,NGROW,0));

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

        Array4<Real> const& tempold = S_old.array(mfi,Temp_comp);
        Array4<Real> const& saltold = S_old.array(mfi,Salt_comp);
        Array4<Real> const& temp    = S_new.array(mfi,Temp_comp);
        Array4<Real> const& salt    = S_new.array(mfi,Salt_comp);

        Array4<Real> const& tempstore = (vec_t3[lev])->array(mfi);
        Array4<Real> const& saltstore = (vec_s3[lev])->array(mfi);
        Array4<Real> const& ru = (mf_ru)->array(mfi);
        Array4<Real> const& rv = (mf_rv)->array(mfi);
        Array4<Real> const& W = (mf_W).array(mfi);
        Array4<Real> const& sustr = (mf_sustr)->array(mfi);
        Array4<Real> const& svstr = (mf_svstr)->array(mfi);
        Array4<Real> const& bustr = (mf_bustr)->array(mfi);
        Array4<Real> const& bvstr = (mf_bvstr)->array(mfi);
        Array4<Real> const& tempcache = (mf_tempcache).array(mfi);
        Array4<Real> const& saltcache = (mf_saltcache).array(mfi);

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
        auto on_u=fab_on_u.array();
        auto om_v=fab_om_v.array();
        auto om_u=fab_om_u.array();
        auto on_v=fab_on_v.array();
        auto om_r=fab_om_r.array();
        auto on_r=fab_on_r.array();
        auto om_p=fab_om_p.array();
        auto on_p=fab_on_p.array();
        auto pmon_u=fab_pmon_u.array();
        auto pnom_u=fab_pnom_u.array();
        auto pmon_v=fab_pmon_v.array();
        auto pnom_v=fab_pnom_v.array();


        ParallelFor(tbxp2D,
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            //Note: are the comment definitons right? Don't seem to match metrics.f90
            pm(i,j,0) = dxi[0];
            pn(i,j,0) = dxi[1];
        });

        prestep_t_3d(bx, gbx, tempold, temp, tempcache, ru, Hz, Akt, Huon, Hvom,
                     pm, pn, W, DC, FC, tempstore, z_r, z_w, h, iic, ntfirst, nnew, nstp, nrhs, N,
                     lambda, dt_lev);
        prestep_t_3d(bx, gbx, saltold, salt, saltcache, ru, Hz, Akt, Huon, Hvom,
                     pm, pn, W, DC, FC, saltstore, z_r, z_w, h, iic, ntfirst, nnew, nstp, nrhs, N,
                     lambda, dt_lev);

        //
        //-----------------------------------------------------------------------
        // prestep_uv_3d
        //-----------------------------------------------------------------------
        //
        //updates u,v,ru,rv (ru and rv have multiple components)
        prestep_uv_3d(bx, gbx, uold, vold, u, v, ru, rv, Hz, Akv,
                          pm, pn, W, DC, FC, z_r, sustr, svstr, bustr, bvstr,
                          iic, ntfirst, nnew, nstp, nrhs, N, lambda, dt_lev);
    }
}
