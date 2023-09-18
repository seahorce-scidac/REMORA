#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

void
ROMSX::prestep (int lev,
                MultiFab& mf_uold, MultiFab& mf_vold,
                MultiFab& mf_u, MultiFab& mf_v,
                std::unique_ptr<MultiFab>& mf_ru,
                std::unique_ptr<MultiFab>& mf_rv,
                MultiFab& mf_tempold, MultiFab& mf_saltold,
                MultiFab& mf_temp, MultiFab& mf_salt,
                std::unique_ptr<MultiFab>& /* mf_Hz */,
                std::unique_ptr<MultiFab>& /* mf_Akv */,
                std::unique_ptr<MultiFab>& /* mf_Huon */,
                std::unique_ptr<MultiFab>& /* mf_Hvom */,
                MultiFab& mf_W, MultiFab& mf_DC,
                /* MF mf_FC? */
                std::unique_ptr<MultiFab>& /* mf_t3 */,
                std::unique_ptr<MultiFab>& /* mf_s3 */,
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

    const BoxArray&            ba = mf_temp.boxArray();
    const DistributionMapping& dm = mf_temp.DistributionMap();

    // Maybe not the best way to do this, but need to cache salt and temp since
    // they get rewritten by prestep_t
    MultiFab mf_saltcache(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_tempcache(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab::Copy(mf_saltcache,mf_salt,0,0,mf_salt.nComp(),IntVect(NGROW,NGROW,0)); //mf_salt.nGrowVect());
    MultiFab::Copy(mf_tempcache,mf_temp,0,0,mf_temp.nComp(),IntVect(NGROW,NGROW,0));

    for ( MFIter mfi(mf_temp, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Array4<Real> const& DC = mf_DC.array(mfi);
        Array4<Real> const& Akv = (vec_Akv[lev])->array(mfi);
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
        Array4<Real> const& tempold = (mf_tempold).array(mfi);
        Array4<Real> const& saltold = (mf_saltold).array(mfi);
        Array4<Real> const& temp = (mf_temp).array(mfi);
        Array4<Real> const& salt = (mf_salt).array(mfi);
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
        FArrayBox fab_FX(gbx2,1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_FE(gbx2,1,amrex::The_Async_Arena()); //3D
        FArrayBox fab_pn(tbxp2D,1,amrex::The_Async_Arena());
        FArrayBox fab_pm(tbxp2D,1,amrex::The_Async_Arena());
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
        FArrayBox fab_fomn(tbxp2D,1,amrex::The_Async_Arena());
        //FArrayBox fab_oHz(gbx11,1,amrex::The_Async_Arena());

        auto FC=fab_FC.array();
        auto FX=fab_FX.array();
        auto FE=fab_FE.array();
        auto pn=fab_pn.array();
        auto pm=fab_pm.array();
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
        auto fomn=fab_fomn.array();

        amrex::ParallelFor(tbxp2D,
        [=] AMREX_GPU_DEVICE (int i, int j, int  )
            {
              pm(i,j,0)=dxi[0];
              pn(i,j,0)=dxi[1];
              //defined UPWELLING
              Real f0=-8.26e-5;
              Real beta=0.0;
              Real Esize=1000*(Mm);
              Real y = prob_lo[1] + (j + 0.5) * dx[1];
              Real f=fomn(i,j,0)=f0+beta*(y-.5*Esize);
              fomn(i,j,0)=f*(1.0/(pm(i,j,0)*pn(i,j,0)));
            });

        amrex::ParallelFor(tbxp2D,
        [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
          //Note: are the comment definitons right? Don't seem to match metrics.f90
          om_v(i,j,0)=1.0/dxi[0]; // 2/(pm(i,j-1)+pm(i,j))
          on_u(i,j,0)=1.0/dxi[1]; // 2/(pm(i,j-1)+pm(i,j))
          om_r(i,j,0)=1.0/dxi[0]; // 1/pm(i,j)
          on_r(i,j,0)=1.0/dxi[1]; // 1/pn(i,j)
          //todo: om_p on_p
          om_p(i,j,0)=1.0/dxi[0]; // 4/(pm(i-1,j-1)+pm(i-1,j)+pm(i,j-1)+pm(i,j))
          on_p(i,j,0)=1.0/dxi[1]; // 4/(pn(i-1,j-1)+pn(i-1,j)+pn(i,j-1)+pn(i,j))
          on_v(i,j,0)=1.0/dxi[1]; // 2/(pn(i-1,j)+pn(i,j))
          om_u(i,j,0)=1.0/dxi[0]; // 2/(pm(i-1,j)+pm(i,j))
          pmon_u(i,j,0)=1.0;        // (pm(i-1,j)+pm(i,j))/(pn(i-1,j)+pn(i,j))
          pnom_u(i,j,0)=1.0;        // (pn(i-1,j)+pn(i,j))/(pm(i-1,j)+pm(i,j))
          pmon_v(i,j,0)=1.0;        // (pm(i,j-1)+pm(i,j))/(pn(i,j-1)+pn(i,j))
          pnom_v(i,j,0)=1.0;        // (pn(i,j-1)+pn(i,j))/(pm(i,j-1)+pm(i,j))
        });
        if (verbose > 2) {
           amrex::PrintToFile("temp_preprestep").SetPrecision(18)<<FArrayBox(temp)<<std::endl;
           amrex::PrintToFile("tempstore_preprestep").SetPrecision(18)<<FArrayBox(tempstore)<<std::endl;
           amrex::PrintToFile("salt_preprestep").SetPrecision(18)<<FArrayBox(salt)<<std::endl;
           amrex::PrintToFile("saltstore_preprestep").SetPrecision(18)<<FArrayBox(saltstore)<<std::endl;
           amrex::PrintToFile("saltold_preprestep").SetPrecision(18)<<FArrayBox(saltold)<<std::endl;
           amrex::PrintToFile("tempold_preprestep").SetPrecision(18)<<FArrayBox(tempold)<<std::endl;
        }

        if (verbose > 1) {
            Print() << "Akv box " << Box(Akv) << std::endl;
        }
        prestep_t_3d(bx, gbx, uold, vold, u, v, tempold, saltold, temp, salt, tempcache,ru, rv, Hz, Akv, on_u, om_v, Huon, Hvom,
                     pm, pn, W, DC, FC, tempstore, saltstore, FX, FE, z_r, z_w, h, iic, ntfirst, nnew, nstp, nrhs, N,
                          lambda, dt_lev);
        prestep_t_3d(bx, gbx, uold, vold, u, v, saltold, saltold, salt, salt, saltcache, ru, rv, Hz, Akv, on_u, om_v, Huon, Hvom,
                     pm, pn, W, DC, FC, saltstore, saltstore, FX, FE, z_r, z_w, h, iic, ntfirst, nnew, nstp, nrhs, N,
                          lambda, dt_lev);

       if (verbose > 2) {
           amrex::PrintToFile("u").SetPrecision(18)<<FArrayBox(u)<<std::endl;
           amrex::PrintToFile("u").SetPrecision(18)<<FArrayBox(u)<<std::endl;
           amrex::PrintToFile("v").SetPrecision(18)<<FArrayBox(v)<<std::endl;
           amrex::PrintToFile("uold").SetPrecision(18)<<FArrayBox(uold)<<std::endl;
           amrex::PrintToFile("vold").SetPrecision(18)<<FArrayBox(vold)<<std::endl;
           amrex::PrintToFile("temp").SetPrecision(18)<<FArrayBox(temp)<<std::endl;
           amrex::PrintToFile("tempstore").SetPrecision(18)<<FArrayBox(tempstore)<<std::endl;
           amrex::PrintToFile("salt").SetPrecision(18)<<FArrayBox(salt)<<std::endl;
           amrex::PrintToFile("saltstore").SetPrecision(18)<<FArrayBox(saltstore)<<std::endl;
       }
        //
        //-----------------------------------------------------------------------
        // prestep_uv_3d
        //-----------------------------------------------------------------------
        //
        //updates u,v,ru,rv (ru and rv have multiple components)
        prestep_uv_3d(bx, gbx, uold, vold, u, v, ru, rv, Hz, Akv, on_u, om_v, Huon, Hvom,
                          pm, pn, W, DC, FC, z_r, sustr, svstr, bustr, bvstr, iic, ntfirst, nnew, nstp, nrhs, N,
                          lambda, dt_lev);
       if (verbose > 2) {
           amrex::PrintToFile("u_after_prestep").SetPrecision(18)<<FArrayBox(u)<<std::endl;
           amrex::PrintToFile("v_after_prestep").SetPrecision(18)<<FArrayBox(v)<<std::endl;
           amrex::PrintToFile("ru_after_prestep").SetPrecision(18)<<FArrayBox(ru)<<std::endl;
           amrex::PrintToFile("rv_after_prestep").SetPrecision(18)<<FArrayBox(rv)<<std::endl;
           amrex::PrintToFile("temp_afterprestep").SetPrecision(18)<<FArrayBox(temp)<<std::endl;
           amrex::PrintToFile("tempstore_afterprestep").SetPrecision(18)<<FArrayBox(tempstore)<<std::endl;
           amrex::PrintToFile("salt_afterprestep").SetPrecision(18)<<FArrayBox(salt)<<std::endl;
           amrex::PrintToFile("saltstore_afterprestep").SetPrecision(18)<<FArrayBox(saltstore)<<std::endl;
           amrex::PrintToFile("saltold_afterprestep").SetPrecision(18)<<FArrayBox(saltold)<<std::endl;
           amrex::PrintToFile("tempold_afterprestep").SetPrecision(18)<<FArrayBox(tempold)<<std::endl;
       }
    }
}
