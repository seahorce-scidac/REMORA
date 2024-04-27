#include <REMORA.H>

using namespace amrex;

void
REMORA::gls_prestep (int lev, MultiFab* mf_gls, MultiFab* mf_tke,
                     MultiFab& mf_W, const int nstp, const int nnew,
                     const int iic, const int ntfirst, const int N, const Real dt_lev)
{
    // temps: grad, gradL, XF, FX, FXL, EF, FE, FEL
    for ( MFIter mfi(*mf_gls, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<Real> const& gls = mf_gls->array(mfi);
        Array4<Real> const& tke = mf_tke->array(mfi);
        Array4<Real> const& W = mf_W.array(mfi);

        Array4<Real> const& Huon = vec_Huon[lev]->array(mfi);
        Array4<Real> const& Hvom = vec_Hvom[lev]->array(mfi);
        Array4<Real> const& Hz = vec_Hz[lev]->array(mfi);
        Array4<Real> const& pm = vec_pm[lev]->array(mfi);
        Array4<Real> const& pn = vec_pn[lev]->array(mfi);

        Box xbx = mfi.nodaltilebox(0);
        Box ybx = mfi.nodaltilebox(1);
        Box bx = mfi.tilebox();

        Box xbx_hi = growHi(xbx,0,1);

        Box ybx_hi = growHi(ybx,0,1);

        const Box& domain = geom[0].Domain();
        const auto dlo = amrex::lbound(domain);
        const auto dhi = amrex::ubound(domain);

        int ncomp = 1;
        Vector<BCRec> bcrs_x(ncomp);
        Vector<BCRec> bcrs_y(ncomp);
        amrex::setBC(xbx,domain,BCVars::xvel_bc,0,1,domain_bcs_type,bcrs_x);
        amrex::setBC(ybx,domain,BCVars::yvel_bc,0,1,domain_bcs_type,bcrs_y);
        auto bcr_x = bcrs_x[0];
        auto bcr_y = bcrs_y[0];

        FArrayBox fab_XF(xbx_hi, 1, amrex::The_Async_Arena()); fab_XF.template setVal<RunOn::Device>(0.);
        FArrayBox fab_FX(xbx_hi, 1, amrex::The_Async_Arena()); fab_FX.template setVal<RunOn::Device>(0.);
        FArrayBox fab_FXL(xbx_hi, 1, amrex::The_Async_Arena()); fab_FXL.template setVal<RunOn::Device>(0.);
        FArrayBox fab_EF(ybx_hi, 1, amrex::The_Async_Arena()); fab_EF.template setVal<RunOn::Device>(0.);
        FArrayBox fab_FE(ybx_hi, 1, amrex::The_Async_Arena()); fab_FE.template setVal<RunOn::Device>(0.);
        FArrayBox fab_FEL(ybx_hi, 1, amrex::The_Async_Arena()); fab_FEL.template setVal<RunOn::Device>(0.);
        FArrayBox fab_Hz_half(bx, 1, amrex::The_Async_Arena()); fab_Hz_half.template setVal<RunOn::Device>(0.);
        FArrayBox fab_CF(bx, 1, amrex::The_Async_Arena()); fab_CF.template setVal<RunOn::Device>(0.);
        FArrayBox fab_FC(bx, 1, amrex::The_Async_Arena()); fab_FC.template setVal<RunOn::Device>(0.);
        FArrayBox fab_FCL(bx, 1, amrex::The_Async_Arena()); fab_FCL.template setVal<RunOn::Device>(0.);

        auto XF  = fab_XF.array();
        auto FX  = fab_FX.array();
        auto FXL = fab_FXL.array();
        auto EF  = fab_EF.array();
        auto FE  = fab_FE.array();
        auto FEL = fab_FEL.array();
        auto Hz_half = fab_Hz_half.array();
        auto CF  = fab_CF.array();
        auto FC  = fab_FC.array();
        auto FCL = fab_FCL.array();

        auto ic_bc_type = solverChoice.ic_bc_type;

        // need XF/FX/FXL from  [xlo to xhi+1] by [ylo to yhi  ]
        ParallelFor(growHi(xbx,0,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real grad_im1 = tke(i-1,j,k,nstp) - tke(i-2,j,k,nstp);
            Real grad_ip1 = tke(i+1,j,k,nstp) - tke(i  ,j,k,nstp);

            Real gradL_im1 = gls(i-1,j,k,nstp) - gls(i-2,j,k,nstp);
            Real gradL_ip1 = gls(i+1,j,k,nstp) - gls(i  ,j,k,nstp);

            // Adjust boundaries
            // TODO: Make sure indices match with what ROMS does
            if (i == dlo.x-1 && (bcr_x.lo(0) == REMORABCType::ext_dir or ic_bc_type==IC_BC_Type::Real)) {
                grad_im1  = tke(i,j,k,nstp) - tke(i-1,j,k,nstp);
                gradL_im1 = gls(i,j,k,nstp) - gls(i-1,j,k,nstp);
            }
            else if (i == dhi.x+1 && (bcr_x.hi(0) == REMORABCType::ext_dir or ic_bc_type==IC_BC_Type::Real)) {
                grad_ip1  = tke(i,j,k,nstp) - tke(i-1,j,k,nstp);
                gradL_ip1 = gls(i,j,k,nstp) - gls(i-1,j,k,nstp);
            }
            Real cff = 1.0_rt/6.0_rt;
            XF(i,j,k) = 0.5_rt * (Huon(i,j,k) + Huon(i,j,k+1));
            FX(i,j,k) = XF(i,j,k) * 0.5_rt * (tke(i-1,j,k,nstp) + tke(i,j,k,nstp) -
                cff * (grad_ip1 - grad_im1));
            FXL(i,j,k) = XF(i,j,k) * 0.5_rt * (gls(i-1,j,k,nstp) + gls(i,j,k,nstp) -
                cff * (gradL_ip1 - gradL_im1));
        });

        // need EF/FE/FEL from  [xlo to xhi  ] by [ylo to yhi+1]
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real grad_jm1 = tke(i,j-1,k,nstp) - tke(i,j-2,k,nstp);
            Real grad_jp1 = tke(i,j+1,k,nstp) - tke(i,j  ,k,nstp);

            Real gradL_jm1 = gls(i,j-1,k,nstp) - gls(i,j-2,k,nstp);
            Real gradL_jp1 = gls(i,j+1,k,nstp) - gls(i,j  ,k,nstp);

            // Adjust boundaries
            // TODO: Make sure indices match with what ROMS does
            if (j == dlo.y-1 && (bcr_y.lo(1) == REMORABCType::ext_dir or ic_bc_type==IC_BC_Type::Real)) {
                grad_jm1  = tke(i,j,k,nstp) - tke(i,j-1,k,nstp);
                gradL_jm1 = gls(i,j,k,nstp) - gls(i,j-1,k,nstp);
            }
            else if (j == dhi.y+1 && (bcr_y.hi(1) == REMORABCType::ext_dir or ic_bc_type==IC_BC_Type::Real)) {
                grad_jp1  = tke(i,j,k,nstp) - tke(i,j-1,k,nstp);
                gradL_jp1 = gls(i,j,k,nstp) - gls(i,j-1,k,nstp);
            }
            Real cff = 1.0_rt/6.0_rt;
            EF(i,j,k) = 0.5_rt * (Hvom(i,j,k) + Hvom(i,j,k+1));
            FE(i,j,k) = EF(i,j,k) * 0.5_rt * (tke(i,j-1,k,nstp) + tke(i,j,k,nstp) -
                cff * (grad_jp1 - grad_jm1));
            FEL(i,j,k) = EF(i,j,k) * 0.5_rt * (gls(i,j-1,k,nstp) + gls(i,j,k,nstp) -
                cff * (gradL_jp1 - gradL_jm1));
        });

        Real gamma = 1.0_rt / 6.0_rt;
        Real cff1, cff2, cff3;
        int indx;
        // Time step horizontal advection
        if (iic == ntfirst) {
            cff1 = 1.0_rt;
            cff2 = 0.0_rt;
            cff3 = 0.5_rt * dt_lev;
            indx = nstp;
        } else {
            cff1 = 0.5_rt + gamma;
            cff2 = 0.5_rt - gamma;
            cff3 = (1.0_rt - gamma) * dt_lev;
            indx = 3 - nstp;
        }

        // update tke, gls from [xlo to xhi  ] by [ylo to yhi  ]
        // need XF/FX/FXL from  [xlo to xhi+1] by [ylo to yhi  ]
        // need EF/FE/FEL from  [xlo to xhi  ] by [ylo to yhi+1]
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff = 0.5_rt * (Hz(i,j,k) + Hz(i,j,k+1));
            Real cff4 = cff3 * pm(i,j,0) * pn(i,j,0);
            Hz_half(i,j,k) = cff - cff4 * (XF(i+1,j,k)-XF(i,j,k)+EF(i,j+1,k)-EF(i,j,k));
            tke(i,j,k,2) = cff * (cff1*tke(i,j,k,nstp) + cff2*tke(i,j,k,indx)) -
                           cff4 * (FX(i+1,j,k)-FX(i,j,k)+FE(i,j+1,k)-FE(i,j,k));
            gls(i,j,k,2) = cff * (cff1 * gls(i,j,k,nstp) + cff2 * gls(i,j,k,indx)) -
                           cff4 * (FXL(i+1,j,k)-FXL(i,j,k)+FEL(i,j+1,k)-FEL(i,j,k));
            tke(i,j,k,nnew) = cff * tke(i,j,k,nstp);
            gls(i,j,k,nnew) = cff * gls(i,j,k,nstp);
        });

        // Will do a FillPatch after this, so don't need to do any ghost zones in x,y
        // Compute vertical advection
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            CF(i,j,k) = 0.5_rt * (W(i,j,k) + W(i,j,k-1));
            // DO k=2,N(ng)-1
            if (k == 0) {
                Real cff1 = 1.0_rt / 3.0_rt;
                Real cff2 = 5.0_rt / 6.0_rt;
                Real cff3 = 1.0_rt / 6.0_rt;
                FC(i,j,k)  = CF(i,j,k) * (cff1 * tke(i,j,0,nstp) +
                                          cff2 * tke(i,j,1,nstp) -
                                          cff3 * tke(i,j,2,nstp));
                FCL(i,j,k) = CF(i,j,k) * (cff1 * gls(i,j,0,nstp) +
                                          cff2 * gls(i,j,1,nstp) -
                                          cff3 * gls(i,j,2,nstp));
            } else if (k == N) {
                Real cff1 = 1.0_rt / 3.0_rt;
                Real cff2 = 5.0_rt / 6.0_rt;
                Real cff3 = 1.0_rt / 6.0_rt;
                FC(i,j,k)  = CF(i,j,k) * (cff1 * tke(i,j,k,  nstp) +
                                          cff2 * tke(i,j,k-1,nstp)-
                                          cff3 * tke(i,j,k-2,nstp));
                FCL(i,j,k) = CF(i,j,k) * (cff1 * gls(i,j,k,  nstp) +
                                          cff2 * gls(i,j,k-1,nstp)-
                                          cff3 * gls(i,j,k-2,nstp));
            } else {
                Real cff1 = 7.0_rt / 12.0_rt;
                Real cff2 = 1.0_rt / 12.0_rt;
                FC(i,j,k)  = CF(i,j,k) * (cff1 * (tke(i,j,k-1,nstp) + tke(i,j,k  ,nstp)) -
                                          cff2 * (tke(i,j,k-2,nstp) + tke(i,j,k+1,nstp)));
                FCL(i,j,k) = CF(i,j,k) * (cff1 * (gls(i,j,k-1,nstp) + gls(i,j,k  ,nstp)) -
                                          cff2 * (gls(i,j,k-2,nstp) + gls(i,j,k+1,nstp)));
            }
        });

        // Time-step vertical advection
        if (iic == ntfirst) {
            cff3 = 0.5_rt * dt_lev;
        } else {
            cff3 = (1.0_rt - gamma) * dt_lev;
        }
        // DO k=1,N-1
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff4 = cff3 * pm(i,j,k) * pn(i,j,k);
            Hz_half(i,j,k) = Hz_half(i,j,k) - cff4 * (CF(i,j,k+1)-CF(i,j,k));
            Real cff1 = 1.0_rt / Hz_half(i,j,k);
            tke(i,j,k,2) = cff1 * (tke(i,j,k,2) - cff4 * (FC (i,j,k+1) - FC (i,j,k)));
            gls(i,j,k,2) = cff1 * (gls(i,j,k,2) - cff4 * (FCL(i,j,k+1) - FCL(i,j,k)));
        });
    }

    FillPatch(lev, t_old[lev], *vec_tke[lev], GetVecOfPtrs(vec_tke), BdyVars::t);
    FillPatch(lev, t_old[lev], *vec_gls[lev], GetVecOfPtrs(vec_gls), BdyVars::t);
}

void
REMORA::gls_corrector(int lev, )
{
}
