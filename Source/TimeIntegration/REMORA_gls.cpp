#include <REMORA.H>

using namespace amrex;

void
REMORA::gls_prestep (int lev, MultiFab* mf_gls, MultiFab* mf_tke,
                     MultiFab& mf_W, MultiFab* mf_msku, MultiFab* mf_mskv,
                     const int nstp, const int nnew,
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
        Array4<Real> const& msku = vec_msku[lev]->array(mfi);
        Array4<Real> const& mskv = vec_mskv[lev]->array(mfi);

        Box bx = mfi.tilebox();
        Box xbx = surroundingNodes(bx,0);
        Box ybx = surroundingNodes(bx,1);

        Box xbx_hi = growHi(xbx,0,1);

        Box ybx_hi = growHi(ybx,0,1);

        const Box& domain = geom[0].Domain();
        const auto dlo = amrex::lbound(domain);
        const auto dhi = amrex::ubound(domain);

        GeometryData const& geomdata = geom[0].data();
        bool is_periodic_in_x = geomdata.isPeriodic(0);
        bool is_periodic_in_y = geomdata.isPeriodic(1);

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
        FArrayBox fab_CF(convert(bx,IntVect(0,0,0)), 1, amrex::The_Async_Arena()); fab_CF.template setVal<RunOn::Device>(0.);
        FArrayBox fab_FC(convert(bx,IntVect(0,0,0)), 1, amrex::The_Async_Arena()); fab_FC.template setVal<RunOn::Device>(0.);
        FArrayBox fab_FCL(convert(bx,IntVect(0,0,0)), 1, amrex::The_Async_Arena()); fab_FCL.template setVal<RunOn::Device>(0.);

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

        // need XF/FX/FXL from  [xlo to xhi] by [ylo to yhi  ] on u points
        ParallelFor(grow(xbx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real grad_im1 = (tke(i-1,j,k,nstp) - tke(i-2,j,k,nstp)) * msku(i-1,j,0);
            Real grad_ip1 = (tke(i+1,j,k,nstp) - tke(i  ,j,k,nstp)) * msku(i+1,j,0);

            Real gradL_im1 = (gls(i-1,j,k,nstp) - gls(i-2,j,k,nstp)) * msku(i-1,j,0);
            Real gradL_ip1 = (gls(i+1,j,k,nstp) - gls(i  ,j,k,nstp)) * msku(i+1,j,0);

            // Adjust boundaries
            // TODO: Make sure indices match with what ROMS does
            if (i == dlo.x-1 && !is_periodic_in_x) {
                grad_im1  = tke(i,j,k,nstp) - tke(i-1,j,k,nstp);
                gradL_im1 = gls(i,j,k,nstp) - gls(i-1,j,k,nstp);
            }
            else if (i == dhi.x+1 && !is_periodic_in_x) {
                grad_ip1  = tke(i,j,k,nstp) - tke(i-1,j,k,nstp);
                gradL_ip1 = gls(i,j,k,nstp) - gls(i-1,j,k,nstp);
            }
            Real cff = 1.0_rt/6.0_rt;
            XF(i,j,k) = 0.5_rt * (Huon(i,j,k) + Huon(i,j,k-1));
            FX(i,j,k) = XF(i,j,k) * 0.5_rt * (tke(i-1,j,k,nstp) + tke(i,j,k,nstp) -
                cff * (grad_ip1 - grad_im1));
            FXL(i,j,k) = XF(i,j,k) * 0.5_rt * (gls(i-1,j,k,nstp) + gls(i,j,k,nstp) -
                cff * (gradL_ip1 - gradL_im1));
        });

        // need EF/FE/FEL from  [xlo to xhi  ] by [ylo to yhi+1]
        ParallelFor(grow(ybx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real grad_jm1 = (tke(i,j-1,k,nstp) - tke(i,j-2,k,nstp)) * mskv(i,j-1,0);
            Real grad_jp1 = (tke(i,j+1,k,nstp) - tke(i,j  ,k,nstp)) * mskv(i,j+1,0);

            Real gradL_jm1 = (gls(i,j-1,k,nstp) - gls(i,j-2,k,nstp)) * mskv(i,j-1,0);
            Real gradL_jp1 = (gls(i,j+1,k,nstp) - gls(i,j  ,k,nstp)) * mskv(i,j+1,0);

            // Adjust boundaries
            // TODO: Make sure indices match with what ROMS does
            if (j == dlo.y-1 && !is_periodic_in_y) {
                grad_jm1  = tke(i,j,k,nstp) - tke(i,j-1,k,nstp);
                gradL_jm1 = gls(i,j,k,nstp) - gls(i,j-1,k,nstp);
            }
            else if (j == dhi.y+1 && !is_periodic_in_y) {
                grad_jp1  = tke(i,j,k,nstp) - tke(i,j-1,k,nstp);
                gradL_jp1 = gls(i,j,k,nstp) - gls(i,j-1,k,nstp);
            }
            Real cff = 1.0_rt/6.0_rt;
            EF(i,j,k) = 0.5_rt * (Hvom(i,j,k) + Hvom(i,j,k-1));
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
            indx = 1 - nstp;
        }

        // update tke, gls from [xlo to xhi  ] by [ylo to yhi  ]
        // need XF/FX/FXL from  [xlo to xhi+1] by [ylo to yhi  ]
        // need EF/FE/FEL from  [xlo to xhi  ] by [ylo to yhi+1]
        ParallelFor(grow(bx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff = 0.5_rt * (Hz(i,j,k) + Hz(i,j,k-1));
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
        ParallelFor(convert(bx,IntVect(0,0,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // CF and FC/FCL are on rho points
            CF(i,j,k) = 0.5_rt * (W(i,j,k+1) + W(i,j,k));
            if (k == 0) {
                Real cff1_vadv = 1.0_rt / 3.0_rt;
                Real cff2_vadv = 5.0_rt / 6.0_rt;
                Real cff3_vadv = 1.0_rt / 6.0_rt;
                FC(i,j,k)  = CF(i,j,k) * (cff1_vadv * tke(i,j,0,nstp) +
                                          cff2_vadv * tke(i,j,1,nstp) -
                                          cff3_vadv * tke(i,j,2,nstp));
                FCL(i,j,k) = CF(i,j,k) * (cff1_vadv * gls(i,j,0,nstp) +
                                          cff2_vadv * gls(i,j,1,nstp) -
                                          cff3_vadv * gls(i,j,2,nstp));
            } else if (k == N) {
                Real cff1_vadv = 1.0_rt / 3.0_rt;
                Real cff2_vadv = 5.0_rt / 6.0_rt;
                Real cff3_vadv = 1.0_rt / 6.0_rt;
                FC(i,j,k)  = CF(i,j,k) * (cff1_vadv * tke(i,j,k+1,  nstp) +
                                          cff2_vadv * tke(i,j,k  ,nstp)-
                                          cff3_vadv * tke(i,j,k-1,nstp));
                FCL(i,j,k) = CF(i,j,k) * (cff1_vadv * gls(i,j,k+1,nstp) +
                                          cff2_vadv * gls(i,j,k  ,nstp)-
                                          cff3_vadv * gls(i,j,k-1,nstp));
            } else {
                Real cff1_vadv = 7.0_rt / 12.0_rt;
                Real cff2_vadv = 1.0_rt / 12.0_rt;
                FC(i,j,k)  = CF(i,j,k) * (cff1_vadv * (tke(i,j,k  ,nstp) + tke(i,j,k+1,nstp)) -
                                          cff2_vadv * (tke(i,j,k-1,nstp) + tke(i,j,k+2,nstp)));
                FCL(i,j,k) = CF(i,j,k) * (cff1_vadv * (gls(i,j,k  ,nstp) + gls(i,j,k+1,nstp)) -
                                          cff2_vadv * (gls(i,j,k-1,nstp) + gls(i,j,k+2,nstp)));
            }
        });

        // Time-step vertical advection
        if (iic == ntfirst) {
            cff3 = 0.5_rt * dt_lev;
        } else {
            cff3 = (1.0_rt - gamma) * dt_lev;
        }
        // DO k=1,N-1
        ParallelFor(grow(bx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff4 = cff3 * pm(i,j,0) * pn(i,j,0);
            Hz_half(i,j,k) = Hz_half(i,j,k) - cff4 * (CF(i,j,k)-CF(i,j,k-1));
            Real cff1_loc = 1.0_rt / Hz_half(i,j,k);
            tke(i,j,k,2) = cff1_loc * (tke(i,j,k,2) - cff4 * (FC (i,j,k) - FC (i,j,k-1)));
            gls(i,j,k,2) = cff1_loc * (gls(i,j,k,2) - cff4 * (FCL(i,j,k) - FCL(i,j,k-1)));
        });
    }

    for (int icomp=0; icomp<3; icomp++) {
        FillPatch(lev, t_old[lev], *vec_tke[lev], GetVecOfPtrs(vec_tke), BCVars::zvel_bc, BdyVars::null, icomp, false, false);
        FillPatch(lev, t_old[lev], *vec_gls[lev], GetVecOfPtrs(vec_gls), BCVars::zvel_bc, BdyVars::null, icomp, false, false);
    }
}

void
REMORA::gls_corrector (int lev, MultiFab* mf_gls, MultiFab* mf_tke,
                       MultiFab& mf_W, MultiFab* mf_Akv, MultiFab* mf_Akt,
                       MultiFab* mf_Akk, MultiFab* mf_Akp,
                       MultiFab* mf_mskr,
                       MultiFab* mf_msku, MultiFab* mf_mskv,
                       const int nstp, const int nnew,
                       const int N, const Real dt_lev)
{
//-----------------------------------------------------------------------
//  Compute several constants.
//-----------------------------------------------------------------------
    bool Lmy25 = ((solverChoice.gls_p == 0.0) &&
                  (solverChoice.gls_n == 1.0) &&
                  (solverChoice.gls_m == 1.0)) ? true : false;

    Real L_sft = solverChoice.vonKar;
    Real gls_sigp_cb = solverChoice.gls_sigp;
    Real ogls_sigp = 1.0_rt/gls_sigp_cb;

    Real gls_c3m = solverChoice.gls_c3m;
    Real gls_c3p = solverChoice.gls_c3p;
    Real gls_cmu0 = solverChoice.gls_cmu0;

    Real gls_m = solverChoice.gls_m;
    Real gls_n = solverChoice.gls_n;
    Real gls_p = solverChoice.gls_p;

    Real gls_Gh0 = solverChoice.gls_Gh0;
    Real gls_Ghcri = solverChoice.gls_Ghcri;
    Real gls_Ghmin = solverChoice.gls_Ghmin;

    Real Akv_bak = solverChoice.Akv_bak;
    Real Akt_bak = solverChoice.Akt_bak;
    Real Akp_bak = solverChoice.Akp_bak;
    Real Akk_bak = solverChoice.Akk_bak;

    Real gls_c1 = solverChoice.gls_c1;
    Real gls_c2 = solverChoice.gls_c2;
    Real gls_E2 = solverChoice.gls_E2;
    Real gls_sigk = solverChoice.gls_sigk;
    Real vonKar = solverChoice.vonKar;
    auto gls_stability_type = solverChoice.gls_stability_type;

    Real sqrt2 = std::sqrt(2.0_rt);
    Real cmu_fac1 = std::pow(solverChoice.gls_cmu0,(-solverChoice.gls_p/solverChoice.gls_n));
    Real cmu_fac2 = std::pow(solverChoice.gls_cmu0,(3.0_rt+solverChoice.gls_p/solverChoice.gls_n));
    Real cmu_fac3 = 1.0_rt/std::pow(solverChoice.gls_cmu0,2.0_rt);
    //Real cmu_fac4 = std::pow(1.5_rt*solverChoice.gls_sigk,(1.0_rt/3.0_rt))/std::pow(solverChoice.gls_cmu0,4.0_rt/3.0_rt);

    //Real gls_fac1 = solverChoice.gls_n*std::pow(solverChoice.gls_cmu0,solverChoice.gls_p+1.0_rt);
    Real gls_fac2 = std::pow(solverChoice.gls_cmu0,solverChoice.gls_p)*solverChoice.gls_n*std::pow(solverChoice.vonKar,solverChoice.gls_n);
    Real gls_fac3 = std::pow(solverChoice.gls_cmu0,solverChoice.gls_p)*solverChoice.gls_n;
    Real gls_fac4 = std::pow(solverChoice.gls_cmu0,solverChoice.gls_p);
    Real gls_fac5 = std::pow(0.56_rt,0.5_rt*solverChoice.gls_n)*std::pow(solverChoice.gls_cmu0,solverChoice.gls_p);
    Real gls_fac6 = 8.0_rt/std::pow(solverChoice.gls_cmu0,6.0_rt);

    Real gls_exp1 = 1.0_rt/solverChoice.gls_n;
    Real tke_exp1 = solverChoice.gls_m/solverChoice.gls_n;
    Real tke_exp2 = 0.5_rt+solverChoice.gls_m/solverChoice.gls_n;
    //Real tke_exp3 = 0.5_rt+solverChoice.gls_m;
    Real tke_exp4 = solverChoice.gls_m+0.5_rt*solverChoice.gls_n;

    Real gls_s0, gls_s1, gls_s2, gls_s4, gls_s5, gls_s6;
    Real gls_b0, gls_b1, gls_b2, gls_b3, gls_b4, gls_b5;
    Real my_Sm2, my_Sh1, my_Sh2, my_Sm3, my_Sm4;

    // Compute parameters for Canuto et al. (2001) stability functions.
    // (Canuto, V.M., Cheng, H.Y., and Dubovikov, M.S., 2001: Ocean
    // turbulence. Part I: One-point closure model - momentum and
    // heat vertical diffusivities, JPO, 1413-1426).

    if (solverChoice.gls_stability_type == GLS_StabilityType::Canuto_A ||
        solverChoice.gls_stability_type == GLS_StabilityType::Canuto_B) {

        gls_s0=3.0_rt/2.0_rt*solverChoice.gls_L1*solverChoice.gls_L5*solverChoice.gls_L5;
        gls_s1=-solverChoice.gls_L4*(solverChoice.gls_L6+solverChoice.gls_L7)
                        +2.0_rt*solverChoice.gls_L4*solverChoice.gls_L5*
                        (solverChoice.gls_L1-1.0_rt/3.0_rt*solverChoice.gls_L2-solverChoice.gls_L3)
                        +3.0_rt/2.0_rt*
                        solverChoice.gls_L1*solverChoice.gls_L5*solverChoice.gls_L8;
        gls_s2=-3.0_rt/8.0_rt*solverChoice.gls_L1
            *(solverChoice.gls_L6*solverChoice.gls_L6-solverChoice.gls_L7*solverChoice.gls_L7);
        gls_s4=2.0_rt*solverChoice.gls_L5;
        gls_s5=2.0_rt*solverChoice.gls_L4;
        gls_s6=2.0_rt/3.0_rt*solverChoice.gls_L5
            *(3.0_rt*solverChoice.gls_L3*solverChoice.gls_L3-solverChoice.gls_L2*solverChoice.gls_L2)-
                    1.0_rt/2.0_rt*solverChoice.gls_L5*solverChoice.gls_L1*(3.0_rt*solverChoice.gls_L3-solverChoice.gls_L2)+
                    3.0_rt/4.0_rt*solverChoice.gls_L1*(solverChoice.gls_L6-solverChoice.gls_L7);
        gls_b0=3.0_rt*solverChoice.gls_L5*solverChoice.gls_L5;
        gls_b1=solverChoice.gls_L5*(7.0_rt*solverChoice.gls_L4+3.0_rt*solverChoice.gls_L8);
        gls_b2=solverChoice.gls_L5*solverChoice.gls_L5*(3.0_rt*solverChoice.gls_L3*solverChoice.gls_L3-solverChoice.gls_L2*solverChoice.gls_L2)-
                    3.0_rt/4.0_rt*(solverChoice.gls_L6*solverChoice.gls_L6-solverChoice.gls_L7*solverChoice.gls_L7);
        gls_b3=solverChoice.gls_L4*(4.0_rt*solverChoice.gls_L4+3.0_rt*solverChoice.gls_L8);
        gls_b5=1.0_rt/4.0_rt*(solverChoice.gls_L2*solverChoice.gls_L2-3.0_rt*solverChoice.gls_L3*solverChoice.gls_L3)*
                    (solverChoice.gls_L6*solverChoice.gls_L6-solverChoice.gls_L7*solverChoice.gls_L7);
        gls_b4=solverChoice.gls_L4*(solverChoice.gls_L2*solverChoice.gls_L6-3.0_rt*solverChoice.gls_L3*solverChoice.gls_L7-
                    solverChoice.gls_L5*(solverChoice.gls_L2*solverChoice.gls_L2-solverChoice.gls_L3*solverChoice.gls_L3))+solverChoice.gls_L5*solverChoice.gls_L8*
                    (3.0_rt*solverChoice.gls_L3*solverChoice.gls_L3-solverChoice.gls_L2*solverChoice.gls_L2);
        my_Sm2 = 0.0_rt;
        my_Sm3 = 0.0_rt;
        my_Sm4 = 0.0_rt;
        my_Sh1 = 0.0_rt;
        my_Sh2 = 0.0_rt;
    } else {
        gls_s0 = 0.0_rt;
        gls_s1 = 0.0_rt;
        gls_s2 = 0.0_rt;
        gls_s4 = 0.0_rt;
        gls_s5 = 0.0_rt;
        gls_s6 = 0.0_rt;
        gls_b0 = 0.0_rt;
        gls_b1 = 0.0_rt;
        gls_b2 = 0.0_rt;
        gls_b3 = 0.0_rt;
        gls_b4 = 0.0_rt;
        gls_b5 = 0.0_rt;
        //my_Sm1=solverChoice.my_A1*solverChoice.my_A2*((solverChoice.my_B2-3.0_rt*solverChoice.my_A2)*
        //                    (1.0_rt-6.0_rt*solverChoice.my_A1/solverChoice.my_B1)-
        //                    3.0_rt*solverChoice.my_C1*(solverChoice.my_B2+6.0_rt*solverChoice.my_A1));
        my_Sm2=9.0_rt*solverChoice.my_A1*solverChoice.my_A2;
        my_Sm3=solverChoice.my_A1*(1.0_rt-3.0_rt*solverChoice.my_C1-6.0_rt*solverChoice.my_A1/solverChoice.my_B1);
        my_Sm4=18.0_rt*solverChoice.my_A1*solverChoice.my_A1+9.0_rt*solverChoice.my_A1*solverChoice.my_A2;
        my_Sh1=solverChoice.my_A2*(1.0_rt-6.0_rt*solverChoice.my_A1/solverChoice.my_B1);
        my_Sh2=3.0_rt*solverChoice.my_A2*(6.0_rt*solverChoice.my_A1+solverChoice.my_B2);
    }

    Real Zos_min = std::max(solverChoice.Zos, 0.0001_rt);
    Real Zob_min = std::max(solverChoice.Zob, 0.0001_rt);
    Real Zos_eff = Zos_min;
    Real Gadv = 1.0_rt/3.0_rt;
    Real eps = 1.0e-10_rt;

    auto ic_bc_type = solverChoice.ic_bc_type;

    const BoxArray&            ba = cons_old[lev]->boxArray();
    const DistributionMapping& dm = cons_old[lev]->DistributionMap();
    MultiFab mf_dU(convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_dV(convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW,NGROW,0));

    MultiFab mf_CF(convert(ba,IntVect(0,0,1)),dm,1,IntVect(NGROW,NGROW,0));

    MultiFab mf_shear2(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_shear2_cached(ba,dm,1,IntVect(NGROW,NGROW,0));

    MultiFab mf_buoy2(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_tmp_buoy(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_tmp_shear(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_curvK(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_curvP(ba,dm,1,IntVect(NGROW,NGROW,0));

    MultiFab mf_FXK(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_FXP(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_FEK(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_FEP(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_FCK(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_FCP(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_BCK(ba,dm,1,IntVect(NGROW,NGROW,0));
    MultiFab mf_BCP(ba,dm,1,IntVect(NGROW,NGROW,0));

    const Box& domain = geom[0].Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    GeometryData const& geomdata = geom[0].data();
    bool is_periodic_in_x = geomdata.isPeriodic(0);
    bool is_periodic_in_y = geomdata.isPeriodic(1);


    for ( MFIter mfi(*mf_gls, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box   bx = mfi.tilebox();
        Box gbx1 = mfi.growntilebox(IntVect(NGROW-1,NGROW-1,0));

        Box bxD = bx;
        bxD.makeSlab(2,0);
        Box gbx1D = gbx1;
        gbx1D.makeSlab(2,0);

        Array4<Real> const& Hz = vec_Hz[lev]->array(mfi);
        Array4<Real> const& u = xvel_old[lev]->array(mfi);
        Array4<Real> const& v = yvel_old[lev]->array(mfi);

        auto dU = mf_dU.array(mfi);
        auto dV = mf_dV.array(mfi);
        auto CF = mf_CF.array(mfi);
        auto shear2_cached = mf_shear2_cached.array(mfi);

        ParallelFor(gbx1D, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            CF(i,j,0) = 0.0_rt;
            dU(i,j,0) = 0.0_rt;
            dV(i,j,0) = 0.0_rt;
            for (int k=1; k<=N; k++) {
                Real cff = 1.0_rt / (2.0_rt * Hz(i,j,k) + Hz(i,j,k-1)*(2.0_rt - CF(i,j,k-1)));
                CF(i,j,k) = cff * Hz(i,j,k);
                dU(i,j,k)=cff*(3.0_rt*(u(i  ,j,k)-u(i,  j,k-1)+
                                       u(i+1,j,k)-u(i+1,j,k-1))-Hz(i,j,k-1)*dU(i,j,k-1));
                dV(i,j,k)=cff*(3.0_rt*(v(i,j  ,k)-v(i,j  ,k-1)+
                                       v(i,j+1,k)-v(i,j+1,k-1))-Hz(i,j,k-1)*dV(i,j,k-1));
            }
            dU(i,j,N+1) = 0.0_rt;
            dV(i,j,N+1) = 0.0_rt;
            for (int k=N; k>=1; k--) {
                dU(i,j,k) = dU(i,j,k) - CF(i,j,k) * dU(i,j,k+1);
                dV(i,j,k) = dV(i,j,k) - CF(i,j,k) * dV(i,j,k+1);
            }
            shear2_cached(i,j,0) = 0.0_rt;
            for (int k=1; k<=N; k++) {
                shear2_cached(i,j,k) = dU(i,j,k) * dU(i,j,k) + dV(i,j,k) * dV(i,j,k);
            }
        });
    }

    // While potentially counterintuitive, this is what ROMS does for handling shear2 at all boundaries, even
    // periodic
    (*physbcs[lev])(mf_shear2_cached,*mf_mskr,0,1,mf_shear2_cached.nGrowVect(),t_new[lev],BCVars::foextrap_bc);
    mf_CF.setVal(0.0_rt);

    for ( MFIter mfi(*mf_gls, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box  bx = mfi.tilebox();
        Box xbx = surroundingNodes(bx,0);
        Box ybx = surroundingNodes(bx,1);
        Box gbx1 = grow(bx,IntVect(NGROW-1,NGROW-1,0));

        Box bx_rho = bx;
        bx_rho.convert(IntVect(0,0,0));

        Box xbx_int = grow(xbx,IntVect(0,0,-1)); // interior points in w of xbx
        Box ybx_int = grow(ybx,IntVect(0,0,-1)); // interior points in w of ybx
        Box bx_growloxy = growLo(growLo(grow(bx,IntVect(0,0,-1)),0,1),1,1);

        Box bxD = bx;
        bxD.makeSlab(2,0);
        Box gbx1D = gbx1;
        gbx1D.makeSlab(2,0);

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

        Array4<Real> const& W = mf_W.array(mfi);
        Array4<Real> const& Hz = vec_Hz[lev]->array(mfi);
        Array4<Real> const& pm = vec_pm[lev]->array(mfi);
        Array4<Real> const& pn = vec_pn[lev]->array(mfi);
        Array4<Real> const& Lscale = vec_Lscale[lev]->array(mfi);

        Array4<Real> const& Huon = vec_Huon[lev]->array(mfi);
        Array4<Real> const& Hvom = vec_Hvom[lev]->array(mfi);
        Array4<Real> const& z_w = vec_z_w[lev]->array(mfi);

        Array4<Real> const& tke = mf_tke->array(mfi);
        Array4<Real> const& gls = mf_gls->array(mfi);

        Array4<Real> const& sustr = vec_sustr[lev]->array(mfi);
        Array4<Real> const& svstr = vec_svstr[lev]->array(mfi);
        Array4<Real> const& bustr = vec_bustr[lev]->array(mfi);
        Array4<Real> const& bvstr = vec_bvstr[lev]->array(mfi);
        Array4<Real> const& msku = vec_msku[lev]->array(mfi);
        Array4<Real> const& mskv = vec_mskv[lev]->array(mfi);

        FArrayBox fab_tmp_buoy(bx_growloxy,1, amrex::The_Async_Arena()); fab_tmp_buoy.template setVal<RunOn::Device>(0.);
        FArrayBox fab_tmp_shear(bx_growloxy,1, amrex::The_Async_Arena()); fab_tmp_shear.template setVal<RunOn::Device>(0.);

        FArrayBox fab_curvK(gbx1,1, amrex::The_Async_Arena()); fab_curvK.template setVal<RunOn::Device>(0.);
        FArrayBox fab_curvP(gbx1,1, amrex::The_Async_Arena()); fab_curvP.template setVal<RunOn::Device>(0.);

        FArrayBox fab_FXK(xbx_int,1,amrex::The_Async_Arena()); fab_FXK.template setVal<RunOn::Device>(0.);
        FArrayBox fab_FXP(xbx_int,1,amrex::The_Async_Arena()); fab_FXP.template setVal<RunOn::Device>(0.);
        FArrayBox fab_FEK(ybx_int,1,amrex::The_Async_Arena()); fab_FEK.template setVal<RunOn::Device>(0.);
        FArrayBox fab_FEP(ybx_int,1,amrex::The_Async_Arena()); fab_FEP.template setVal<RunOn::Device>(0.);

        FArrayBox fab_FCK(bx_rho,1,amrex::The_Async_Arena()); fab_FCK.template setVal<RunOn::Device>(0.);
        FArrayBox fab_FCP(bx_rho,1,amrex::The_Async_Arena()); fab_FCP.template setVal<RunOn::Device>(0.);
        FArrayBox fab_BCK(bx_rho,1,amrex::The_Async_Arena()); fab_BCK.template setVal<RunOn::Device>(0.);
        FArrayBox fab_BCP(bx_rho,1,amrex::The_Async_Arena()); fab_BCP.template setVal<RunOn::Device>(0.);

        auto CF = mf_CF.array(mfi);
        auto shear2 = mf_shear2.array(mfi);
        auto shear2_cached = mf_shear2_cached.array(mfi);
        auto buoy2 = mf_buoy2.array(mfi);
        Array4<Real> const& bvf = vec_bvf[lev]->array(mfi);

        auto tmp_buoy = fab_tmp_buoy.array();
        auto tmp_shear = fab_tmp_shear.array();
        auto curvK = fab_curvK.array();
        auto curvP = fab_curvP.array();
        auto FXK = fab_FXK.array();
        auto FXP = fab_FXP.array();
        auto FEK = fab_FEK.array();
        auto FEP = fab_FEP.array();
        auto FCK = fab_FCK.array();
        auto FCP = fab_FCP.array();
        auto BCK = fab_BCK.array();
        auto BCP = fab_BCP.array();

        auto Akt = mf_Akt->array(mfi);
        auto Akv = mf_Akv->array(mfi);
        auto Akp = mf_Akp->array(mfi);
        auto Akk = mf_Akk->array(mfi);

        ParallelFor(bx_growloxy, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            tmp_buoy(i,j,k)=0.25_rt * (bvf(i,j,k) + bvf(i+1,j,k) + bvf(i,j+1,k)+bvf(i+1,j+1,k));
            tmp_shear(i,j,k)=0.25_rt * (shear2_cached(i,j,k) + shear2_cached(i+1,j,k) + shear2_cached(i,j+1,k)+shear2_cached(i+1,j+1,k));
        });

        ParallelFor(grow(bx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            buoy2(i,j,k)=0.25_rt * (tmp_buoy(i,j,k) + tmp_buoy(i-1,j,k) + tmp_buoy(i,j-1,k)+tmp_buoy(i-1,j-1,k));
            shear2(i,j,k)=0.25_rt * (tmp_shear(i,j,k) + tmp_shear(i-1,j,k) + tmp_shear(i,j-1,k)+tmp_shear(i-1,j-1,k));
        });

        //Time step advective terms
        ParallelFor(growLo(grow(xbx,IntVect(0,0,-1)),0,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real gradK, gradK_ip1, gradP, gradP_ip1;

            if (i == dlo.x-1 && !is_periodic_in_x) {
                gradK_ip1 = tke(i+1,j,k,2)-tke(i  ,j,k,2);
                gradK = gradK_ip1;
                gradP_ip1 = gls(i+1,j,k,2)-gls(i  ,j,k,2);
                gradP = gradP_ip1;
            } else if (i == dhi.x+1 && !is_periodic_in_x) {
                gradK = tke(i  ,j,k,2)-tke(i-1,j,k,2);
                gradK_ip1 = gradK;
                gradP = gls(i  ,j,k,2)-gls(i-1,j,k,2);
                gradP_ip1 = gradP;
            } else {
                gradK     = (tke(i  ,j,k,2)-tke(i-1,j,k,2)) * msku(i  ,j,0);
                gradK_ip1 = (tke(i+1,j,k,2)-tke(i  ,j,k,2)) * msku(i+1,j,0);
                gradP     = (gls(i  ,j,k,2)-gls(i-1,j,k,2)) * msku(i  ,j,0);
                gradP_ip1 = (gls(i+1,j,k,2)-gls(i  ,j,k,2)) * msku(i+1,j,0);
            }

            curvK(i,j,k) = gradK_ip1 - gradK;
            curvP(i,j,k) = gradP_ip1 - gradP;
        });
        ParallelFor(grow(xbx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff = 0.5_rt * (Huon(i,j,k) + Huon(i,j,k-1));
            Real cff1 = (cff > 0.0) ? curvK(i-1,j,k) : curvK(i,j,k);
            Real cff2 = (cff > 0.0) ? curvP(i-1,j,k) : curvP(i,j,k);

            FXK(i,j,k) = cff * 0.5_rt * (tke(i-1,j,k,2)+tke(i,j,k,2)-Gadv*cff1);
            FXP(i,j,k) = cff * 0.5_rt * (gls(i-1,j,k,2)+gls(i,j,k,2)-Gadv*cff2);
        });

        //Time step advective terms
        ParallelFor(growLo(grow(ybx,IntVect(0,0,-1)),1,1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real gradK     = (tke(i,j  ,k,2)-tke(i,j-1,k,2)) * mskv(i,j  ,0);
            Real gradK_jp1 = (tke(i,j+1,k,2)-tke(i,j  ,k,2)) * mskv(i,j+1,0);
            Real gradP     = (gls(i,j  ,k,2)-gls(i,j-1,k,2)) * mskv(i,j  ,0);
            Real gradP_jp1 = (gls(i,j+1,k,2)-gls(i,j  ,k,2)) * mskv(i,j+1,0);

            if (j == dlo.y-1 && !is_periodic_in_y) {
                gradK = gradK_jp1;
                gradP = gradP_jp1;
            }
            else if (j == dhi.y+1 && !is_periodic_in_y) {
                gradK_jp1 = gradK;
                gradP_jp1 = gradP;
            }

            curvK(i,j,k) = gradK_jp1 - gradK;
            curvP(i,j,k) = gradP_jp1 - gradP;
        });
        ParallelFor(grow(ybx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff = 0.5_rt * (Hvom(i,j,k) + Hvom(i,j,k-1));
            Real cff1 = (cff > 0.0) ? curvK(i,j-1,k) : curvK(i,j,k);
            Real cff2 = (cff > 0.0) ? curvP(i,j-1,k) : curvP(i,j,k);

            FEK(i,j,k) = cff * 0.5_rt * (tke(i,j-1,k,2)+tke(i,j,k,2)-Gadv*cff1);
            FEP(i,j,k) = cff * 0.5_rt * (gls(i,j-1,k,2)+gls(i,j,k,2)-Gadv*cff2);
        });

        Real gls_Kmin = solverChoice.gls_Kmin;
        Real gls_Pmin = solverChoice.gls_Pmin;
        ParallelFor(grow(bx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff = dt_lev * pm(i,j,0) * pn(i,j,0);
            tke(i,j,k,nnew) = tke(i,j,k,nnew) - cff * (FXK(i+1,j  ,k)-FXK(i,j,k)+
                                                       FEK(i  ,j+1,k)-FEK(i,j,k));
            tke(i,j,k,nnew) = std::max(tke(i,j,k,nnew), gls_Kmin);

            gls(i,j,k,nnew) = gls(i,j,k,nnew) - cff * (FXP(i+1,j  ,k)-FXP(i,j,k)+
                                                       FEP(i  ,j+1,k)-FEP(i,j,k));
            gls(i,j,k,nnew) = std::max(gls(i,j,k,nnew), gls_Pmin);
        });


        // Vertical advection
        ParallelFor(bxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Real cff1 = 7.0_rt / 12.0_rt;
            Real cff2 = 1.0_rt / 12.0_rt;
            for (int k=1; k<=N-1; k++) {
                Real cff = 0.5_rt * (W(i,j,k+1)+W(i,j,k));
                FCK(i,j,k) = cff * (cff1 * (tke(i,j,k  ,2)+tke(i,j,k+1,2))-
                                    cff2 * (tke(i,j,k-1,2)+tke(i,j,k+2,2)));
                FCP(i,j,k) = cff * (cff1 * (gls(i,j,k  ,2)+gls(i,j,k+1,2))-
                                    cff2 * (gls(i,j,k-1,2)+gls(i,j,k+2,2)));
            }
            cff1 = 1.0_rt/3.0_rt;
            cff2 = 5.0_rt/6.0_rt;
            Real cff3 = 1.0_rt / 6.0_rt;
            Real cff = 0.5_rt * (W(i,j,0)+W(i,j,1));
            FCK(i,j,0) = cff * (cff1 * tke(i,j,0,2)+cff2 * tke(i,j,1,2)-cff3 * tke(i,j,2,2));
            FCP(i,j,0) = cff * (cff1 * gls(i,j,0,2)+cff2 * gls(i,j,1,2)-cff3 * gls(i,j,2,2));
            cff = 0.5_rt * (W(i,j,N+1)+W(i,j,N));
            FCK(i,j,N) = cff * (cff1 * tke(i,j,N+1,2)+cff2*tke(i,j,N,2)-cff3*tke(i,j,N-1,2));
            FCP(i,j,N) = cff * (cff1 * gls(i,j,N+1,2)+cff2*gls(i,j,N,2)-cff3*gls(i,j,N-1,2));
        });
        ParallelFor(grow(bx,2,-1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff = dt_lev * pm(i,j,0) * pn(i,j,0);
            tke(i,j,k,nnew) = tke(i,j,k,nnew) - cff*(FCK(i,j,k  )-FCK(i,j,k-1));
            tke(i,j,k,nnew) = std::max(tke(i,j,k,nnew),gls_Kmin);
            gls(i,j,k,nnew) = gls(i,j,k,nnew) - cff*(FCP(i,j,k  )-FCP(i,j,k-1));
            gls(i,j,k,nnew) = std::max(gls(i,j,k,nnew),gls_Pmin);
        });

        // Compute vertical mixing, turbulent production and turbulent
        // dissipation.
        //
        Real cff = -0.5 * dt_lev;
        ParallelFor(convert(bx,IntVect(0,0,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (k==0 or k==N) {
                FCK(i,j,k) = 0.0_rt;
                FCP(i,j,k) = 0.0_rt;
            } else {
                FCK(i,j,k) = cff * (Akk(i,j,k) + Akk(i,j,k+1)) / Hz(i,j,k);
                FCP(i,j,k) = cff * (Akp(i,j,k) + Akp(i,j,k+1)) / Hz(i,j,k);
            }
        });
        // Compute production and dissipation terms.
        ParallelFor(grow(bx,2,-1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // Compute shear and buoyant production of turbulent energy (m3/s3)
            // at W-points (ignore small negative values of buoyancy).
            Real strat2 = buoy2(i,j,k);
            Real gls_c3 = (strat2 > 0.0) ? gls_c3m : gls_c3p;
            Real Kprod = shear2(i,j,k) * (Akv(i,j,k)-Akv_bak) -
                         strat2 * (Akt(i,j,k,Temp_comp)-Akt_bak);
            Real Pprod = gls_c1 * shear2(i,j,k) * (Akv(i,j,k)-Akv_bak) -
                         gls_c3 * strat2 * (Akt(i,j,k,Temp_comp)-Akt_bak);

            // If negative production terms, then add buoyancy to dissipation terms
            // (BCK and BCP) below, using "cff1" and "cff2" as the on/off switch.
            Real cff1 = (Kprod < 0.0_rt) ? 0.0_rt : 1.0_rt;
            Real cff2 = (Pprod < 0.0_rt) ? 0.0_rt : 1.0_rt;
            if (Kprod < 0.0_rt) {
                Kprod = Kprod + strat2*(Akt(i,j,k,Temp_comp)-Akt_bak);
            }
            if (Pprod < 0.0_rt) {
                Pprod = Pprod + gls_c3*strat2*(Akt(i,j,k,Temp_comp)-Akt_bak);
            }
            // Time-step shear and buoyancy production terms.
            Real cff_Hz = 0.5_rt * (Hz(i,j,k) + Hz(i,j,k-1));
            tke(i,j,k,nnew) = tke(i,j,k,nnew)+dt_lev * cff_Hz * Kprod;
            gls(i,j,k,nnew) = gls(i,j,k,nnew)+dt_lev
                                *cff_Hz*Pprod*gls(i,j,k,nstp) / std::max(tke(i,j,k,nstp),gls_Kmin);

            // Compute dissipation of turbulent energy (m3/s3).
            Real wall_fac = 1.0_rt;
            if (Lmy25) {
                // Parabolic wall function,  L = ds db / (ds + db).
                wall_fac=1.0_rt+gls_E2/(vonKar*vonKar)*
                        std::pow(std::pow(gls(i,j,k,nstp),( gls_exp1))*cmu_fac1*
                         std::pow(tke(i,j,k,nstp),-tke_exp1)*
                         (1.0_rt/ (z_w(i,j,k)-z_w(i,j,0))),2)+
                        0.25_rt/(vonKar*vonKar)*
                        std::pow(std::pow(gls(i,j,k,nstp), gls_exp1)*cmu_fac1*
                         std::pow(tke(i,j,k,nstp),-tke_exp1)*
                         (1.0_rt/ (z_w(i,j,N+1)-z_w(i,j,k))),2);
            }
            BCK(i,j,k)=cff_Hz*(1.0_rt+dt_lev*
                          std::pow(gls(i,j,k,nstp),(-gls_exp1))*cmu_fac2*
                          std::pow(tke(i,j,k,nstp), tke_exp2)+
                          dt_lev*(1.0_rt-cff1)*strat2*
                          (Akt(i,j,k,Temp_comp)-Akt_bak)/
                          tke(i,j,k,nstp))-
                          FCK(i,j,k)-FCK(i,j,k-1);
            BCP(i,j,k)=cff_Hz*(1.0_rt+dt_lev*gls_c2*wall_fac*
                          std::pow(gls(i,j,k,nstp),-gls_exp1)*cmu_fac2*
                          std::pow(tke(i,j,k,nstp), tke_exp2)+
                          dt_lev*(1.0_rt-cff2)*gls_c3*strat2*
                          (Akt(i,j,k,Temp_comp)-Akt_bak)/
                          tke(i,j,k,nstp))-
                          FCP(i,j,k)-FCP(i,j,k-1);
        });
        // Compute production and dissipation terms.
        ParallelFor(bxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            //----------------------------------------------------------------------
            // Time-step dissipation and vertical diffusion terms implicitly.
            //----------------------------------------------------------------------
            //
            // Set Dirichlet surface and bottom boundary conditions. Compute
            // surface roughness from wind stress (Charnok) and set Craig and
            // Banner wave breaking surface flux, if appropriate.


            tke(i,j,N+1,nnew)=std::max(cmu_fac3*0.5_rt*
                                     std::sqrt((sustr(i,j,0)+sustr(i+1,j,0))*(sustr(i,j,0)+sustr(i+1,j,0))+
                                          (svstr(i,j,0)+svstr(i,j+1,0))*(svstr(i,j,0)+svstr(i,j+1,0))),
                                     gls_Kmin);
            tke(i,j,0,nnew)=std::max(cmu_fac3*0.5_rt*
                                 std::sqrt((bustr(i,j,0)+bustr(i+1,j,0))*(bustr(i,j,0)+bustr(i+1,j,0))+
                                      (bvstr(i,j,0)+bvstr(i,j+1,0))*(bvstr(i,j,0)+bvstr(i,j+1,0))),
                                        gls_Kmin);

            gls(i,j,N+1,nnew)=std::max(std::pow(gls_cmu0,gls_p)*
                                    std::pow(tke(i,j,N+1,nnew),gls_m)*
                                    std::pow(L_sft*Zos_eff,gls_n), gls_Pmin);
            Real cff_gls = gls_fac4*std::pow(vonKar*Zob_min,gls_n);
            gls(i,j,0,nnew)=std::max(cff_gls*std::pow(tke(i,j,0,nnew),(gls_m)), gls_Pmin);

            // Solve tri-diagonal system for turbulent kinetic energy.
            // Might be N instead of N-1?
            Real tke_fluxt = 0.0_rt;
            Real tke_fluxb = 0.0_rt;
            Real cff_BCK = 1.0_rt/BCK(i,j,N);
            CF(i,j,N)=cff_BCK*FCK(i,j,N-1);
            tke(i,j,N,nnew)=cff_BCK*(tke(i,j,N,nnew)+tke_fluxt);
            for (int k=N-1;k>=1;k--) {
                cff_BCK = 1.0_rt / (BCK(i,j,k)-CF(i,j,k+1)*FCK(i,j,k));
                CF(i,j,k) = cff_BCK * FCK(i,j,k-1);
                tke(i,j,k,nnew) = cff_BCK * (tke(i,j,k,nnew) - FCK(i,j,k) * tke(i,j,k+1,nnew));
            }
            tke(i,j,1,nnew) = tke(i,j,1,nnew) - cff_BCK * tke_fluxb;
            tke(i,j,1,nnew) = std::max(tke(i,j,1,nnew),gls_Kmin);
            for (int k=2;k<=N;k++) {
                tke(i,j,k,nnew) = tke(i,j,k,nnew) - CF(i,j,k) * tke(i,j,k-1,nnew);
                tke(i,j,k,nnew) = std::max(tke(i,j,k,nnew), gls_Kmin);
            }

            // Solve tri-diagonal system for generic statistical field.
            Real cff_tke = 0.5_rt * (tke(i,j,N+1,nnew) + tke(i,j,N,nnew));
            Real gls_fluxt = dt_lev*gls_fac3*std::pow(cff_tke,gls_m)*
                             std::pow(L_sft,(gls_n))*
                             std::pow(Zos_eff+0.5_rt*Hz(i,j,N),gls_n-1.0_rt)*
                             0.5_rt*(Akp(i,j,N+1)+Akp(i,j,N));
            cff_tke=0.5_rt*(tke(i,j,0,nnew)+tke(i,j,1,nnew));
            Real gls_fluxb = dt_lev*gls_fac2*std::pow(cff_tke,gls_m)*
                              std::pow(0.5_rt*Hz(i,j,0)+Zob_min,gls_n-1.0_rt)*
                              0.5_rt*(Akp(i,j,0)+Akp(i,j,1));
            Real cff_BCP = 1.0_rt / BCP(i,j,N);
            CF(i,j,N) = cff_BCP * FCP(i,j,N-1);
            gls(i,j,N,nnew)=cff_BCP*(gls(i,j,N,nnew)-gls_fluxt);
            for (int k=N-1;k>=1;k--) {
                cff_BCP = 1.0_rt / (BCP(i,j,k)-CF(i,j,k+1)*FCP(i,j,k));
                CF(i,j,k) = cff_BCP * FCP(i,j,k-1);
                gls(i,j,k,nnew) = cff_BCP * (gls(i,j,k,nnew) - FCP(i,j,k)*gls(i,j,k+1,nnew));
            }
            gls(i,j,1,nnew) = gls(i,j,1,nnew)-cff_BCP*gls_fluxb;
            for (int k=2; k<=N; k++) {
                gls(i,j,k,nnew) = gls(i,j,k,nnew) - CF(i,j,k) * gls(i,j,k-1,nnew);
            }
        });

        // Compute vertical mixing coefficients (m2/s).
        ParallelFor(grow(bx,2,-1), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            tke(i,j,k,nnew) = std::max(tke(i,j,k,nnew),gls_Kmin);
            gls(i,j,k,nnew) = std::max(gls(i,j,k,nnew),gls_Pmin);
            if (gls_n >= 0.0_rt) {
                gls(i,j,k,nnew)=std::min(gls(i,j,k,nnew),gls_fac5*
                                    std::pow(tke(i,j,k,nnew),tke_exp4)*
                                    std::pow(std::sqrt(std::max(0.0_rt,
                                          buoy2(i,j,k)))+eps,-gls_n));
            } else {
                gls(i,j,k,nnew)=std::max(gls(i,j,k,nnew), gls_fac5*
                                 std::pow(tke(i,j,k,nnew),(tke_exp4))*
                                 std::pow((std::sqrt(std::max(0.0_rt,
                                       buoy2(i,j,k)))+eps),(-gls_n)));
            }
            Real Ls_lmt;
            Real Ls_unlmt=std::max(eps,
                                   std::pow(gls(i,j,k,nnew),( gls_exp1))*cmu_fac1*
                                   std::pow(tke(i,j,k,nnew),(-tke_exp1)));
            if (buoy2(i,j,k) > 0.0_rt) {
                Ls_lmt=std::min(Ls_unlmt,
                                std::sqrt(0.56_rt*tke(i,j,k,nnew)/
                                (std::max(0.0_rt,buoy2(i,j,k))+eps)));
            } else {
                Ls_lmt = Ls_unlmt;
            }
            //
            //  Recompute gls based on limited length scale
            //
            gls(i,j,k,nnew)=std::max(std::pow(gls_cmu0,gls_p)*
                                           std::pow(tke(i,j,k,nnew),gls_m)*
                                           std::pow(Ls_lmt,gls_n), gls_Pmin);

            //   Compute nondimensional stability functions for tracers (Sh) and
            //   momentum (Sm).
            Real Sh, Sm;
            Real Gh=std::min(gls_Gh0,-buoy2(i,j,k)*Ls_lmt*Ls_lmt/
                            (2.0_rt*tke(i,j,k,nnew)));
            Gh=std::min(Gh,Gh-(Gh-gls_Ghcri)*(Gh-gls_Ghcri)/
                       (Gh+gls_Gh0-2.0_rt*gls_Ghcri));
            Gh=std::max(Gh,gls_Ghmin);

            if (gls_stability_type == GLS_StabilityType::Canuto_A ||
                gls_stability_type == GLS_StabilityType::Canuto_B) {
                //
                //   Canuto stability: Compute shear number.
                //
                Real Gm=(gls_b0/gls_fac6-gls_b1*Gh+gls_b3*gls_fac6*(Gh*Gh))/
                             (gls_b2-gls_b4*gls_fac6*Gh);
                Gm=std::min(Gm,shear2(i,j,k)*Ls_lmt*Ls_lmt/
                            (2.0_rt*tke(i,j,k,nnew)));
                /////Gm=std::min(Gm,(gls_s1*gls_fac6*Gh-gls_s0)/(gls_s2*gls_fac6));
                //
                //  Compute stability functions
                //
                Real stab_cff=gls_b0-gls_b1*gls_fac6*Gh+gls_b2*gls_fac6*Gm+
                    gls_b3*gls_fac6*gls_fac6*Gh*Gh-gls_b4*gls_fac6*gls_fac6*Gh*Gm+
                    gls_b5*gls_fac6*gls_fac6*Gm*Gm;
                Sm=(gls_s0-gls_s1*gls_fac6*Gh+gls_s2*gls_fac6*Gm)/stab_cff;
                Sh=(gls_s4-gls_s5*gls_fac6*Gh+gls_s6*gls_fac6*Gm)/stab_cff;
                Sm=std::max(Sm,0.0_rt);
                Sh=std::max(Sh,0.0_rt);

                //
                //  Relate Canuto stability to ROMS notation
                //
                Real gls_cmu0_cube = gls_cmu0 * gls_cmu0 * gls_cmu0;
                Sm=Sm*sqrt2/(gls_cmu0_cube);
                Sh=Sh*sqrt2/gls_cmu0_cube;
            } else if (gls_stability_type == GLS_StabilityType::Galperin) {
                Real cff_galperin = 1.0_rt - my_Sh2*Gh;
                Sh = my_Sh1 / cff_galperin;
                Sm = (my_Sm3+Sh*Gh*my_Sm4)/(1.0_rt-my_Sm2*Gh);
            }

            //  Compute vertical mixing (m2/s) coefficients of momentum and
            //  tracers.  Average ql over the two timesteps rather than using
            //  the new Lscale and just averaging tke.

            Real ql=sqrt2*0.5_rt*(Ls_lmt*std::sqrt(tke(i,j,k,nnew))+
                                  Lscale(i,j,k)*std::sqrt(tke(i,j,k,nstp)));
            Akv(i,j,k)=Akv_bak+Sm*ql;
            for (int n=0; n<NCONS; n++) {
                Akt(i,j,k,n)=Akt_bak+Sh*ql;
            }

            //  Compute vertical mixing (m2/s) coefficients of turbulent kinetic
            //  energy and generic statistical field.

            Akk(i,j,k)=Akk_bak+Sm*ql/gls_sigk;
            Akp(i,j,k)=Akp_bak+Sm*ql*ogls_sigp;

            //  Save limited length scale.
            Lscale(i,j,k)=Ls_lmt;
        });
        ParallelFor(bxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Akv(i,j,N+1)=Akv_bak+L_sft*Zos_eff*gls_cmu0*
                          std::sqrt(tke(i,j,N+1,nnew));
            Akv(i,j,0)=Akv_bak+vonKar*Zob_min*gls_cmu0*
                      std::sqrt(tke(i,j,0,nnew));

            Akk(i,j,N+1)=Akk_bak+Akv(i,j,N+1)/gls_sigk;
            Akk(i,j,0)=Akk_bak+Akv(i,j,0)/gls_sigk;
            Akp(i,j,N+1)=Akp_bak+Akv(i,j,N+1)*ogls_sigp;
            Akp(i,j,0)=Akp_bak+Akv(i,j,0)/gls_sigp_cb;

            for (int n=0; n<NCONS; n++) {
                Akt(i,j,N+1,n)  = Akt_bak;
                Akt(i,j,0,n) = Akt_bak;
            }
        });
    }

    for (int icomp=0; icomp<3; icomp++) {
        FillPatch(lev, t_old[lev], *mf_tke, GetVecOfPtrs(vec_tke), BCVars::zvel_bc, BdyVars::null, icomp, false, false);
        FillPatch(lev, t_old[lev], *mf_gls, GetVecOfPtrs(vec_gls), BCVars::zvel_bc, BdyVars::null, icomp, false, false);
    }
    for (int icomp=0; icomp<NCONS; icomp++) {
        FillPatch(lev, t_old[lev], *mf_Akt, GetVecOfPtrs(vec_Akt), BCVars::zvel_bc, BdyVars::null, icomp, false, false);
    }
    FillPatchNoBC(lev, t_old[lev], *mf_Akv, GetVecOfPtrs(vec_Akv), BdyVars::null);
    FillPatchNoBC(lev, t_old[lev], *mf_Akp, GetVecOfPtrs(vec_Akp), BdyVars::null);
    FillPatchNoBC(lev, t_old[lev], *mf_Akk, GetVecOfPtrs(vec_Akk), BdyVars::null);
}
