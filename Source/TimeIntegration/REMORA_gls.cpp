#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] lev            level to operate on
 * @param[inout] mf_gls         turbulent generic length scale
 * @param[inout] mf_tke         turbulent kinetic energy
 * @param[in   ] mf_W           vertical velocity
 * @param[in   ] mf_msku        land-sea mask on u points
 * @param[in   ] mf_mskv        land-sea mask on v points
 * @param[in   ] nstp           index of last time step in gls and tke MultiFabs
 * @param[in   ] nnew           index of time step to update in gls and tke MultiFabs
 * @param[in   ] iic            which time step we're on
 * @param[in   ] ntfirst        what is the first time step?
 * @param[in   ] N              number of vertical levels
 * @param[in   ] dt_lev         time step at this level
 */
void
REMORA::gls_prestep (int lev, MultiFab* mf_gls, MultiFab* mf_tke,
                     MultiFab& mf_W, MultiFab* mf_msku, MultiFab* mf_mskv,
                     const int nstp, const int nnew,
                     const int iic, const int ntfirst, const int N, const Real dt_lev)
{
    BL_PROFILE("REMORA::gls_prestep()");
    // temps: grad, gradL, XF, FX, FXL, EF, FE, FEL
    for ( MFIter mfi(*mf_gls, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<Real> const& gls = mf_gls->array(mfi);
        Array4<Real> const& tke = mf_tke->array(mfi);
        Array4<Real const> const& W = mf_W.const_array(mfi);

        Array4<Real const> const& Huon = vec_Huon[lev]->const_array(mfi);
        Array4<Real const> const& Hvom = vec_Hvom[lev]->const_array(mfi);
        Array4<Real const> const& Hz = vec_Hz[lev]->const_array(mfi);
        Array4<Real const> const& pm = vec_pm[lev]->const_array(mfi);
        Array4<Real const> const& pn = vec_pn[lev]->const_array(mfi);
        Array4<Real const> const& msku = mf_msku->const_array(mfi);
        Array4<Real const> const& mskv = mf_mskv->const_array(mfi);

        Box bx = mfi.tilebox();
        Box xbx = surroundingNodes(bx,0);
        Box ybx = surroundingNodes(bx,1);

        Box xbx_hi = growHi(xbx,0,1);

        Box ybx_hi = growHi(ybx,0,1);

        const Box& domain = geom[lev].Domain();
        const auto dlo = amrex::lbound(domain);
        const auto dhi = amrex::ubound(domain);

        GeometryData const& geomdata = geom[0].data();
        bool is_periodic_in_x = geomdata.isPeriodic(0);
        bool is_periodic_in_y = geomdata.isPeriodic(1);

        int ncomp = 1;
        Vector<BCRec> bcrs_x(ncomp);
        Vector<BCRec> bcrs_y(ncomp);
        amrex::setBC(xbx,domain,xvel_bc(),0,1,domain_bcs_type,bcrs_x);
        amrex::setBC(ybx,domain,yvel_bc(),0,1,domain_bcs_type,bcrs_y);

        FArrayBox fab_XF(xbx_hi, 1, amrex::The_Async_Arena()); fab_XF.template setVal<RunOn::Device>(zero);
        FArrayBox fab_FX(xbx_hi, 1, amrex::The_Async_Arena()); fab_FX.template setVal<RunOn::Device>(zero);
        FArrayBox fab_FXL(xbx_hi, 1, amrex::The_Async_Arena()); fab_FXL.template setVal<RunOn::Device>(zero);
        FArrayBox fab_EF(ybx_hi, 1, amrex::The_Async_Arena()); fab_EF.template setVal<RunOn::Device>(zero);
        FArrayBox fab_FE(ybx_hi, 1, amrex::The_Async_Arena()); fab_FE.template setVal<RunOn::Device>(zero);
        FArrayBox fab_FEL(ybx_hi, 1, amrex::The_Async_Arena()); fab_FEL.template setVal<RunOn::Device>(zero);
        FArrayBox fab_Hz_half(bx, 1, amrex::The_Async_Arena()); fab_Hz_half.template setVal<RunOn::Device>(zero);
        FArrayBox fab_CF(convert(bx,IntVect(0,0,0)), 1, amrex::The_Async_Arena()); fab_CF.template setVal<RunOn::Device>(zero);
        FArrayBox fab_FC(convert(bx,IntVect(0,0,0)), 1, amrex::The_Async_Arena()); fab_FC.template setVal<RunOn::Device>(zero);
        FArrayBox fab_FCL(convert(bx,IntVect(0,0,0)), 1, amrex::The_Async_Arena()); fab_FCL.template setVal<RunOn::Device>(zero);

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
            Real cff = one/Real(6.0);
            XF(i,j,k) = Real(0.5) * (Huon(i,j,k) + Huon(i,j,k-1));
            FX(i,j,k) = XF(i,j,k) * Real(0.5) * (tke(i-1,j,k,nstp) + tke(i,j,k,nstp) -
                cff * (grad_ip1 - grad_im1));
            FXL(i,j,k) = XF(i,j,k) * Real(0.5) * (gls(i-1,j,k,nstp) + gls(i,j,k,nstp) -
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
            Real cff = one/Real(6.0);
            EF(i,j,k) = Real(0.5) * (Hvom(i,j,k) + Hvom(i,j,k-1));
            FE(i,j,k) = EF(i,j,k) * Real(0.5) * (tke(i,j-1,k,nstp) + tke(i,j,k,nstp) -
                cff * (grad_jp1 - grad_jm1));
            FEL(i,j,k) = EF(i,j,k) * Real(0.5) * (gls(i,j-1,k,nstp) + gls(i,j,k,nstp) -
                cff * (gradL_jp1 - gradL_jm1));
        });

        Real gamma = one / Real(6.0);
        Real cff1, cff2, cff3;
        int indx;
        // Time step horizontal advection
        if (iic == ntfirst) {
            cff1 = one;
            cff2 = zero;
            cff3 = Real(0.5) * dt_lev;
            indx = nstp;
        } else {
            cff1 = Real(0.5) + gamma;
            cff2 = Real(0.5) - gamma;
            cff3 = (one - gamma) * dt_lev;
            indx = 1 - nstp;
        }

        // update tke, gls from [xlo to xhi  ] by [ylo to yhi  ]
        // need XF/FX/FXL from  [xlo to xhi+1] by [ylo to yhi  ]
        // need EF/FE/FEL from  [xlo to xhi  ] by [ylo to yhi+1]
        ParallelFor(grow(bx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff = Real(0.5) * (Hz(i,j,k) + Hz(i,j,k-1));
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
            CF(i,j,k) = Real(0.5) * (W(i,j,k+1) + W(i,j,k));
            if (k == 0) {
                Real cff1_vadv = one / Real(3.0);
                Real cff2_vadv = Real(5.0) / Real(6.0);
                Real cff3_vadv = one / Real(6.0);
                FC(i,j,k)  = CF(i,j,k) * (cff1_vadv * tke(i,j,0,nstp) +
                                          cff2_vadv * tke(i,j,1,nstp) -
                                          cff3_vadv * tke(i,j,2,nstp));
                FCL(i,j,k) = CF(i,j,k) * (cff1_vadv * gls(i,j,0,nstp) +
                                          cff2_vadv * gls(i,j,1,nstp) -
                                          cff3_vadv * gls(i,j,2,nstp));
            } else if (k == N) {
                Real cff1_vadv = one / Real(3.0);
                Real cff2_vadv = Real(5.0) / Real(6.0);
                Real cff3_vadv = one / Real(6.0);
                FC(i,j,k)  = CF(i,j,k) * (cff1_vadv * tke(i,j,k+1,  nstp) +
                                          cff2_vadv * tke(i,j,k  ,nstp)-
                                          cff3_vadv * tke(i,j,k-1,nstp));
                FCL(i,j,k) = CF(i,j,k) * (cff1_vadv * gls(i,j,k+1,nstp) +
                                          cff2_vadv * gls(i,j,k  ,nstp)-
                                          cff3_vadv * gls(i,j,k-1,nstp));
            } else {
                Real cff1_vadv = Real(7.0) / Real(12.0);
                Real cff2_vadv = one / Real(12.0);
                FC(i,j,k)  = CF(i,j,k) * (cff1_vadv * (tke(i,j,k  ,nstp) + tke(i,j,k+1,nstp)) -
                                          cff2_vadv * (tke(i,j,k-1,nstp) + tke(i,j,k+2,nstp)));
                FCL(i,j,k) = CF(i,j,k) * (cff1_vadv * (gls(i,j,k  ,nstp) + gls(i,j,k+1,nstp)) -
                                          cff2_vadv * (gls(i,j,k-1,nstp) + gls(i,j,k+2,nstp)));
            }
        });

        // Time-step vertical advection
        if (iic == ntfirst) {
            cff3 = Real(0.5) * dt_lev;
        } else {
            cff3 = (one - gamma) * dt_lev;
        }
        // DO k=1,N-1
        ParallelFor(grow(bx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real cff4 = cff3 * pm(i,j,0) * pn(i,j,0);
            Hz_half(i,j,k) = Hz_half(i,j,k) - cff4 * (CF(i,j,k)-CF(i,j,k-1));
            Real cff1_loc = one / Hz_half(i,j,k);
            tke(i,j,k,2) = cff1_loc * (tke(i,j,k,2) - cff4 * (FC (i,j,k) - FC (i,j,k-1)));
            gls(i,j,k,2) = cff1_loc * (gls(i,j,k,2) - cff4 * (FCL(i,j,k) - FCL(i,j,k-1)));
        });
    }

    for (int icomp=0; icomp<3; icomp++) {
        FillPatch(lev, t_old[lev], *vec_tke[lev], GetVecOfPtrs(vec_tke), zvel_bc(), BdyVars::null, icomp, false, false);
        FillPatch(lev, t_old[lev], *vec_gls[lev], GetVecOfPtrs(vec_gls), zvel_bc(), BdyVars::null, icomp, false, false);
    }
}

/**
 * @param[in   ] lev            level to operate on
 * @param[inout] mf_gls         turbulent generic length scale
 * @param[inout] mf_tke         turbulent kinetic energy
 * @param[in   ] mf_W           vertical velocity
 * @param[inout] mf_Akv         vertical viscosity coefficient
 * @param[inout] mf_Akt         vertical diffusivity coefficients
 * @param[inout] mf_Akk         turbulent kinetic energy vertical diffusion coefficient
 * @param[inout] mf_Akp         turbulent length scale vertical diffusion coefficient
 * @param[in   ] mf_mskr        land-sea mask on rho points
 * @param[in   ] mf_msku        land-sea mask on u points
 * @param[in   ] mf_mskv        land-sea mask on v points
 * @param[in   ] nstp           index of last time step in gls and tke MultiFabs
 * @param[in   ] nnew           index of time step to update in gls and tke MultiFabs
 * @param[in   ] N              number of vertical levels
 * @param[in   ] dt_lev         time step at this level
 */
void
REMORA::gls_corrector (int lev, MultiFab* mf_gls, MultiFab* mf_tke,
                       MultiFab& mf_W, MultiFab* mf_Akv, MultiFab* mf_Akt,
                       MultiFab* mf_Akk, MultiFab* mf_Akp,
                       MultiFab* mf_mskr,
                       MultiFab* mf_msku, MultiFab* mf_mskv,
                       const int nstp, const int nnew,
                       const int N, const Real dt_lev)
{
    BL_PROFILE("REMORA::gls_corrector()");
//-----------------------------------------------------------------------
//  Compute several constants.
//-----------------------------------------------------------------------
    bool Lmy25 = ((solverChoice.gls_p == zero) &&
                  (solverChoice.gls_n == one) &&
                  (solverChoice.gls_m == one)) ? true : false;

    Real L_sft = vonKar;
    Real gls_sigp_cb = solverChoice.gls_sigp;
    Real ogls_sigp = one/gls_sigp_cb;

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
    Real Akp_bak = solverChoice.Akp_bak;
    Real Akk_bak = solverChoice.Akk_bak;

    // Akt_bak has one entry per active tracer, so it has to reach the device as an array.
    // The stratification terms below are specifically the temperature ones, and use
    // Akt_bak[Temp_comp] to match the Akt component they are paired with.
    Gpu::DeviceVector<Real> Akt_bak_d(NAT);
    Gpu::copy(Gpu::hostToDevice, solverChoice.Akt_bak.begin(),
              solverChoice.Akt_bak.begin() + NAT, Akt_bak_d.begin());
    Real const* Akt_bak = Akt_bak_d.data();

    Real gls_c1 = solverChoice.gls_c1;
    Real gls_c2 = solverChoice.gls_c2;
    Real gls_E2 = solverChoice.gls_E2;
    Real gls_sigk = solverChoice.gls_sigk;
    auto gls_stability_type = solverChoice.gls_stability_type;

    Real sqrt2 = std::sqrt(two);
    Real cmu_fac1 = std::pow(solverChoice.gls_cmu0,(-solverChoice.gls_p/solverChoice.gls_n));
    Real cmu_fac2 = std::pow(solverChoice.gls_cmu0,(Real(3.0)+solverChoice.gls_p/solverChoice.gls_n));
    Real cmu_fac3 = one/std::pow(solverChoice.gls_cmu0,two);

    Real gls_fac2 = std::pow(solverChoice.gls_cmu0,solverChoice.gls_p)*solverChoice.gls_n*std::pow(vonKar,solverChoice.gls_n);
    Real gls_fac3 = std::pow(solverChoice.gls_cmu0,solverChoice.gls_p)*solverChoice.gls_n;
    Real gls_fac4 = std::pow(solverChoice.gls_cmu0,solverChoice.gls_p);
    Real gls_fac5 = std::pow(Real(0.56),Real(0.5)*solverChoice.gls_n)*std::pow(solverChoice.gls_cmu0,solverChoice.gls_p);
    Real gls_fac6 = Real(8.0)/std::pow(solverChoice.gls_cmu0,Real(6.0));

    Real gls_exp1 = one/solverChoice.gls_n;
    Real tke_exp1 = solverChoice.gls_m/solverChoice.gls_n;
    Real tke_exp2 = Real(0.5)+solverChoice.gls_m/solverChoice.gls_n;
    Real tke_exp4 = solverChoice.gls_m+Real(0.5)*solverChoice.gls_n;

    Real cmu0_exp_p = std::pow(gls_cmu0, gls_p);
    Real gls_cmu0_cube = gls_cmu0 * gls_cmu0 * gls_cmu0;

    Real gls_s0, gls_s1, gls_s2, gls_s4, gls_s5, gls_s6;
    Real gls_b0, gls_b1, gls_b2, gls_b3, gls_b4, gls_b5;
    Real my_Sm2, my_Sh1, my_Sh2, my_Sm3, my_Sm4;

    // Compute parameters for Canuto et al. (2001) stability functions.
    // (Canuto, V.M., Cheng, H.Y., and Dubovikov, M.S., 2001: Ocean
    // turbulence. Part I: One-point closure model - momentum and
    // heat vertical diffusivities, JPO, 1413-1426).

    if (solverChoice.gls_stability_type == GLS_StabilityType::Canuto_A ||
        solverChoice.gls_stability_type == GLS_StabilityType::Canuto_B) {

        gls_s0=Real(3.0)/Real(2.0)*solverChoice.gls_L1*solverChoice.gls_L5*solverChoice.gls_L5;
        gls_s1=-solverChoice.gls_L4*(solverChoice.gls_L6+solverChoice.gls_L7)
                        +two*solverChoice.gls_L4*solverChoice.gls_L5*
                        (solverChoice.gls_L1-one/Real(3.0)*solverChoice.gls_L2-solverChoice.gls_L3)
                        +Real(3.0)/Real(2.0)*
                        solverChoice.gls_L1*solverChoice.gls_L5*solverChoice.gls_L8;
        gls_s2=Real(-3.0)/Real(8.0)*solverChoice.gls_L1
            *(solverChoice.gls_L6*solverChoice.gls_L6-solverChoice.gls_L7*solverChoice.gls_L7);
        gls_s4=two*solverChoice.gls_L5;
        gls_s5=two*solverChoice.gls_L4;
        gls_s6=two/Real(3.0)*solverChoice.gls_L5
            *(Real(3.0)*solverChoice.gls_L3*solverChoice.gls_L3-solverChoice.gls_L2*solverChoice.gls_L2)-
                    Real(0.5)*solverChoice.gls_L5*solverChoice.gls_L1*(Real(3.0)*solverChoice.gls_L3-solverChoice.gls_L2)+
                    Real(3.0)/Real(4.0)*solverChoice.gls_L1*(solverChoice.gls_L6-solverChoice.gls_L7);
        gls_b0=Real(3.0)*solverChoice.gls_L5*solverChoice.gls_L5;
        gls_b1=solverChoice.gls_L5*(Real(7.0)*solverChoice.gls_L4+Real(3.0)*solverChoice.gls_L8);
        gls_b2=solverChoice.gls_L5*solverChoice.gls_L5*(Real(3.0)*solverChoice.gls_L3*solverChoice.gls_L3-solverChoice.gls_L2*solverChoice.gls_L2)-
                    Real(3.0)/Real(4.0)*(solverChoice.gls_L6*solverChoice.gls_L6-solverChoice.gls_L7*solverChoice.gls_L7);
        gls_b3=solverChoice.gls_L4*(Real(4.0)*solverChoice.gls_L4+Real(3.0)*solverChoice.gls_L8);
        gls_b5=one/Real(4.0)*(solverChoice.gls_L2*solverChoice.gls_L2-Real(3.0)*solverChoice.gls_L3*solverChoice.gls_L3)*
                    (solverChoice.gls_L6*solverChoice.gls_L6-solverChoice.gls_L7*solverChoice.gls_L7);
        gls_b4=solverChoice.gls_L4*(solverChoice.gls_L2*solverChoice.gls_L6-Real(3.0)*solverChoice.gls_L3*solverChoice.gls_L7-
                    solverChoice.gls_L5*(solverChoice.gls_L2*solverChoice.gls_L2-solverChoice.gls_L3*solverChoice.gls_L3))+solverChoice.gls_L5*solverChoice.gls_L8*
                    (Real(3.0)*solverChoice.gls_L3*solverChoice.gls_L3-solverChoice.gls_L2*solverChoice.gls_L2);
        my_Sm2 = zero;
        my_Sm3 = zero;
        my_Sm4 = zero;
        my_Sh1 = zero;
        my_Sh2 = zero;
    } else {
        gls_s0 = zero;
        gls_s1 = zero;
        gls_s2 = zero;
        gls_s4 = zero;
        gls_s5 = zero;
        gls_s6 = zero;
        gls_b0 = zero;
        gls_b1 = zero;
        gls_b2 = zero;
        gls_b3 = zero;
        gls_b4 = zero;
        gls_b5 = zero;
        my_Sm2=Real(9.0)*solverChoice.my_A1*solverChoice.my_A2;
        my_Sm3=solverChoice.my_A1*(one-Real(3.0)*solverChoice.my_C1-Real(6.0)*solverChoice.my_A1/solverChoice.my_B1);
        my_Sm4=Real(18.0)*solverChoice.my_A1*solverChoice.my_A1+Real(9.0)*solverChoice.my_A1*solverChoice.my_A2;
        my_Sh1=solverChoice.my_A2*(one-Real(6.0)*solverChoice.my_A1/solverChoice.my_B1);
        my_Sh2=Real(3.0)*solverChoice.my_A2*(Real(6.0)*solverChoice.my_A1+solverChoice.my_B2);
    }

    Real Zos_min = std::max(solverChoice.Zos, Real(0.0001));
    Real Zos_eff = Zos_min;
    Real Gadv = one/Real(3.0);
    Real eps = Real(1.0e-10);

    const BoxArray&            ba = cons_old[lev]->boxArray();
    const DistributionMapping& dm = cons_old[lev]->DistributionMap();

    int ncomp_w = 0;
    int dU_comp = ncomp_w++;
    int dV_comp = ncomp_w++;
    int CF_comp = ncomp_w++;

    int ncomp = 0;
    int shear2_comp = ncomp++;
    int shear2_cache_comp = ncomp++;
    int buoy2_comp = ncomp++;

    MultiFab mf_w(convert(ba, IntVect(0,0,1)),dm,ncomp_w,IntVect(NGROW,NGROW,0));
    MultiFab mf(ba,dm,ncomp,IntVect(NGROW,NGROW,0));

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

        auto dU = mf_w.array(mfi,dU_comp);
        auto dV = mf_w.array(mfi,dV_comp);
        auto CF = mf_w.array(mfi,CF_comp);
        auto shear2_cached = mf.array(mfi,shear2_cache_comp);

        ParallelFor(gbx1D, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            CF(i,j,0) = zero;
            dU(i,j,0) = zero;
            dV(i,j,0) = zero;
            for (int k=1; k<=N; k++) {
                Real cff = one / (two * Hz(i,j,k) + Hz(i,j,k-1)*(two - CF(i,j,k-1)));
                CF(i,j,k) = cff * Hz(i,j,k);
                dU(i,j,k)=cff*(Real(3.0)*(u(i  ,j,k)-u(i,  j,k-1)+
                                       u(i+1,j,k)-u(i+1,j,k-1))-Hz(i,j,k-1)*dU(i,j,k-1));
                dV(i,j,k)=cff*(Real(3.0)*(v(i,j  ,k)-v(i,j  ,k-1)+
                                       v(i,j+1,k)-v(i,j+1,k-1))-Hz(i,j,k-1)*dV(i,j,k-1));
            }
            dU(i,j,N+1) = zero;
            dV(i,j,N+1) = zero;
            for (int k=N; k>=1; k--) {
                dU(i,j,k) = dU(i,j,k) - CF(i,j,k) * dU(i,j,k+1);
                dV(i,j,k) = dV(i,j,k) - CF(i,j,k) * dV(i,j,k+1);
            }
            shear2_cached(i,j,0) = zero;
            for (int k=1; k<=N; k++) {
                shear2_cached(i,j,k) = dU(i,j,k) * dU(i,j,k) + dV(i,j,k) * dV(i,j,k);
            }
        });
    }

    // While potentially counterintuitive, this is what ROMS does for handling shear2 at all boundaries, even
    // periodic
    (*physbcs[lev])(mf,*mf_mskr,shear2_cache_comp,1,mf.nGrowVect(),t_new[lev],foextrap_bc());
    mf.setVal(zero,CF_comp,1);

    int ncomp_fab = 0;
    int tmp_buoy_comp  = ncomp_fab++;
    int tmp_shear_comp = ncomp_fab++;
    int curvK_comp = ncomp_fab++;
    int curvP_comp = ncomp_fab++;
    int FXK_comp = ncomp_fab++;
    int FXP_comp = ncomp_fab++;
    int FEK_comp = ncomp_fab++;
    int FEP_comp = ncomp_fab++;
    int FCK_comp = ncomp_fab++;
    int FCP_comp = ncomp_fab++;
    int BCK_comp = ncomp_fab++;
    int BCP_comp = ncomp_fab++;

    for ( MFIter mfi(*mf_gls, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Box  bx = mfi.tilebox();
        Box xbx = surroundingNodes(bx,0);
        Box ybx = surroundingNodes(bx,1);
        Box gbx1 = grow(bx,IntVect(NGROW-1,NGROW-1,0));

        Box bx_rho = bx;
        bx_rho.convert(IntVect(0,0,0));
        Box bx_growloxy = growLo(growLo(grow(bx,IntVect(0,0,-1)),0,1),1,1);

        Box bxD = bx;
        bxD.makeSlab(2,0);
        Box gbx1D = gbx1;
        gbx1D.makeSlab(2,0);

        int ncompbc = 1;
        Vector<BCRec> bcrs_x(ncompbc);
        Vector<BCRec> bcrs_y(ncompbc);
        amrex::setBC(xbx,domain,xvel_bc(),0,1,domain_bcs_type,bcrs_x);
        amrex::setBC(ybx,domain,yvel_bc(),0,1,domain_bcs_type,bcrs_y);

        Array4<Real const> const& W = mf_W.const_array(mfi);
        Array4<Real> const& Hz = vec_Hz[lev]->array(mfi);
        Array4<Real> const& pm = vec_pm[lev]->array(mfi);
        Array4<Real> const& pn = vec_pn[lev]->array(mfi);
        Array4<Real> const& Lscale = vec_Lscale[lev]->array(mfi);

        Array4<Real> const& Huon = vec_Huon[lev]->array(mfi);
        Array4<Real> const& Hvom = vec_Hvom[lev]->array(mfi);
        Array4<Real> const& z_w = vec_z_w[lev]->array(mfi);

        Array4<Real> const& tke = mf_tke->array(mfi);
        Array4<Real> const& gls = mf_gls->array(mfi);

        Array4<Real const> const& sustr = vec_sustr[lev]->const_array(mfi);
        Array4<Real const> const& svstr = vec_svstr[lev]->const_array(mfi);
        Array4<Real const> const& bustr = vec_bustr[lev]->const_array(mfi);
        Array4<Real const> const& bvstr = vec_bvstr[lev]->const_array(mfi);
        Array4<Real const> const& msku = mf_msku->const_array(mfi);
        Array4<Real const> const& mskv = mf_mskv->const_array(mfi);

        Array4<Real> const& ZoBot = vec_ZoBot[lev]->array(mfi);

        FArrayBox fab(gbx1,ncomp_fab, amrex::The_Async_Arena()); fab.template setVal<RunOn::Device>(zero);

        auto CF = mf_w.array(mfi,CF_comp);
        auto shear2 = mf.array(mfi,shear2_comp);
        auto shear2_cached = mf.array(mfi,shear2_cache_comp);
        auto buoy2 = mf.array(mfi,buoy2_comp);
        Array4<Real> const& bvf = vec_bvf[lev]->array(mfi);

        auto tmp_buoy = fab.array(tmp_buoy_comp);
        auto tmp_shear = fab.array(tmp_shear_comp);
        auto curvK = fab.array(curvK_comp);
        auto curvP = fab.array(curvP_comp);
        auto FXK = fab.array(FXK_comp);
        auto FXP = fab.array(FXP_comp);
        auto FEK = fab.array(FEK_comp);
        auto FEP = fab.array(FEP_comp);
        auto FCK = fab.array(FCK_comp);
        auto FCP = fab.array(FCP_comp);
        auto BCK = fab.array(BCK_comp);
        auto BCP = fab.array(BCP_comp);

        auto Akt = mf_Akt->array(mfi);
        auto Akv = mf_Akv->array(mfi);
        auto Akp = mf_Akp->array(mfi);
        auto Akk = mf_Akk->array(mfi);

        ParallelFor(bx_growloxy, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            tmp_buoy(i,j,k)=Real(0.25) * (bvf(i,j,k) + bvf(i+1,j,k) + bvf(i,j+1,k)+bvf(i+1,j+1,k));
            tmp_shear(i,j,k)=Real(0.25) * (shear2_cached(i,j,k) + shear2_cached(i+1,j,k) + shear2_cached(i,j+1,k)+shear2_cached(i+1,j+1,k));
        });

        ParallelFor(grow(bx,IntVect(0,0,-1)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            buoy2(i,j,k)=Real(0.25) * (tmp_buoy(i,j,k) + tmp_buoy(i-1,j,k) + tmp_buoy(i,j-1,k)+tmp_buoy(i-1,j-1,k));
            shear2(i,j,k)=Real(0.25) * (tmp_shear(i,j,k) + tmp_shear(i-1,j,k) + tmp_shear(i,j-1,k)+tmp_shear(i-1,j-1,k));
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
            Real cff = Real(0.5) * (Huon(i,j,k) + Huon(i,j,k-1));
            Real cff1 = (cff > zero) ? curvK(i-1,j,k) : curvK(i,j,k);
            Real cff2 = (cff > zero) ? curvP(i-1,j,k) : curvP(i,j,k);

            FXK(i,j,k) = cff * Real(0.5) * (tke(i-1,j,k,2)+tke(i,j,k,2)-Gadv*cff1);
            FXP(i,j,k) = cff * Real(0.5) * (gls(i-1,j,k,2)+gls(i,j,k,2)-Gadv*cff2);
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
            Real cff = Real(0.5) * (Hvom(i,j,k) + Hvom(i,j,k-1));
            Real cff1 = (cff > zero) ? curvK(i,j-1,k) : curvK(i,j,k);
            Real cff2 = (cff > zero) ? curvP(i,j-1,k) : curvP(i,j,k);

            FEK(i,j,k) = cff * Real(0.5) * (tke(i,j-1,k,2)+tke(i,j,k,2)-Gadv*cff1);
            FEP(i,j,k) = cff * Real(0.5) * (gls(i,j-1,k,2)+gls(i,j,k,2)-Gadv*cff2);
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
            Real cff1 = Real(7.0) / Real(12.0);
            Real cff2 = one / Real(12.0);
            for (int k=1; k<=N-1; k++) {
                Real cff = Real(0.5) * (W(i,j,k+1)+W(i,j,k));
                FCK(i,j,k) = cff * (cff1 * (tke(i,j,k  ,2)+tke(i,j,k+1,2))-
                                    cff2 * (tke(i,j,k-1,2)+tke(i,j,k+2,2)));
                FCP(i,j,k) = cff * (cff1 * (gls(i,j,k  ,2)+gls(i,j,k+1,2))-
                                    cff2 * (gls(i,j,k-1,2)+gls(i,j,k+2,2)));
            }
            cff1 = one/Real(3.0);
            cff2 = Real(5.0)/Real(6.0);
            Real cff3 = one / Real(6.0);
            Real cff = Real(0.5) * (W(i,j,0)+W(i,j,1));
            FCK(i,j,0) = cff * (cff1 * tke(i,j,0,2)+cff2 * tke(i,j,1,2)-cff3 * tke(i,j,2,2));
            FCP(i,j,0) = cff * (cff1 * gls(i,j,0,2)+cff2 * gls(i,j,1,2)-cff3 * gls(i,j,2,2));
            cff = Real(0.5) * (W(i,j,N+1)+W(i,j,N));
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
        Real cff = -Real(0.5) * dt_lev;
        ParallelFor(convert(bx,IntVect(0,0,0)), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (k==0 or k==N) {
                FCK(i,j,k) = zero;
                FCP(i,j,k) = zero;
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
            Real gls_c3 = (strat2 > zero) ? gls_c3m : gls_c3p;
            Real Kprod = shear2(i,j,k) * (Akv(i,j,k)-Akv_bak) -
                         strat2 * (Akt(i,j,k,Temp_comp)-Akt_bak[Temp_comp]);
            Real Pprod = gls_c1 * shear2(i,j,k) * (Akv(i,j,k)-Akv_bak) -
                         gls_c3 * strat2 * (Akt(i,j,k,Temp_comp)-Akt_bak[Temp_comp]);

            // If negative production terms, then add buoyancy to dissipation terms
            // (BCK and BCP) below, using "cff1" and "cff2" as the on/off switch.
            Real cff1 = (Kprod < zero) ? zero : one;
            Real cff2 = (Pprod < zero) ? zero : one;
            Kprod = (Kprod < zero) ? Kprod + strat2*(Akt(i,j,k,Temp_comp)-Akt_bak[Temp_comp]) : Kprod;
            Pprod = (Pprod < zero) ? Pprod + gls_c3*strat2*(Akt(i,j,k,Temp_comp)-Akt_bak[Temp_comp]) : Pprod;
            // Time-step shear and buoyancy production terms.
            Real cff_Hz = Real(0.5) * (Hz(i,j,k) + Hz(i,j,k-1));
            tke(i,j,k,nnew) = tke(i,j,k,nnew)+dt_lev * cff_Hz * Kprod;
            gls(i,j,k,nnew) = gls(i,j,k,nnew)+dt_lev
                                *cff_Hz*Pprod*gls(i,j,k,nstp) / std::max(tke(i,j,k,nstp),gls_Kmin);

            Real gls_exp_exp1 = std::pow(gls(i,j,k,nstp),gls_exp1);
            Real gls_exp_mexp1 = one / (gls_exp_exp1);
            Real tke_exp_mexp1 = std::pow(tke(i,j,k,nstp),-tke_exp1);
            Real tke_exp_exp2 = std::pow(tke(i,j,k,nstp),tke_exp2);

            // Compute dissipation of turbulent energy (m3/s3).
            Real wall_fac = one;
            if (Lmy25) {
                wall_fac=one+gls_E2/(vonKar*vonKar)*
                        std::pow(gls_exp_exp1*cmu_fac1*
                         tke_exp_mexp1*
                         (one/ (z_w(i,j,k)-z_w(i,j,0))),2)+
                        Real(0.25)/(vonKar*vonKar)*
                        std::pow(gls_exp_exp1*cmu_fac1*
                         tke_exp_mexp1*
                         (one/ (z_w(i,j,N+1)-z_w(i,j,k))),2);
            }
            BCK(i,j,k)=cff_Hz*(one+dt_lev*
                          gls_exp_mexp1*cmu_fac2*
                          tke_exp_exp2+
                          dt_lev*(one-cff1)*strat2*
                          (Akt(i,j,k,Temp_comp)-Akt_bak[Temp_comp])/
                          tke(i,j,k,nstp))-
                          FCK(i,j,k)-FCK(i,j,k-1);
            BCP(i,j,k)=cff_Hz*(one+dt_lev*gls_c2*wall_fac*
                          gls_exp_mexp1*cmu_fac2*
                          tke_exp_exp2+
                          dt_lev*(one-cff2)*gls_c3*strat2*
                          (Akt(i,j,k,Temp_comp)-Akt_bak[Temp_comp])/
                          tke(i,j,k,nstp))-
                          FCP(i,j,k)-FCP(i,j,k-1);
        });

        // Compute production and dissipation terms.
        ParallelFor(bxD, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            Real Zob_min = std::max(ZoBot(i,j,0), Real(0.0001));
            //----------------------------------------------------------------------
            // Time-step dissipation and vertical diffusion terms implicitly.
            //----------------------------------------------------------------------
            //
            // Set Dirichlet surface and bottom boundary conditions. Compute
            // surface roughness from wind stress (Charnok) and set Craig and
            // Banner wave breaking surface flux, if appropriate.


            tke(i,j,N+1,nnew)=std::max(cmu_fac3*Real(0.5)*
                                     std::sqrt((sustr(i,j,0)+sustr(i+1,j,0))*(sustr(i,j,0)+sustr(i+1,j,0))+
                                          (svstr(i,j,0)+svstr(i,j+1,0))*(svstr(i,j,0)+svstr(i,j+1,0))),
                                     gls_Kmin);
            tke(i,j,0,nnew)=std::max(cmu_fac3*Real(0.5)*
                                 std::sqrt((bustr(i,j,0)+bustr(i+1,j,0))*(bustr(i,j,0)+bustr(i+1,j,0))+
                                      (bvstr(i,j,0)+bvstr(i,j+1,0))*(bvstr(i,j,0)+bvstr(i,j+1,0))),
                                        gls_Kmin);

            gls(i,j,N+1,nnew)=std::max(cmu0_exp_p*
                                    std::pow(tke(i,j,N+1,nnew),gls_m)*
                                    std::pow(L_sft*Zos_eff,gls_n), gls_Pmin);
            Real cff_gls = gls_fac4*std::pow(vonKar*Zob_min,gls_n);
            gls(i,j,0,nnew)=std::max(cff_gls*std::pow(tke(i,j,0,nnew),(gls_m)), gls_Pmin);

            // Solve tri-diagonal system for turbulent kinetic energy.
            // Might be N instead of N-1?
            Real tke_fluxt = zero;
            Real tke_fluxb = zero;
            Real cff_BCK = one/BCK(i,j,N);
            CF(i,j,N)=cff_BCK*FCK(i,j,N-1);
            tke(i,j,N,nnew)=cff_BCK*(tke(i,j,N,nnew)+tke_fluxt);
            for (int k=N-1;k>=1;k--) {
                cff_BCK = one / (BCK(i,j,k)-CF(i,j,k+1)*FCK(i,j,k));
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
            Real cff_tke = Real(0.5) * (tke(i,j,N+1,nnew) + tke(i,j,N,nnew));
            Real gls_fluxt = dt_lev*gls_fac3*std::pow(cff_tke,gls_m)*
                             std::pow(L_sft,(gls_n))*
                             std::pow(Zos_eff+Real(0.5)*Hz(i,j,N),gls_n-one)*
                             Real(0.5)*(Akp(i,j,N+1)+Akp(i,j,N));
            cff_tke=Real(0.5)*(tke(i,j,0,nnew)+tke(i,j,1,nnew));
            Real gls_fluxb = dt_lev*gls_fac2*std::pow(cff_tke,gls_m)*
                              std::pow(Real(0.5)*Hz(i,j,0)+Zob_min,gls_n-one)*
                              Real(0.5)*(Akp(i,j,0)+Akp(i,j,1));
            Real cff_BCP = one / BCP(i,j,N);
            CF(i,j,N) = cff_BCP * FCP(i,j,N-1);
            gls(i,j,N,nnew)=cff_BCP*(gls(i,j,N,nnew)-gls_fluxt);
            for (int k=N-1;k>=1;k--) {
                cff_BCP = one / (BCP(i,j,k)-CF(i,j,k+1)*FCP(i,j,k));
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
            Real gls_comparison = gls_fac5 *
                                    std::pow(tke(i,j,k,nnew),tke_exp4)*
                                    std::pow(std::sqrt(std::max(Real(0.0),
                                          buoy2(i,j,k)))+eps,-gls_n);
            gls(i,j,k,nnew) = (gls_n >= Real(0.0)) ? std::min(gls(i,j,k,nnew),gls_comparison) : std::max(gls(i,j,k,nnew),gls_comparison);
            Real Ls_lmt;
            Real Ls_unlmt=std::max(eps,
                                   std::pow(gls(i,j,k,nnew),( gls_exp1))*cmu_fac1*
                                   std::pow(tke(i,j,k,nnew),(-tke_exp1)));
            // Some problems are very sensitive to this condition (ultimate cause of
            // some discrepancies in BoundaryLayer test between CPU and GPU)
            Ls_lmt = (buoy2(i,j,k) > Real(0.0)) ? std::min(Ls_unlmt,
                                                std::sqrt(Real(0.56)*tke(i,j,k,nnew)/
                                                (std::max(Real(0.0),buoy2(i,j,k))+eps))) : Ls_unlmt;
            //
            //  Recompute gls based on limited length scale
            //
            gls(i,j,k,nnew)=std::max(cmu0_exp_p*
                                           std::pow(tke(i,j,k,nnew),gls_m)*
                                           std::pow(Ls_lmt,gls_n), gls_Pmin);

            //   Compute nondimensional stability functions for tracers (Sh) and
            //   momentum (Sm).
            Real Sh, Sm;
            Real Gh=std::min(gls_Gh0,-buoy2(i,j,k)*Ls_lmt*Ls_lmt/
                            (two*tke(i,j,k,nnew)));
            Gh=std::min(Gh,Gh-(Gh-gls_Ghcri)*(Gh-gls_Ghcri)/
                       (Gh+gls_Gh0-two*gls_Ghcri));
            Gh=std::max(Gh,gls_Ghmin);

            if (gls_stability_type == GLS_StabilityType::Canuto_A ||
                gls_stability_type == GLS_StabilityType::Canuto_B) {
                //
                //   Canuto stability: Compute shear number.
                //
                Real Gm=(gls_b0/gls_fac6-gls_b1*Gh+gls_b3*gls_fac6*(Gh*Gh))/
                             (gls_b2-gls_b4*gls_fac6*Gh);
                Gm=std::min(Gm,shear2(i,j,k)*Ls_lmt*Ls_lmt/
                            (two*tke(i,j,k,nnew)));
                /////Gm=std::min(Gm,(gls_s1*gls_fac6*Gh-gls_s0)/(gls_s2*gls_fac6));
                //
                //  Compute stability functions
                //
                Real stab_cff=gls_b0-gls_b1*gls_fac6*Gh+gls_b2*gls_fac6*Gm+
                    gls_b3*gls_fac6*gls_fac6*Gh*Gh-gls_b4*gls_fac6*gls_fac6*Gh*Gm+
                    gls_b5*gls_fac6*gls_fac6*Gm*Gm;
                Sm=(gls_s0-gls_s1*gls_fac6*Gh+gls_s2*gls_fac6*Gm)/stab_cff;
                Sh=(gls_s4-gls_s5*gls_fac6*Gh+gls_s6*gls_fac6*Gm)/stab_cff;
                Sm=std::max(Sm,Real(0.0));
                Sh=std::max(Sh,Real(0.0));

                //
                //  Relate Canuto stability to ROMS notation
                //
                Sm=Sm*sqrt2/(gls_cmu0_cube);
                Sh=Sh*sqrt2/gls_cmu0_cube;
            } else if (gls_stability_type == GLS_StabilityType::Galperin) {
                Real cff_galperin = one - my_Sh2*Gh;
                Sh = my_Sh1 / cff_galperin;
                Sm = (my_Sm3+Sh*Gh*my_Sm4)/(one-my_Sm2*Gh);
            } else {
                Sh = zero;
                Sm = zero;
            }

            //  Compute vertical mixing (m2/s) coefficients of momentum and
            //  tracers.  Average ql over the two timesteps rather than using
            //  the new Lscale and just averaging tke.

            Real ql=sqrt2*Real(0.5)*(Ls_lmt*std::sqrt(tke(i,j,k,nnew))+
                                  Lscale(i,j,k)*std::sqrt(tke(i,j,k,nstp)));
            Akv(i,j,k)=Akv_bak+Sm*ql;
            for (int n=0; n<NAT; n++) {
                Akt(i,j,k,n)=Akt_bak[n]+Sh*ql;
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
            Real Zob_min = std::max(ZoBot(i,j,0), Real(0.0001));
            Akv(i,j,N+1)=Akv_bak+L_sft*Zos_eff*gls_cmu0*
                          std::sqrt(tke(i,j,N+1,nnew));
            Akv(i,j,0)=Akv_bak+vonKar*Zob_min*gls_cmu0*
                      std::sqrt(tke(i,j,0,nnew));

            Akk(i,j,N+1)=Akk_bak+Akv(i,j,N+1)/gls_sigk;
            Akk(i,j,0)=Akk_bak+Akv(i,j,0)/gls_sigk;
            Akp(i,j,N+1)=Akp_bak+Akv(i,j,N+1)*ogls_sigp;
            Akp(i,j,0)=Akp_bak+Akv(i,j,0)/gls_sigp_cb;

            for (int n=0; n<NAT; n++) {
                Akt(i,j,N+1,n)  = Akt_bak[n];
                Akt(i,j,0,n) = Akt_bak[n];
            }
        });
    }

    for (int icomp=0; icomp<3; icomp++) {
        FillPatch(lev, t_old[lev], *mf_tke, GetVecOfPtrs(vec_tke), zvel_bc(), BdyVars::null, icomp, false, false);
        FillPatch(lev, t_old[lev], *mf_gls, GetVecOfPtrs(vec_gls), zvel_bc(), BdyVars::null, icomp, false, false);
    }
    for (int icomp=0; icomp<NAT; icomp++) {
        FillPatch(lev, t_old[lev], *mf_Akt, GetVecOfPtrs(vec_Akt), zvel_bc(), BdyVars::null, icomp, false, false);
    }
    FillPatchNoBC(lev, t_old[lev], *mf_Akv, GetVecOfPtrs(vec_Akv), BdyVars::null);
    FillPatchNoBC(lev, t_old[lev], *mf_Akp, GetVecOfPtrs(vec_Akp), BdyVars::null);
    FillPatchNoBC(lev, t_old[lev], *mf_Akk, GetVecOfPtrs(vec_Akk), BdyVars::null);
}
