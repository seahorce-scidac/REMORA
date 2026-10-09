#include <REMORA.H>

using namespace amrex;
using namespace REMORASeaIce;

// The ice-ocean coupling and thermodynamics of the myROMS ice model:
// ice_spdiw.F, ice_vbc.F (Mellor-Kantha branch), ice_mk.h and ice_limit.F.
//
// Index offsets: REMORA cell i is ROMS rho point i+1, REMORA face i is ROMS
// u (or v) point i+1, REMORA level k is ROMS level k+1, and REMORA z_w(k) is
// ROMS z_w(k). ROMS loops over a tile's interior (Istr..Iend) become loops
// over the valid box; boundary rows and faces then come from the BC routines,
// after an exchange, as described in REMORA_SeaIce_bc.cpp.

/**
 * ROMS ice_spdiw: velocity averaged over the top 10 m (ROMS uses the KPP
 * boundary layer depth when it is deeper; REMORA has no KPP) and the speed of
 * the ice relative to it. The ice velocity is read at liuol, before the EVP
 * substeps swap the velocity levels: from the second step on that is the
 * velocity of two steps back, as in ROMS.
 */
void
REMORA::ice_spdiw (int lev)
{
    BL_PROFILE("REMORA::ice_spdiw()");

    const int nz = geom[lev].Domain().length(2);
    const int liuol = ice_lvl.liuol;
    const int cU = ucomp(IceU::Uice, liuol);

    MultiFab& fir = *vec_ice_fi_r[lev];
    MultiFab& fiu = *vec_ice_fi_u[lev];
    MultiFab& fiv = *vec_ice_fi_v[lev];

    for (MFIter mfi(fir, false); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        Array4<Real const> const& z_r = vec_z_r[lev]->const_array(mfi);
        Array4<Real const> const& z_w = vec_z_w[lev]->const_array(mfi);
        Array4<Real const> const& u   = xvel_old[lev]->const_array(mfi);
        Array4<Real const> const& v   = yvel_old[lev]->const_array(mfi);
        Array4<Real const> const& ui  = vec_ice_si_u[lev]->const_array(mfi, cU);
        Array4<Real const> const& vi  = vec_ice_si_v[lev]->const_array(mfi, cU);
        Array4<Real      > const& fr  = fir.array(mfi);
        Array4<Real      > const& fu  = fiu.array(mfi);
        Array4<Real      > const& fv  = fiv.array(mfi);

        // Level (ROMS 1-based) just below the mixed layer depth
        const Box nbx = amrex::grow(bx, IntVect(2,2,0));
        IArrayBox nbot_fab(nbx, 1, The_Async_Arena());
        Array4<int> const& nbot = nbot_fab.array();
        ParallelFor(nbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real mlio = Real(-10.0);
            int nb = 1;
            for (int k = nz; k >= 1; --k) {
                if (z_r(i,j,k-1) < mlio) {
                    nb = std::max(std::min(k, nz), 1);
                    break;
                }
            }
            nbot(i,j,0) = nb;
        });

        // Mixed-layer averaged velocities on every face of the box
        const Box xbx = amrex::surroundingNodes(bx, 0);
        const Box ybx = amrex::surroundingNodes(bx, 1);
        FArrayBox uw_fab(xbx, 1, The_Async_Arena());
        FArrayBox vw_fab(ybx, 1, The_Async_Arena());
        Array4<Real> const& uw = uw_fab.array();
        Array4<Real> const& vw = vw_fab.array();

        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            int nbotu = static_cast<int>(std::lround(Real(0.5)*(nbot(i-1,j,0)+nbot(i,j,0))));
            nbotu = std::max(std::min(nbotu, nz), 1);
            Real sum = zero;
            Real totml = zero;
            for (int k = nz; k >= nbotu; --k) {
                const Real dml = Real(0.5)*(z_w(i-1,j,k)-z_w(i-1,j,k-1)+
                                            z_w(i  ,j,k)-z_w(i  ,j,k-1));
                sum += u(i,j,k-1)*dml;
                totml += dml;
            }
            uw(i,j,0) = sum/totml;
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            int nbotv = static_cast<int>(std::lround(Real(0.5)*(nbot(i,j-1,0)+nbot(i,j,0))));
            nbotv = std::max(std::min(nbotv, nz), 1);
            Real sum = zero;
            Real totml = zero;
            for (int k = nz; k >= nbotv; --k) {
                const Real dml = Real(0.5)*(z_w(i,j-1,k)-z_w(i,j-1,k-1)+
                                            z_w(i,j  ,k)-z_w(i,j  ,k-1));
                sum += v(i,j,k-1)*dml;
                totml += dml;
            }
            vw(i,j,0) = sum/totml;
        });

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real du = uw(i,j,0)-ui(i,j,0)+uw(i+1,j,0)-ui(i+1,j,0);
            const Real dv = vw(i,j,0)-vi(i,j,0)+vw(i,j+1,0)-vi(i,j+1,0);
            fr(i,j,0,IceFR::IOvs) = Real(0.5)*std::sqrt(du*du+dv*dv);
        });
        ParallelFor(interior_faces(xbx, 0, geom[lev]), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            fu(i,j,0,IceFU::Uavg) = uw(i,j,0);
        });
        ParallelFor(interior_faces(ybx, 1, geom[lev]), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            fv(i,j,0,IceFU::Uavg) = vw(i,j,0);
        });
    }

    ice_exchange(lev, fir, IceFR::IOvs, 1);
    ice_exchange(lev, fiu, IceFU::Uavg, 1);
    ice_exchange(lev, fiv, IceFU::Uavg, 1);
    ice_bc_r2d(lev, fir, IceFR::IOvs);
    ice_bc_u2d(lev, fiu, IceFU::Uavg);
    ice_bc_v2d(lev, fiv, IceFU::Uavg);
}

/**
 * ROMS ice_vbc, Mellor-Kantha branch: the ice-water transfer coefficient, the
 * air-ice stress on the ice, and the ocean surface stress, which becomes the
 * ice-water stress over the ice-covered fraction and the bulk-flux air-sea
 * stress over the rest. The ice velocity is at liuol, the ice cover at liold.
 */
void
REMORA::ice_vbc (int lev)
{
    BL_PROFILE("REMORA::ice_vbc()");

    const int nz = geom[lev].Domain().length(2);
    const int liold = ice_lvl.liold;
    const int liuol = ice_lvl.liuol;
    const int cA = rcomp(IceR::Aice, liold);
    const int cH = rcomp(IceR::Hice, liold);
    const int cU = ucomp(IceU::Uice, liuol);
    const bool mom_bulk = seaice_params.mom_bulk;
    const Real Cd_io = seaice_params.Cd_io;
    const Real IceRho = seaice_params.IceRho;
    const Real min_ai = seaice_params.min_ai;
    constexpr Real kappa = Real(0.4);
    constexpr Real z0ii = Real(0.01);
    // ROMS writes 1.0E-20 without _r8, a default (single precision) real
    const Real eps = static_cast<Real>(1.0e-20f);
    const Real pi = PI;

    MultiFab& fir = *vec_ice_fi_r[lev];
    MultiFab& fiu = *vec_ice_fi_u[lev];
    MultiFab& fiv = *vec_ice_fi_v[lev];

    for (MFIter mfi(fir, false); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        Array4<Real const> const& si   = vec_ice_si_r[lev]->const_array(mfi);
        Array4<Real const> const& ui   = vec_ice_si_u[lev]->const_array(mfi, cU);
        Array4<Real const> const& vi   = vec_ice_si_v[lev]->const_array(mfi, cU);
        Array4<Real const> const& z_r  = vec_z_r[lev]->const_array(mfi);
        Array4<Real const> const& z_w  = vec_z_w[lev]->const_array(mfi);
        Array4<Real const> const& Zt_avg1 = vec_Zt_avg1[lev]->const_array(mfi);
        Array4<Real      > const& fr   = fir.array(mfi);
        Array4<Real      > const& fu   = fiu.array(mfi);
        Array4<Real      > const& fv   = fiv.array(mfi);
        Array4<Real      > const& sustr = vec_sustr[lev]->array(mfi);
        Array4<Real      > const& svstr = vec_svstr[lev]->array(mfi);

        // ROMS loops over Istr-1..Iend, Jstr-1..Jend; one more cell is harmless
        const Box cbx = amrex::grow(bx, IntVect(1,1,0));
        FArrayBox chuiw_fab(cbx, 3, The_Async_Arena());
        Array4<Real> const& chuiw  = chuiw_fab.array(0);
        Array4<Real> const& utauiw = chuiw_fab.array(1);
        Array4<Real> const& chuai  = chuiw_fab.array(2);

        ParallelFor(cbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real spd = std::max(fr(i,j,0,IceFR::IOvs), Real(0.15));
            if (mom_bulk) {
                utauiw(i,j,0) = spd;
                chuiw(i,j,0) = Cd_io*spd;
                chuai(i,j,0) = zero;
            } else {
                const Real thic = si(i,j,0,cH)/std::max(si(i,j,0,cA), std::max(min_ai, eps));
                Real z0 = std::max(z0ii*thic, Real(0.01));
                z0 = std::min(z0, Real(0.1));
                const Real dztop = z_w(i,j,nz)-z_r(i,j,nz-1);
                Real zdz0 = dztop/z0;
                if (zdz0 < Real(6.0)) { zdz0 = Real(6.0); }
                utauiw(i,j,0) = std::sqrt(fr(i,j,0,IceFR::IOmt)*spd);
                utauiw(i,j,0) = std::max(utauiw(i,j,0), Real(1.0e-4));
                chuiw(i,j,0) = kappa*utauiw(i,j,0)/std::log(zdz0);
                chuai(i,j,0) = one-std::cos(one*pi*std::min(thic+Real(0.05), one));
            }
        });

        ParallelFor(interior_faces(amrex::surroundingNodes(bx,0), 0, geom[lev]),
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real aix  = Real(0.5)*(si(i,j,0,cA)+si(i-1,j,0,cA));
            const Real chux = Real(0.5)*(chuiw(i,j,0)+chuiw(i-1,j,0));
            const Real tau_ai = fr(i,j,0,IceFR::sustr_ai)+fr(i-1,j,0,IceFR::sustr_ai);
            if (mom_bulk) {
                fu(i,j,0,IceFU::AIus) = Real(0.5)*aix*tau_ai/IceRho;
            } else {
                const Real chuax = Real(0.5)*(chuai(i,j,0)+chuai(i-1,j,0));
                fu(i,j,0,IceFU::AIus) = Real(0.5)*aix*chuax*tau_ai/IceRho;
            }
            sustr(i,j,0) = aix*chux*(ui(i,j,0)-fu(i,j,0,IceFU::Uavg))+
                           (one-aix)*fu(i,j,0,IceFU::sustr_ao);
        });
        ParallelFor(interior_faces(amrex::surroundingNodes(bx,1), 1, geom[lev]),
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real aiy  = Real(0.5)*(si(i,j,0,cA)+si(i,j-1,0,cA));
            const Real chuy = Real(0.5)*(chuiw(i,j,0)+chuiw(i,j-1,0));
            const Real tau_ai = fr(i,j,0,IceFR::svstr_ai)+fr(i,j-1,0,IceFR::svstr_ai);
            if (mom_bulk) {
                fv(i,j,0,IceFU::AIus) = Real(0.5)*aiy*tau_ai/IceRho;
            } else {
                const Real chuay = Real(0.5)*(chuai(i,j,0)+chuai(i,j-1,0));
                fv(i,j,0,IceFU::AIus) = Real(0.5)*aiy*chuay*tau_ai/IceRho;
            }
            svstr(i,j,0) = aiy*chuy*(vi(i,j,0)-fv(i,j,0,IceFU::Uavg))+
                           (one-aiy)*fv(i,j,0,IceFU::sustr_ao);
        });

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            fr(i,j,0,IceFR::Hsse) = Zt_avg1(i,j,0);
            fr(i,j,0,IceFR::IOmt) = chuiw(i,j,0);
            fr(i,j,0,IceFR::IOfv) = utauiw(i,j,0);
        });
    }

    ice_exchange(lev, fir, 0, fir.nComp());
    ice_exchange(lev, fiu, IceFU::AIus, 1);
    ice_exchange(lev, fiv, IceFU::AIus, 1);
    vec_sustr[lev]->FillBoundary(geom[lev].periodicity());
    vec_svstr[lev]->FillBoundary(geom[lev].periodicity());

    ice_bc_r2d(lev, fir, IceFR::Hsse);
    ice_bc_u2d(lev, fiu, IceFU::AIus);
    ice_bc_v2d(lev, fiv, IceFU::AIus);
    ice_bc_r2d(lev, fir, IceFR::IOmt);
    ice_bc_r2d(lev, fir, IceFR::IOfv);
    ice_bc_u2d(lev, *vec_sustr[lev], 0);
    ice_bc_v2d(lev, *vec_svstr[lev], 0);
}

/**
 * ROMS ice_mk.h: Mellor and Kantha (1989) ice thermodynamics, on level linew
 * in place. Updates the ice surface and interior temperatures, the melt-pond
 * and snow budgets, the molecular sublayer under the ice, the ice thickness,
 * concentration and age, and sets the ocean surface heat flux (overwriting
 * it) and salt flux (adding to it).
 */
void
REMORA::ice_thermo (int lev, Real dt_lev)
{
    BL_PROFILE("REMORA::ice_thermo()");

    const int nz = geom[lev].Domain().length(2);
    const int kN = nz-1;
    const int linew = ice_lvl.linew;
    const int cA = rcomp(IceR::Aice, linew);
    const int cH = rcomp(IceR::Hice, linew);
    const int cS = rcomp(IceR::Hsno, linew);
    const int cM = rcomp(IceR::Hmel, linew);
    const int cT = rcomp(IceR::Tice, linew);
    const int cG = rcomp(IceR::Iage, linew);

    const Real dtice = dt_lev;
    const Real min_ai = seaice_params.min_ai;
    const Real max_ai = seaice_params.max_ai;
    const Real min_hi = seaice_params.min_hi;
    const Real max_hmelt = seaice_params.max_hmelt;
    const Real IceRho = seaice_params.IceRho;
    const Real SnowDryRho = seaice_params.SnowDryRho;
    const Real SnowWetRho = seaice_params.SnowWetRho;
    const bool convsnow = seaice_params.convsnow;

    constexpr Real AlphIc = Real(2.034);         // [W m-1 K-1]
    constexpr Real AlphSn = Real(0.31);          // [W m-1 K-1]
    constexpr Real Cp_i = Real(2093.0);          // [J kg-1 K-1]
    constexpr Real Cp_w = Real(3990.0);          // [J kg-1 K-1]
    constexpr Real eps = Real(1.0e-4);           // zero division
    constexpr Real frln = Real(-0.0543);         // [psu C-1]
    constexpr Real hfus = Real(3.347e5);         // [J kg-1]
    constexpr Real kappa = Real(0.4);            // von Karman
    constexpr Real nu = Real(1.8e-6);            // m2/s
    constexpr Real prs = Real(2432.0);           // S Schmidt number
    constexpr Real prt = Real(13.0);             // T Prandtl number
    constexpr Real RhoCpr = Real(0.2442754e-6);  // [m s2 K kg-1]
    constexpr Real RhoSW = Real(1026.0);         // [kg m-3]
    constexpr Real sice_ref = Real(3.2);         // [psu]
    constexpr Real tpr = Real(0.85);             // turbulent Prandtl number
    const Real ykf = static_cast<Real>(3.14f);   // Yaglom/Kader; ROMS writes 3.14 without _r8
    constexpr Real z0ii = Real(0.02);            // ice roughness
    constexpr Real sec2day = Real(1.0)/Real(86400.0);

    MultiFab& sir = *vec_ice_si_r[lev];
    MultiFab& fir = *vec_ice_fi_r[lev];

    for (MFIter mfi(sir, false); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        Array4<Real      > const& si    = sir.array(mfi);
        Array4<Real      > const& fi    = fir.array(mfi);
        Array4<Real const> const& t     = cons_old[lev]->const_array(mfi);
        Array4<Real const> const& sustr = vec_sustr[lev]->const_array(mfi);
        Array4<Real const> const& svstr = vec_svstr[lev]->const_array(mfi);
        Array4<Real const> const& rain  = vec_rain[lev]->const_array(mfi);
        Array4<Real const> const& z_r   = vec_z_r[lev]->const_array(mfi);
        Array4<Real const> const& z_w   = vec_z_w[lev]->const_array(mfi);
        Array4<Real const> const& rmask = vec_mskr[lev]->const_array(mfi);
        Array4<Real      > const& stflx = vec_stflx[lev]->array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            // Ocean surface temperature and salinity
            const Real temp_top = t(i,j,kN,Temp_comp);
            Real salt_top = t(i,j,kN,Salt_comp);
            salt_top = std::min(std::max(zero, salt_top), Real(40.0));

            // Square root of the surface stress magnitude
            const Real tx = Real(0.5)*(sustr(i,j,0)+sustr(i+1,j,0));
            const Real ty = Real(0.5)*(svstr(i,j,0)+svstr(i,j+1,0));
            Real utau = std::sqrt(std::sqrt(tx*tx+ty*ty));
            utau = std::max(utau, Real(1.0e-4));

            // Snow and ice thickness, brine fraction, thermal conductivity
            const Real sice = std::min(sice_ref, salt_top);
            const Real ice_thick = Real(0.05)+si(i,j,0,cH)/(si(i,j,0,cA)+eps);
            const Real snow_thick = si(i,j,0,cS)/(si(i,j,0,cA)+eps);
            Real brnfr = frln*sice/(si(i,j,0,cT)-eps);
            brnfr = std::min(brnfr, Real(0.2));
            brnfr = std::max(brnfr, zero);
            Real alph = AlphIc*std::max(one-Real(1.2)*brnfr, Real(0.25));
            {
                const Real cff = (si(i,j,0,cH)/one)*(si(i,j,0,cH)/one);
                const Real corfac = one/(Real(0.5)*(one+std::exp(-cff)));
                alph = alph*corfac;
            }
            const Real Coa = Real(2.0)*alph*snow_thick/(AlphSn*ice_thick);

            // Temperature at the top of the ice
            const Real b2d = Real(2.0)*alph/(ice_thick*(one+Coa));
            fi(i,j,0,IceFR::Qcon) += b2d;
            if (si(i,j,0,cA) > min_ai) {
                fi(i,j,0,IceFR::Qrhs) += b2d*(si(i,j,0,cT)+Real(273.15));
                fi(i,j,0,IceFR::Isst) = (fi(i,j,0,IceFR::Qrhs)/fi(i,j,0,IceFR::Qcon))-Real(273.15);
                fi(i,j,0,IceFR::Isst) = std::min(std::max(fi(i,j,0,IceFR::Isst), Real(-45.0)), zero);
            } else {
                fi(i,j,0,IceFR::Isst) = temp_top;
            }

            // New interior ice temperature
            if (si(i,j,0,cA) > min_ai) {
                const Real tmin = std::min(si(i,j,0,cT), frln*sice_ref);
                const Real cot = -frln*sice*hfus/(tmin*tmin)+Cp_i;
                const Real cff1 = IceRho*cot*ice_thick*ice_thick;
                const Real cff2 = fi(i,j,0,IceFR::Isst)-(Real(2.0)+Coa)*si(i,j,0,cT);
                const Real cff3 = one+Coa;
                si(i,j,0,cT) = si(i,j,0,cT)+dtice*(Real(2.0)*alph/cff1*
                                                    (fi(i,j,0,IceFR::T0mk)+cff2/cff3));
                si(i,j,0,cT) = std::max(si(i,j,0,cT), Real(-35.0));
                si(i,j,0,cT) = std::min(si(i,j,0,cT), frln*sice_ref);
                si(i,j,0,cT) = std::min(si(i,j,0,cT),
                                        std::max(fi(i,j,0,IceFR::Isst), fi(i,j,0,IceFR::T0mk)));
            } else {
                si(i,j,0,cT) = temp_top;
            }

            // Associated heat fluxes
            Real Qi2 = zero;
            Real Qio = zero;
            Real Qai;
            if (si(i,j,0,cA) > min_ai) {
                const Real hicehinv = one/(Real(0.5)*ice_thick);
                const Real t2 = (fi(i,j,0,IceFR::Isst)+Coa*si(i,j,0,cT))/(one+Coa);
                Qi2 = alph*(si(i,j,0,cT)-t2)*hicehinv;
                Qio = alph*(fi(i,j,0,IceFR::T0mk)-si(i,j,0,cT))*hicehinv;
            }
            Qai = fi(i,j,0,IceFR::Qnet_ai);

            // Open water: no ice fluxes
            if (si(i,j,0,cA) <= min_ai) {
                fi(i,j,0,IceFR::Isst) = fi(i,j,0,IceFR::T0mk);
                si(i,j,0,cT) = Real(-2.0);
                fi(i,j,0,IceFR::Isst) *= rmask(i,j,0);
                si(i,j,0,cT) *= rmask(i,j,0);
                Qi2 = zero;
                Qai = zero;
                Qio = zero;
                si(i,j,0,cS) = zero;
                si(i,j,0,cM) = zero;
            }

            // Surface water accumulation: ice melting
            const Real ws = std::max(fi(i,j,0,IceFR::snow), zero);
            Real hfus1 = zero;
            {
                const Real tfrz = frln*sice;
                fi(i,j,0,IceFR::W_ai) = zero;
                fi(i,j,0,IceFR::W_ro) = zero;

                if (si(i,j,0,cA) > min_ai) {
                    const Real cff = one-brnfr;
                    hfus1 = hfus*cff+fi(i,j,0,IceFR::Isst)*Cp_w-
                            (cff*Cp_i+brnfr*Cp_w)*si(i,j,0,cT);
                    Qai = fi(i,j,0,IceFR::Qnet_ai);
                    Qi2 = b2d*(si(i,j,0,cT)-fi(i,j,0,IceFR::Isst));

                    Real Qsur;
                    if ((si(i,j,0,cS) <= eps) && (si(i,j,0,cM) <= eps)) {
                        Qsur = -(Qai-Qi2)/(hfus1*RhoSW);
                    } else if ((si(i,j,0,cS) <= eps) && (si(i,j,0,cM) > eps)) {
                        Qsur = -(Qai-Qi2)/(hfus1*Real(1003.1));
                    } else {
                        Qsur = -(Qai-Qi2)/(hfus*SnowWetRho);
                    }

                    if ((si(i,j,0,cS) > eps) && (fi(i,j,0,IceFR::Isst) >= zero)) {
                        si(i,j,0,cS) = si(i,j,0,cS)-si(i,j,0,cA)*std::max(Qsur, zero)*dtice;
                        si(i,j,0,cM) = si(i,j,0,cM)+si(i,j,0,cA)*std::max(Qsur, zero)*
                                       SnowWetRho/RhoSW*dtice;
                    } else if ((si(i,j,0,cM) > eps) && (fi(i,j,0,IceFR::Isst) <= tfrz)) {
                        fi(i,j,0,IceFR::W_ai) = std::min(Qsur, zero);
                        si(i,j,0,cM) = si(i,j,0,cM)+si(i,j,0,cA)*std::min(Qsur, zero)*dtice;
                    } else if ((si(i,j,0,cS) <= eps) && (si(i,j,0,cM) >= eps) &&
                               (fi(i,j,0,IceFR::Isst) > tfrz)) {
                        fi(i,j,0,IceFR::W_ai) = std::max(Qsur, zero);
                        si(i,j,0,cM) = si(i,j,0,cM)+si(i,j,0,cA)*std::max(Qsur, zero)*dtice;
                    } else if ((si(i,j,0,cS) < eps) && (si(i,j,0,cM) < eps) &&
                               (fi(i,j,0,IceFR::Isst) > tfrz)) {
                        fi(i,j,0,IceFR::W_ai) = std::max(Qsur, zero);
                        si(i,j,0,cM) = si(i,j,0,cM)+si(i,j,0,cA)*std::max(Qsur, zero)*dtice;
                    }

                    const Real rn = rain(i,j,0);
                    if (rn <= zero) {
                        si(i,j,0,cS) = si(i,j,0,cS)+si(i,j,0,cA)*ws*dtice;
                    } else if ((si(i,j,0,cS) > zero) && (si(i,j,0,cM) == zero)) {
                        si(i,j,0,cS) = std::max(zero, si(i,j,0,cS)-si(i,j,0,cA)*rn/SnowDryRho);
                        fi(i,j,0,IceFR::W_ai) = fi(i,j,0,IceFR::W_ai)-
                                                Real(2.0)*si(i,j,0,cA)*rn/IceRho;
                    } else if ((si(i,j,0,cS) > zero) && (si(i,j,0,cM) > zero)) {
                        si(i,j,0,cS) = std::max(zero, si(i,j,0,cS)-
                                                Real(0.5)*si(i,j,0,cA)*rn/SnowDryRho);
                        fi(i,j,0,IceFR::W_ai) = fi(i,j,0,IceFR::W_ai)-
                                                Real(0.5)*si(i,j,0,cA)*rn/IceRho;
                        si(i,j,0,cM) = si(i,j,0,cM)+si(i,j,0,cA)*Real(0.5)*rn/RhoSW*dtice;
                    } else if (si(i,j,0,cM) > zero) {
                        si(i,j,0,cM) = si(i,j,0,cM)+si(i,j,0,cA)*rn/RhoSW*dtice;
                    } else {
                        fi(i,j,0,IceFR::W_ai) = fi(i,j,0,IceFR::W_ai)-si(i,j,0,cA)*rn/IceRho;
                    }

                    // Limit the surface water by the smaller of a maximum and the ice thickness
                    const Real hmax = std::min(max_hmelt, si(i,j,0,cH));
                    if (si(i,j,0,cM) > hmax) {
                        fi(i,j,0,IceFR::W_ro) = (si(i,j,0,cM)-hmax)/dtice;
                        si(i,j,0,cM) = hmax;
                    }
                }
            }

            // Molecular sublayer under ice: Yaglom and Kader (1974) roughness scales
            Real cht, chs;
            {
                Real z0 = std::max(z0ii*ice_thick, Real(0.01));
                z0 = std::min(z0, Real(0.1));
                const Real dztop = z_w(i,j,nz)-z_r(i,j,kN);
                Real zdz0 = dztop/z0;
                zdz0 = std::max(zdz0, Real(3.0));
                const Real rno = utau*Real(0.09)/nu;
                const Real termt = ykf*std::sqrt(rno)*std::pow(prt, Real(0.666667));
                const Real terms = ykf*std::sqrt(rno)*std::pow(prs, Real(0.666667));
                cht = utau/(tpr*(std::log(zdz0)/kappa+termt));
                chs = utau/(tpr*(std::log(zdz0)/kappa+terms));
            }

            // Temperature and salinity of the molecular sublayer under ice
            const Real tfz = frln*salt_top;
            fi(i,j,0,IceFR::W_ao) = zero;
            fi(i,j,0,IceFR::W_io) = zero;
            const Real xwai = std::max(zero, fi(i,j,0,IceFR::W_ai));
            {
                const Real cff = one-brnfr;
                hfus1 = hfus*cff+fi(i,j,0,IceFR::T0mk)*Cp_w-(cff*Cp_i+brnfr*Cp_w)*si(i,j,0,cT);
            }
            const Real Qnet_ao = fi(i,j,0,IceFR::Qnet_ao);
            if (((temp_top <= tfz) && (Qnet_ao > zero)) ||
                ((temp_top >= tfz) && (Qnet_ao < zero) && (si(i,j,0,cA) > zero))) {
                fi(i,j,0,IceFR::W_ao) = Qnet_ao/(hfus1*RhoSW);
            }
            Real xtot;
            if ((si(i,j,0,cA) <= min_ai) || (si(i,j,0,cH) <= min_hi)) {
                fi(i,j,0,IceFR::S0mk) = salt_top;
                fi(i,j,0,IceFR::T0mk) = temp_top;
                fi(i,j,0,IceFR::W_ai) = zero;
                xtot = (one-si(i,j,0,cA))*fi(i,j,0,IceFR::W_ao);
            } else {
                // MK89 version
                fi(i,j,0,IceFR::W_io) = (Qio/RhoSW+Cp_w*cht*(fi(i,j,0,IceFR::T0mk)-temp_top))/hfus1;
                xtot = si(i,j,0,cA)*fi(i,j,0,IceFR::W_io)+(one-si(i,j,0,cA))*fi(i,j,0,IceFR::W_ao);
                fi(i,j,0,IceFR::S0mk) = (chs*salt_top+(xwai-fi(i,j,0,IceFR::W_io))*sice)/
                                        (chs+xwai+fi(i,j,0,IceFR::W_ro)-fi(i,j,0,IceFR::W_io));
                fi(i,j,0,IceFR::S0mk) = std::max(fi(i,j,0,IceFR::S0mk), zero);
                fi(i,j,0,IceFR::S0mk) = std::min(fi(i,j,0,IceFR::S0mk), Real(40.0));
                fi(i,j,0,IceFR::T0mk) = frln*fi(i,j,0,IceFR::S0mk);
            }

            // Ocean surface heat flux (overwritten, W/m2 then Celsius m/s) and salt flux (added)
            Real shf;
            if (si(i,j,0,cA) <= min_ai) {
                shf = Qnet_ao;
            } else {
                shf = (one-si(i,j,0,cA))*Qnet_ao+si(i,j,0,cA)*Qio-xtot*hfus1*RhoSW;
            }
            shf = -shf*RhoCpr;
            stflx(i,j,0,Temp_comp) = shf*rmask(i,j,0);

            stflx(i,j,0,Salt_comp) = stflx(i,j,0,Salt_comp)+
                                     ((si(i,j,0,cA)*(fi(i,j,0,IceFR::W_io)-fi(i,j,0,IceFR::W_ai))+
                                       (one-si(i,j,0,cA))*fi(i,j,0,IceFR::W_ao)+
                                       fi(i,j,0,IceFR::W_fr)))*(salt_top-sice)-
                                     si(i,j,0,cA)*(fi(i,j,0,IceFR::W_ro)-xwai)*salt_top;
            stflx(i,j,0,Salt_comp) *= rmask(i,j,0);

            // Ice-ocean mass flux
            fi(i,j,0,IceFR::IOmf) = xtot-si(i,j,0,cA)*xwai-si(i,j,0,cA)*fi(i,j,0,IceFR::W_ro)+
                                    fi(i,j,0,IceFR::W_fr);
            fi(i,j,0,IceFR::IOmf) *= rmask(i,j,0);

            // Update ice properties
            const Real mi_old = si(i,j,0,cH);
            const Real phi = (fi(i,j,0,IceFR::W_ao) < zero) ? Real(0.5) : Real(3.0);
            si(i,j,0,cH) = si(i,j,0,cH)+dtice*(si(i,j,0,cA)*(fi(i,j,0,IceFR::W_io)-fi(i,j,0,IceFR::W_ai))+
                                               (one-si(i,j,0,cA))*fi(i,j,0,IceFR::W_ao)+
                                               fi(i,j,0,IceFR::W_fr));
            const Real ai_tmp = si(i,j,0,cA);
            si(i,j,0,cA) = si(i,j,0,cA)+dtice*(one-si(i,j,0,cA))*
                                        (phi*fi(i,j,0,IceFR::W_ao)+fi(i,j,0,IceFR::W_fr));
            si(i,j,0,cA) = std::min(si(i,j,0,cA), max_ai);
            if (si(i,j,0,cA) < ai_tmp) {
                si(i,j,0,cS) = si(i,j,0,cS)*si(i,j,0,cA)/std::max(ai_tmp, eps);
            }

            if (convsnow) {
                // Raise a snow base below sea level by converting snow to ice
                const Real hstar = si(i,j,0,cS)-(RhoSW-IceRho)*si(i,j,0,cH)/SnowDryRho;
                if (hstar > zero) {
                    const Real cff = hstar/RhoSW;
                    si(i,j,0,cS) = si(i,j,0,cS)-IceRho*cff;
                    si(i,j,0,cH) = si(i,j,0,cH)+SnowDryRho*cff;
                }
            }

            // Ice age: new ice production lowers the parcel's age
            if ((si(i,j,0,cG) <= zero) && (si(i,j,0,cH) > min_hi)) {
                si(i,j,0,cG) = dtice*sec2day;
            } else if ((si(i,j,0,cG) > zero) && (si(i,j,0,cH) > min_hi)) {
                const Real delta_mi = std::min(std::max(si(i,j,0,cH)-mi_old, zero)/si(i,j,0,cH), one);
                si(i,j,0,cG) = si(i,j,0,cG)+dtice*sec2day-si(i,j,0,cG)*delta_mi;
            } else {
                si(i,j,0,cG) = zero;
            }

            si(i,j,0,cA) *= rmask(i,j,0);
            si(i,j,0,cH) *= rmask(i,j,0);

            // Limit the values
            si(i,j,0,cA) = std::min(si(i,j,0,cA), max_ai);
            si(i,j,0,cA) = std::max(si(i,j,0,cA), zero);
            si(i,j,0,cH) = std::max(si(i,j,0,cH), zero);
            si(i,j,0,cS) = std::max(si(i,j,0,cS), zero);
            si(i,j,0,cM) = std::max(si(i,j,0,cM), zero);
            si(i,j,0,cT) = std::max(si(i,j,0,cT), Real(-70.0));
            if (si(i,j,0,cH) <= zero) { si(i,j,0,cA) = zero; }
            if (si(i,j,0,cA) <= zero) { si(i,j,0,cH) = zero; }
        });
    }

    // The fields set above all lie between IOmf and W_ro. The state's halos and boundary
    // rows are left to ice_limit, which always follows and clamps cell by cell first.
    ice_exchange(lev, fir, IceFR::IOmf, IceFR::W_ro - IceFR::IOmf + 1);
    vec_stflx[lev]->FillBoundary(geom[lev].periodicity());

    ice_bc_r2d(lev, fir, IceFR::Isst);
    ice_bc_r2d(lev, fir, IceFR::Qcon);
    ice_bc_r2d(lev, fir, IceFR::Qrhs);
    ice_bc_r2d(lev, *vec_stflx[lev], Salt_comp);
    ice_bc_r2d(lev, *vec_stflx[lev], Temp_comp);
}

/**
 * ROMS ice_limit: positive concentration and thicknesses, bounded
 * temperature, and a cell with no ice concentration or thickness cleared,
 * on level linew; then the lateral BCs.
 */
void
REMORA::ice_limit (int lev)
{
    BL_PROFILE("REMORA::ice_limit()");

    const int linew = ice_lvl.linew;
    const int cA = rcomp(IceR::Aice, linew);
    const int cH = rcomp(IceR::Hice, linew);
    const int cS = rcomp(IceR::Hsno, linew);
    const int cM = rcomp(IceR::Hmel, linew);
    const int cT = rcomp(IceR::Tice, linew);
    const int cG = rcomp(IceR::Iage, linew);
    const Real max_ai = seaice_params.max_ai;

    MultiFab& sir = *vec_ice_si_r[lev];

    for (MFIter mfi(sir, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<Real> const& si = sir.array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            si(i,j,0,cA) = std::min(si(i,j,0,cA), max_ai);
            si(i,j,0,cA) = std::max(si(i,j,0,cA), zero);
            si(i,j,0,cH) = std::max(si(i,j,0,cH), zero);
            si(i,j,0,cS) = std::max(si(i,j,0,cS), zero);
            si(i,j,0,cM) = std::max(si(i,j,0,cM), zero);
            si(i,j,0,cT) = std::max(si(i,j,0,cT), Real(-70.0));
            si(i,j,0,cT) = std::min(si(i,j,0,cT), zero);
            si(i,j,0,cG) = std::max(si(i,j,0,cG), zero);
            if ((si(i,j,0,cA) <= zero) || (si(i,j,0,cH) <= zero)) {
                si(i,j,0,cA) = zero;
                si(i,j,0,cH) = zero;
                si(i,j,0,cM) = zero;
                si(i,j,0,cS) = zero;
                si(i,j,0,cG) = zero;
            }
        });
    }

    ice_exchange(lev, sir, rcomp(0, linew), IceR::NumVars);   // only linew changed
    ice_bc2d(lev, sir, cA);
    ice_bc2d(lev, sir, cH);
    ice_bc2d(lev, sir, cS);
    ice_bc2d(lev, sir, cM);
    ice_bc2d(lev, sir, cG);
    ice_tibc(lev);
}
