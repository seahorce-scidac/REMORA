#include <REMORA.H>

using namespace amrex;
using namespace REMORASeaIce;

// Ice momentum with elastic-viscous-plastic rheology (Hunke and Dukowicz
// 1997): ROMS ice_evp.F, ice_evp_sig.F and ice_elastic.F, run nEVP times per
// ice step from advance_seaice. Index offsets as in REMORA_SeaIce_thermo.cpp.

/**
 * ROMS ice_evp: on the first substep, load the ice velocity (at liuol, the
 * latest after the swap) into both EVP velocity levels; then strain rates
 * from the EVP velocity at lieol, ice strength, and bulk and shear viscosity.
 *
 * @param[in] lev  level to operate on
 * @param[in] iEVP substep, 1-based as in ROMS
 */
void
REMORA::ice_evp (int lev, int iEVP)
{
    BL_PROFILE("REMORA::ice_evp()");

    const int liold = ice_lvl.liold;
    const int liuol = ice_lvl.liuol;
    const int lieol = ice_lvl.lieol;
    const int cA = rcomp(IceR::Aice, liold);
    const int cH = rcomp(IceR::Hice, liold);
    const int cUe = ucomp(IceU::Uevp, lieol);
    const Real e2r = one/seaice_params.ellip_sq;
    const Real Pstar = seaice_params.Pstar;
    const Real Astrength = seaice_params.Astrength;
    const Real zetaMax = seaice_params.zetaMax;
    const Real zetaMin = seaice_params.zetaMin;
    constexpr Real epso = Real(1.0e-12);

    MultiFab& siu = *vec_ice_si_u[lev];
    MultiFab& siv = *vec_ice_si_v[lev];
    MultiFab& fir = *vec_ice_fi_r[lev];

    if (iEVP == 1) {
        // ROMS loads the interior, boundary and first ghost faces; the rest of the
        // ghost faces are never read
        for (int l = 0; l < 2; ++l) {
            MultiFab::Copy(siu, siu, ucomp(IceU::Uice, liuol), ucomp(IceU::Uevp, l), 1, siu.nGrowVect());
            MultiFab::Copy(siv, siv, ucomp(IceU::Uice, liuol), ucomp(IceU::Uevp, l), 1, siv.nGrowVect());
        }
    }

    for (MFIter mfi(fir, false); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        Array4<Real const> const& si = vec_ice_si_r[lev]->const_array(mfi);
        Array4<Real const> const& ue = siu.const_array(mfi, cUe);
        Array4<Real const> const& ve = siv.const_array(mfi, cUe);
        Array4<Real const> const& pm = vec_pm[lev]->const_array(mfi);
        Array4<Real const> const& pn = vec_pn[lev]->const_array(mfi);
        Array4<Real      > const& fi = fir.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real eps_xx = (ue(i+1,j,0)-ue(i,j,0))*pm(i,j,0);
            const Real eps_yy = (ve(i,j+1,0)-ve(i,j,0))*pn(i,j,0);
            const Real epx = Real(0.25)*(ve(i+1,j+1,0)+ve(i+1,j,0)-
                                         ve(i-1,j+1,0)-ve(i-1,j,0))*pm(i,j,0);
            const Real epy = Real(0.25)*(ue(i+1,j+1,0)+ue(i,j+1,0)-
                                         ue(i+1,j-1,0)-ue(i,j-1,0))*pn(i,j,0);
            const Real eps_xy = Real(0.5)*(epx+epy);

            const Real eone = eps_xx+eps_yy;
            const Real etwos = (eps_xx-eps_yy)*(eps_xx-eps_yy)+Real(4.0)*eps_xy*eps_xy;
            Real delta = std::abs(eone*eone+e2r*etwos);
            delta = std::max(std::sqrt(delta), epso);

            // ICE_STRENGTH_QUAD is not ported: ROMS never sets its gridded strength,
            // so it would give zero strength
            fi(i,j,0,IceFR::Pice) = Pstar*si(i,j,0,cH)*std::exp(-Astrength*(one-si(i,j,0,cA)));
            fi(i,j,0,IceFR::Bvis) = fi(i,j,0,IceFR::Pice)/(Real(2.0)*delta);
            fi(i,j,0,IceFR::Bvis) = std::min(fi(i,j,0,IceFR::Bvis), zetaMax);
            fi(i,j,0,IceFR::Bvis) = std::max(fi(i,j,0,IceFR::Bvis), zetaMin);
            fi(i,j,0,IceFR::Svis) = e2r*fi(i,j,0,IceFR::Bvis);
        });
    }

    // Only the three fields set above changed; they are not adjacent components
    ice_exchange(lev, fir, IceFR::Pice, 1);
    ice_exchange(lev, fir, IceFR::Bvis, 1);
    ice_exchange(lev, fir, IceFR::Svis, 1);
    ice_bc_r2d(lev, fir, IceFR::Pice);
    ice_bc_r2d(lev, fir, IceFR::Bvis);
    ice_bc_r2d(lev, fir, IceFR::Svis);
}

/**
 * ROMS ice_evp_sig: the internal ice stress at lienw from the stress and
 * strain rates at lieol, with an elastic term (Young modulus scaled by E0 =
 * 0.25) where the ice is thicker than 1 cm, and the viscous-plastic stress
 * where it is thinner.
 *
 * @param[in] lev   level to operate on
 * @param[in] dtevp EVP substep
 */
void
REMORA::ice_evp_sig (int lev, Real dtevp)
{
    BL_PROFILE("REMORA::ice_evp_sig()");

    const int liold = ice_lvl.liold;
    const int lieol = ice_lvl.lieol;
    const int lienw = ice_lvl.lienw;
    const int cH = rcomp(IceR::Hice, liold);
    const int cUe = ucomp(IceU::Uevp, lieol);
    const int cXXo = rcomp(IceR::ISxx, lieol), cXXn = rcomp(IceR::ISxx, lienw);
    const int cXYo = rcomp(IceR::ISxy, lieol), cXYn = rcomp(IceR::ISxy, lienw);
    const int cYYo = rcomp(IceR::ISyy, lieol), cYYn = rcomp(IceR::ISyy, lienw);
    const Real IceRho = seaice_params.IceRho;
    constexpr Real E0 = Real(0.25);

    MultiFab& sir = *vec_ice_si_r[lev];

    for (MFIter mfi(sir, false); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        Array4<Real      > const& si = sir.array(mfi);
        Array4<Real const> const& fi = vec_ice_fi_r[lev]->const_array(mfi);
        Array4<Real const> const& ue = vec_ice_si_u[lev]->const_array(mfi, cUe);
        Array4<Real const> const& ve = vec_ice_si_v[lev]->const_array(mfi, cUe);
        Array4<Real const> const& pm = vec_pm[lev]->const_array(mfi);
        Array4<Real const> const& pn = vec_pn[lev]->const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real eps_xx = (ue(i+1,j,0)-ue(i,j,0))*pm(i,j,0);
            const Real eps_yy = (ve(i,j+1,0)-ve(i,j,0))*pn(i,j,0);
            const Real epx = Real(0.25)*(ve(i+1,j+1,0)+ve(i+1,j,0)-
                                         ve(i-1,j+1,0)-ve(i-1,j,0))*pm(i,j,0);
            const Real epy = Real(0.25)*(ue(i+1,j+1,0)+ue(i,j+1,0)-
                                         ue(i+1,j-1,0)-ue(i,j-1,0))*pn(i,j,0);
            const Real eps_xy = Real(0.5)*(epx+epy);

            const Real Pice = fi(i,j,0,IceFR::Pice);
            const Real Bvis = fi(i,j,0,IceFR::Bvis);
            const Real Svis = fi(i,j,0,IceFR::Svis);

            if (si(i,j,0,cH) > Real(0.01)) {
                const Real pmdt = pm(i,j,0)*dtevp;
                const Real E   = Real(2.0)*E0*IceRho*si(i,j,0,cH)/(pmdt*pmdt);
                const Real ep  = E*Pice/(Real(4.0)*Bvis+Real(1.0e-8));
                const Real ee  = E/(Real(2.0)*Svis+Real(1.0e-8));
                const Real ees = E*(Svis-Bvis)/(Real(4.0)*Svis*Bvis+Real(1.0e-8));

                const Real cff   = one/dtevp;
                const Real alfa  = cff+ee+ees;
                const Real beta  = ees;
                const Real gamma = one/dtevp+ee;
                const Real f1 = E*eps_xx+cff*si(i,j,0,cXXo)-ep;
                const Real f2 = E*eps_yy+cff*si(i,j,0,cYYo)-ep;
                const Real f3 = E*eps_xy+cff*si(i,j,0,cXYo);
                const Real f4 = one/(alfa*alfa-beta*beta);
                si(i,j,0,cXXn) = f4*(alfa*f1-beta*f2);
                si(i,j,0,cXYn) = f3/gamma;
                si(i,j,0,cYYn) = f4*(alfa*f2-beta*f1);
            } else {
                si(i,j,0,cXXn) = Real(2.0)*Svis*eps_xx+(Bvis-Svis)*(eps_xx+eps_yy)-Real(0.5)*Pice;
                si(i,j,0,cYYn) = Real(2.0)*Svis*eps_yy+(Bvis-Svis)*(eps_xx+eps_yy)-Real(0.5)*Pice;
                si(i,j,0,cXYn) = Real(2.0)*Svis*eps_xy;
            }
        });
    }

    // Only the new stresses changed: ISxx, ISxy and ISyy are adjacent
    static_assert(IceR::ISxy == IceR::ISxx + 1 && IceR::ISyy == IceR::ISxx + 2);
    ice_exchange(lev, sir, cXXn, 3);
    ice_bc2d(lev, sir, cXXn);
    ice_bc2d(lev, sir, cXYn);
    ice_bc2d(lev, sir, cYYn);
}

/**
 * ROMS ice_elastic: one EVP substep of the ice momentum equations, implicit in
 * the ice-water drag. Forces: divergence of the internal stress (the
 * cross-derivative of Sxy falls back to a two-point average next to land),
 * the air-ice stress, the sea surface tilt, Coriolis on the mask-weighted
 * average of the other velocity component, and the drag of the mixed-layer
 * ocean velocity. No turning angle is applied, as in ROMS. On the last
 * substep the EVP velocity becomes the ice velocity at liunw.
 *
 * @param[in] lev   level to operate on
 * @param[in] iEVP  substep, 1-based
 * @param[in] dtevp EVP substep
 */
void
REMORA::ice_elastic (int lev, int iEVP, Real dtevp)
{
    BL_PROFILE("REMORA::ice_elastic()");

    const int liold = ice_lvl.liold;
    const int liunw = ice_lvl.liunw;
    const int lieol = ice_lvl.lieol;
    const int lienw = ice_lvl.lienw;
    const int cA = rcomp(IceR::Aice, liold);
    const int cH = rcomp(IceR::Hice, liold);
    const int cXX = rcomp(IceR::ISxx, lienw);
    const int cXY = rcomp(IceR::ISxy, lienw);
    const int cYY = rcomp(IceR::ISyy, lienw);
    const int cUeo = ucomp(IceU::Uevp, lieol);
    const int cUen = ucomp(IceU::Uevp, lienw);
    const Real IceRho = seaice_params.IceRho;
    const Real l_g = solverChoice.g;
    const Real rho0 = solverChoice.rho0;

    MultiFab& siu = *vec_ice_si_u[lev];
    MultiFab& siv = *vec_ice_si_v[lev];

    for (MFIter mfi(siu, false); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();   // cell-centered valid box of this grid
        const Box cbx = amrex::enclosedCells(bx, 0);
        Array4<Real const> const& si   = vec_ice_si_r[lev]->const_array(mfi);
        Array4<Real const> const& fi   = vec_ice_fi_r[lev]->const_array(mfi);
        Array4<Real const> const& fu   = vec_ice_fi_u[lev]->const_array(mfi);
        Array4<Real const> const& fv   = vec_ice_fi_v[lev]->const_array(mfi);
        Array4<Real const> const& ueo  = siu.const_array(mfi, cUeo);
        Array4<Real const> const& veo  = siv.const_array(mfi, cUeo);
        Array4<Real      > const& uen  = siu.array(mfi, cUen);
        Array4<Real      > const& ven  = siv.array(mfi, cUen);
        Array4<Real const> const& pm   = vec_pm[lev]->const_array(mfi);
        Array4<Real const> const& pn   = vec_pn[lev]->const_array(mfi);
        Array4<Real const> const& f    = vec_fcor[lev]->const_array(mfi);
        Array4<Real const> const& umask = vec_msku[lev]->const_array(mfi);
        Array4<Real const> const& vmask = vec_mskv[lev]->const_array(mfi);

        // U-momentum
        ParallelFor(interior_faces(amrex::surroundingNodes(cbx,0), 0, geom[lev]),
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real om_u = Real(2.0)/(pm(i-1,j,0)+pm(i,j,0));
            const Real on_u = Real(2.0)/(pn(i-1,j,0)+pn(i,j,0));
            const Real area = om_u*on_u;
            const Real chux = Real(0.5)*(fi(i-1,j,0,IceFR::IOmt)+fi(i,j,0,IceFR::IOmt));
            Real masu = Real(0.5)*(si(i,j,0,cH)+si(i-1,j,0,cH));
            masu = std::max(masu, Real(0.1));
            masu = masu*IceRho;
            const Real auf = std::max(Real(0.5)*(si(i-1,j,0,cA)+si(i,j,0,cA)), Real(0.01));
            const Real pmu = Real(0.5)*(pm(i,j,0)+pm(i-1,j,0));
            const Real pnu = Real(0.5)*(pn(i,j,0)+pn(i-1,j,0));

            // Forces from the ice rheology
            const Real s1 = (si(i,j,0,cXX)-si(i-1,j,0,cXX))*pmu;
            Real f1, f2;
            if ((umask(i,j,0) > zero) && (umask(i,j+1,0) < one)) {
                f1 = Real(0.5)*(si(i,j,0,cXY)+si(i-1,j,0,cXY));
            } else {
                f1 = Real(0.25)*(si(i,j,0,cXY)+si(i,j+1,0,cXY)+si(i-1,j+1,0,cXY)+si(i-1,j,0,cXY));
            }
            if ((umask(i,j,0) > zero) && (umask(i,j-1,0) < one)) {
                f2 = Real(0.5)*(si(i,j,0,cXY)+si(i-1,j,0,cXY));
            } else {
                f2 = Real(0.25)*(si(i,j,0,cXY)+si(i-1,j,0,cXY)+si(i-1,j-1,0,cXY)+si(i,j-1,0,cXY));
            }
            const Real s2 = (f1-f2)*pnu;
            Real uforce = (s1+s2)*area;

            // Wind forcing
            uforce += IceRho*fu(i,j,0,IceFU::AIus)*area;

            // Tilt of the ocean surface
            uforce -= l_g*masu*on_u*(fi(i,j,0,IceFR::Hsse)-fi(i-1,j,0,IceFR::Hsse));

            // Coriolis
            Real fakt = zero;
            const Real dsum = vmask(i-1,j,0)+vmask(i,j,0)+vmask(i-1,j+1,0)+vmask(i,j+1,0);
            if (dsum > zero) { fakt = one/dsum; }
            const Real mfv11 = Real(0.5)*(si(i-1,j-1,0,cH)*f(i-1,j-1,0)+si(i-1,j,0,cH)*f(i-1,j,0))*
                               veo(i-1,j,0)*vmask(i-1,j,0);
            const Real mfv21 = Real(0.5)*(si(i,j-1,0,cH)*f(i,j-1,0)+si(i,j,0,cH)*f(i,j,0))*
                               veo(i,j,0)*vmask(i,j,0);
            const Real mfv12 = Real(0.5)*(si(i-1,j,0,cH)*f(i-1,j,0)+si(i-1,j+1,0,cH)*f(i-1,j+1,0))*
                               veo(i-1,j+1,0)*vmask(i-1,j+1,0);
            const Real mfv22 = Real(0.5)*(si(i,j,0,cH)*f(i,j,0)+si(i,j+1,0,cH)*f(i,j+1,0))*
                               veo(i,j+1,0)*vmask(i,j+1,0);
            uforce += fakt*IceRho*area*(mfv11+mfv21+mfv12+mfv22);

            // Stress from the ocean current
            uforce = uforce/area+auf*chux*rho0*fu(i,j,0,IceFU::Uavg);

            const Real alfa = masu+dtevp*auf*rho0*chux;
            uen(i,j,0) = (masu*ueo(i,j,0)+dtevp*uforce)/alfa;
            uen(i,j,0) *= umask(i,j,0);
        });

        // V-momentum
        ParallelFor(interior_faces(amrex::surroundingNodes(cbx,1), 1, geom[lev]),
        [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real om_v = Real(2.0)/(pm(i,j-1,0)+pm(i,j,0));
            const Real on_v = Real(2.0)/(pn(i,j-1,0)+pn(i,j,0));
            const Real area = om_v*on_v;
            Real masv = Real(0.5)*(si(i,j,0,cH)+si(i,j-1,0,cH));
            masv = std::max(masv, Real(0.1));
            masv = masv*IceRho;
            const Real avf = std::max(Real(0.5)*(si(i,j-1,0,cA)+si(i,j,0,cA)), Real(0.01));
            const Real chuy = Real(0.5)*(fi(i,j-1,0,IceFR::IOmt)+fi(i,j,0,IceFR::IOmt));
            const Real pmv = Real(0.5)*(pm(i,j,0)+pm(i,j-1,0));
            const Real pnv = Real(0.5)*(pn(i,j,0)+pn(i,j-1,0));

            const Real s1 = (si(i,j,0,cYY)-si(i,j-1,0,cYY))*pnv;
            Real f1, f2;
            if ((vmask(i,j,0) > zero) && (vmask(i+1,j,0) < one)) {
                f1 = Real(0.5)*(si(i,j,0,cXY)+si(i,j-1,0,cXY));
            } else {
                f1 = Real(0.25)*(si(i,j,0,cXY)+si(i+1,j,0,cXY)+si(i+1,j-1,0,cXY)+si(i,j-1,0,cXY));
            }
            if ((vmask(i,j,0) > zero) && (vmask(i-1,j,0) < one)) {
                f2 = Real(0.5)*(si(i,j,0,cXY)+si(i,j-1,0,cXY));
            } else {
                f2 = Real(0.25)*(si(i,j,0,cXY)+si(i-1,j,0,cXY)+si(i-1,j-1,0,cXY)+si(i,j-1,0,cXY));
            }
            const Real s2 = (f1-f2)*pmv;
            Real vforce = (s1+s2)*area;

            vforce += IceRho*fv(i,j,0,IceFU::AIus)*area;

            vforce -= l_g*masv*om_v*(fi(i,j,0,IceFR::Hsse)-fi(i,j-1,0,IceFR::Hsse));

            Real fakt = zero;
            const Real dsum = umask(i,j-1,0)+umask(i+1,j-1,0)+umask(i,j,0)+umask(i+1,j,0);
            if (dsum > zero) { fakt = one/dsum; }
            const Real mfu11 = Real(0.5)*(si(i-1,j-1,0,cH)*f(i-1,j-1,0)+si(i,j-1,0,cH)*f(i,j-1,0))*
                               ueo(i,j-1,0)*umask(i,j-1,0);
            const Real mfu21 = Real(0.5)*(si(i,j-1,0,cH)*f(i,j-1,0)+si(i+1,j-1,0,cH)*f(i+1,j-1,0))*
                               ueo(i+1,j-1,0)*umask(i+1,j-1,0);
            const Real mfu12 = Real(0.5)*(si(i-1,j,0,cH)*f(i-1,j,0)+si(i,j,0,cH)*f(i,j,0))*
                               ueo(i,j,0)*umask(i,j,0);
            const Real mfu22 = Real(0.5)*(si(i,j,0,cH)*f(i,j,0)+si(i+1,j,0,cH)*f(i+1,j,0))*
                               ueo(i+1,j,0)*umask(i+1,j,0);
            vforce -= fakt*IceRho*area*(mfu11+mfu21+mfu12+mfu22);

            vforce = vforce/area+avf*chuy*rho0*fv(i,j,0,IceFU::Uavg);

            const Real alfa = masv+dtevp*avf*rho0*chuy;
            ven(i,j,0) = (masv*veo(i,j,0)+dtevp*vforce)/alfa;
            ven(i,j,0) *= vmask(i,j,0);
        });
    }

    ice_exchange(lev, siu, cUen, 1);
    ice_exchange(lev, siv, cUen, 1);
    ice_uibc(lev, siu, cUen);
    ice_vibc(lev, siv, cUen);

    if (iEVP == seaice_params.nEVP) {
        // ROMS copies the interior, boundary and first ghost faces
        MultiFab::Copy(siu, siu, cUen, ucomp(IceU::Uice, liunw), 1, IntVect(1,1,0));
        MultiFab::Copy(siv, siv, cUen, ucomp(IceU::Uice, liunw), 1, IntVect(1,1,0));
        ice_exchange(lev, siu, ucomp(IceU::Uice, liunw), 1);
        ice_exchange(lev, siv, ucomp(IceU::Uice, liunw), 1);
    }
}
