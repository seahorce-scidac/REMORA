#include <REMORA.H>

using namespace amrex;
using namespace REMORASeaIce;

// Ice advection: ROMS ice_enthalpy.F and ice_smolar.h (Smolarkiewicz MPDATA,
// or first-order upwind with ICE_UPWIND). Fields are advected from liold to
// linew with the ice velocity at liunw. Index offsets as in
// REMORA_SeaIce_thermo.cpp.

/**
 * ROMS ice_enthalpy: enthalpy as thickness times interior temperature, on
 * both levels, zero on the old level where the old thickness is at most
 * min_hi.
 */
void
REMORA::ice_enthalpy (int lev)
{
    const int liold = ice_lvl.liold;
    const int linew = ice_lvl.linew;
    const int cHo = rcomp(IceR::Hice, liold), cHn = rcomp(IceR::Hice, linew);
    const int cTo = rcomp(IceR::Tice, liold), cTn = rcomp(IceR::Tice, linew);
    const int cEo = rcomp(IceR::Enth, liold), cEn = rcomp(IceR::Enth, linew);
    const Real min_hi = seaice_params.min_hi;

    // Pointwise over the ghost cells too, from thickness and temperature whose halos are
    // current, so the enthalpy needs no exchange of its own
    MultiFab& sir = *vec_ice_si_r[lev];
    for (MFIter mfi(sir, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<Real> const& si = sir.array(mfi);
        ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            si(i,j,0,cEo) = si(i,j,0,cHo)*si(i,j,0,cTo);
            si(i,j,0,cEn) = si(i,j,0,cHn)*si(i,j,0,cTn);
            if (si(i,j,0,cHo) <= min_hi) { si(i,j,0,cEo) = zero; }
        });
    }
}

/**
 * ROMS ice_mpdata_tile: advect one rho-point ice field from component cold to
 * cnew. Upstream fluxes give a first guess; unless ice.advect_upwind, one
 * antidiffusive MPDATA correction follows, with no flux limiter (ROMS's
 * nonoscillatory shift is declared but unused). Only the interior is written;
 * the caller exchanges and applies the BCs.
 */
void
REMORA::ice_mpdata (int lev, int cold, int cnew, Real dt_lev)
{
    const int liunw = ice_lvl.liunw;
    const int cU = ucomp(IceU::Uice, liunw);
    const bool upwind = seaice_params.advect_upwind;
    const Real dtice = dt_lev;
    constexpr Real epsil = Real(1.0e-15);

    const Box& dom = geom[lev].Domain();
    const bool xper = geom[lev].isPeriodic(0);
    const bool yper = geom[lev].isPeriodic(1);

    MultiFab& sir = *vec_ice_si_r[lev];

    for (MFIter mfi(sir, false); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        Array4<Real> const& fld = sir.array(mfi);
        Array4<Real const> const& ui = vec_ice_si_u[lev]->const_array(mfi, cU);
        Array4<Real const> const& vi = vec_ice_si_v[lev]->const_array(mfi, cU);
        Array4<Real const> const& pm = vec_pm[lev]->const_array(mfi);
        Array4<Real const> const& pn = vec_pn[lev]->const_array(mfi);
        Array4<Real const> const& rmask = vec_mskr[lev]->const_array(mfi);
        Array4<Real const> const& umask = vec_msku[lev]->const_array(mfi);
        Array4<Real const> const& vmask = vec_mskv[lev]->const_array(mfi);

        // ROMS Imin..Imax: the interior for upwind; one cell more for MPDATA,
        // clipped to the domain interior (ROMS 1..Lm) unless periodic
        Box abx = bx;
        if (!upwind) {
            abx.grow(IntVect(1,1,0));
            if (!xper) {
                abx.setSmall(0, std::max(abx.smallEnd(0), dom.smallEnd(0)));
                abx.setBig  (0, std::min(abx.bigEnd(0),   dom.bigEnd(0)));
            }
            if (!yper) {
                abx.setSmall(1, std::max(abx.smallEnd(1), dom.smallEnd(1)));
                abx.setBig  (1, std::min(abx.bigEnd(1),   dom.bigEnd(1)));
            }
        }

        // Scratch over one cell around the box: aif, and the face fluxes
        const Box sbx = amrex::grow(bx, IntVect(1,1,0));
        FArrayBox aif_fab(sbx, 1, The_Async_Arena());
        FArrayBox fx_fab(amrex::surroundingNodes(sbx,0), 2, The_Async_Arena());
        FArrayBox fy_fab(amrex::surroundingNodes(sbx,1), 2, The_Async_Arena());
        Array4<Real> const& aif   = aif_fab.array();
        Array4<Real> const& aflxu = fx_fab.array(0);
        Array4<Real> const& FX    = fx_fab.array(1);
        Array4<Real> const& aflxv = fy_fab.array(0);
        Array4<Real> const& FE    = fy_fab.array(1);
        aif_fab.template setVal<RunOn::Device>(zero);

        // Upstream fluxes on the faces of abx
        Box ufaces = amrex::surroundingNodes(abx, 0);
        Box vfaces = amrex::surroundingNodes(abx, 1);
        ParallelFor(ufaces, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real on_u = Real(2.0)/(pn(i-1,j,0)+pn(i,j,0));
            const Real cff1 = std::max(zero, ui(i,j,0));
            const Real cff2 = std::min(zero, ui(i,j,0));
            aflxu(i,j,0) = on_u*(cff1*fld(i-1,j,0,cold)+cff2*fld(i,j,0,cold));
        });
        ParallelFor(vfaces, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real om_v = Real(2.0)/(pm(i,j-1,0)+pm(i,j,0));
            const Real cff1 = std::max(zero, vi(i,j,0));
            const Real cff2 = std::min(zero, vi(i,j,0));
            aflxv(i,j,0) = om_v*(cff1*fld(i,j-1,0,cold)+cff2*fld(i,j,0,cold));
        });

        // Step 1: first guess
        ParallelFor(abx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real omn = one/(pm(i,j,0)*pn(i,j,0));   // ROMS metrics
            const Real ar = one/omn;
            Real a = fld(i,j,0,cold)-dtice*(aflxu(i+1,j,0)-aflxu(i,j,0)+
                                            aflxv(i,j+1,0)-aflxv(i,j,0))*ar;
            aif(i,j,0) = a*rmask(i,j,0);
        });

        // Values on the boundary rows: zero gradient
        const bool west  = !xper && bx.smallEnd(0) == dom.smallEnd(0);
        const bool east  = !xper && bx.bigEnd(0)   == dom.bigEnd(0);
        const bool south = !yper && bx.smallEnd(1) == dom.smallEnd(1);
        const bool north = !yper && bx.bigEnd(1)   == dom.bigEnd(1);
        const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0);
        const int jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);
        const int jmn = abx.smallEnd(1), jmx = abx.bigEnd(1);
        const int imn = abx.smallEnd(0), imx = abx.bigEnd(0);
        if (west) {
            ParallelFor(Box(IntVect(ilo-1,jmn,0), IntVect(ilo-1,jmx,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                aif(i,j,0) = aif(i+1,j,0);
            });
        }
        if (east) {
            ParallelFor(Box(IntVect(ihi+1,jmn,0), IntVect(ihi+1,jmx,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                aif(i,j,0) = aif(i-1,j,0);
            });
        }
        if (south) {
            ParallelFor(Box(IntVect(imn,jlo-1,0), IntVect(imx,jlo-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                aif(i,j,0) = aif(i,j+1,0);
            });
        }
        if (north) {
            ParallelFor(Box(IntVect(imn,jhi+1,0), IntVect(imx,jhi+1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                aif(i,j,0) = aif(i,j-1,0);
            });
        }
        if (south && !xper) {
            if (west) {
                ParallelFor(Box(IntVect(ilo-1,jlo-1,0), IntVect(ilo-1,jlo-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    aif(i,j,0) = aif(i+1,j+1,0);
                });
            }
            if (east) {
                ParallelFor(Box(IntVect(ihi+1,jlo-1,0), IntVect(ihi+1,jlo-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    aif(i,j,0) = aif(i-1,j+1,0);
                });
            }
        }
        if (north && !xper) {
            if (west) {
                ParallelFor(Box(IntVect(ilo-1,jhi+1,0), IntVect(ilo-1,jhi+1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    aif(i,j,0) = aif(i+1,j-1,0);
                });
            }
            if (east) {
                ParallelFor(Box(IntVect(ihi+1,jhi+1,0), IntVect(ihi+1,jhi+1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    aif(i,j,0) = aif(i-1,j-1,0);
                });
            }
        }

        if (!upwind) {
            // Step 2: antidiffusive velocities
            ParallelFor(Box(IntVect(ilo-1,jlo,0), IntVect(ihi+1,jhi+1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                FE(i,j,0) = Real(0.5)*vmask(i,j,0)*(aif(i,j,0)-aif(i,j-1,0));
            });
            ParallelFor(Box(IntVect(ilo,jlo-1,0), IntVect(ihi+1,jhi+1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                FX(i,j,0) = Real(0.5)*umask(i,j,0)*(aif(i,j,0)-aif(i-1,j,0));
            });

            ParallelFor(Box(IntVect(ilo,jlo,0), IntVect(ihi+1,jhi,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                const Real rateu = (aif(i,j,0)-aif(i-1,j,0))/std::max(epsil, aif(i,j,0)+aif(i-1,j,0));
                const Real rateyiu = (FE(i,j+1,0)+FE(i,j,0)+FE(i-1,j+1,0)+FE(i-1,j,0))/
                                     (std::max(epsil, aif(i,j,0)+FE(i,j+1,0)-FE(i,j,0)+
                                                      aif(i-1,j,0)+FE(i-1,j+1,0)-FE(i-1,j,0)));
                const Real Cu = Real(0.5)*dtice*(pm(i,j,0)+pm(i-1,j,0))*ui(i,j,0);
                const Real Cu_crss = Real(0.5)*dtice*Real(0.0625)*
                                     (pn(i-1,j+1,0)+pn(i,j+1,0)+pn(i-1,j-1,0)+pn(i,j-1,0))*
                                     (vi(i-1,j+1,0)+vi(i,j+1,0)+vi(i-1,j,0)+vi(i,j,0));
                const Real uspeed = rateu*(std::abs(ui(i,j,0))-Cu*ui(i,j,0))-rateyiu*Cu_crss*ui(i,j,0);
                const Real cff1 = std::max(zero, uspeed);
                const Real cff2 = std::min(zero, uspeed);
                const Real on_u = Real(2.0)/(pn(i-1,j,0)+pn(i,j,0));
                aflxu(i,j,0) = on_u*(cff1*aif(i-1,j,0)+cff2*aif(i,j,0));
            });
            ParallelFor(Box(IntVect(ilo,jlo,0), IntVect(ihi,jhi+1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                const Real ratev = (aif(i,j,0)-aif(i,j-1,0))/std::max(epsil, aif(i,j,0)+aif(i,j-1,0));
                const Real ratexiv = (FX(i+1,j,0)+FX(i,j,0)+FX(i+1,j-1,0)+FX(i,j-1,0))/
                                     (std::max(epsil, aif(i,j,0)+FX(i+1,j,0)-FX(i,j,0)+
                                                      aif(i,j-1,0)+FX(i+1,j-1,0)-FX(i,j-1,0)));
                const Real Cu = Real(0.5)*dtice*(pn(i,j,0)+pn(i,j-1,0))*vi(i,j,0);
                const Real Cu_crss = Real(0.5)*dtice*Real(0.0625)*
                                     (pm(i+1,j,0)+pm(i+1,j-1,0)+pm(i-1,j,0)+pm(i-1,j-1,0))*
                                     (ui(i,j,0)+ui(i+1,j,0)+ui(i,j-1,0)+ui(i+1,j-1,0));
                const Real vspeed = ratev*(std::abs(vi(i,j,0))-Cu*vi(i,j,0))-ratexiv*Cu_crss*vi(i,j,0);
                const Real cff1 = std::max(zero, vspeed);
                const Real cff2 = std::min(zero, vspeed);
                const Real om_v = Real(2.0)/(pm(i,j-1,0)+pm(i,j,0));
                aflxv(i,j,0) = om_v*(cff1*aif(i,j-1,0)+cff2*aif(i,j,0));
            });

            // Correct by the antidiffusive fluxes, then load
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                Real a = aif(i,j,0)-dtice*pm(i,j,0)*pn(i,j,0)*
                                    (aflxu(i+1,j,0)-aflxu(i,j,0)+aflxv(i,j+1,0)-aflxv(i,j,0));
                fld(i,j,0,cnew) = a*rmask(i,j,0);
            });
        } else {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                fld(i,j,0,cnew) = aif(i,j,0);
            });
        }
    }

    ice_exchange(lev, sir, cnew, 1);
}

/**
 * ROMS ice_advect (ice_smolar.h): advect concentration and thickness, then,
 * with thermodynamics on, snow, melt pond, enthalpy and age-weighted
 * thickness, recovering the interior temperature and age from them.
 *
 * ROMS negates the whole enthalpy array inside each tile's call, so that
 * MPDATA sees a positive field; under OpenMP the tiles race on it. Here it is
 * negated once, for the whole level, before and after.
 */
void
REMORA::ice_advect (int lev, Real dt_lev)
{
    BL_PROFILE("REMORA::ice_advect()");

    const int liold = ice_lvl.liold;
    const int linew = ice_lvl.linew;
    const Real min_ai = seaice_params.min_ai;
    const Real min_hi = seaice_params.min_hi;
    MultiFab& sir = *vec_ice_si_r[lev];
    MultiFab& fir = *vec_ice_fi_r[lev];

    auto advect = [&] (int var) {
        ice_mpdata(lev, rcomp(var, liold), rcomp(var, linew), dt_lev);
    };

    advect(IceR::Aice);
    ice_bc2d(lev, sir, rcomp(IceR::Aice, linew));

    advect(IceR::Hice);
    ice_bc2d(lev, sir, rcomp(IceR::Hice, linew));

    // Rate of ice divergence
    {
        const int cHo = rcomp(IceR::Hice, liold), cHn = rcomp(IceR::Hice, linew);
        for (MFIter mfi(sir, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Array4<Real const> const& si = sir.const_array(mfi);
            Array4<Real      > const& fi = fir.array(mfi);
            ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                fi(i,j,0,IceFR::Wdiv) = (si(i,j,0,cHn)-si(i,j,0,cHo))/dt_lev;
            });
        }
    }

    if (seaice_params.do_thermo) {
        advect(IceR::Hsno);
        ice_bc2d(lev, sir, rcomp(IceR::Hsno, linew));

        advect(IceR::Hmel);
        ice_bc2d(lev, sir, rcomp(IceR::Hmel, linew));

        // Interior temperature, through the (negated) enthalpy. The ghost cells are
        // in the horizontal only, so negate over growntilebox rather than with
        // MultiFab::mult, whose scalar ghost count would also grow the 2D boxes in z.
        auto negate_enthalpy = [&] () {
            const int cE0 = rcomp(IceR::Enth, 0), cE1 = rcomp(IceR::Enth, 1);
            for (MFIter mfi(sir, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Array4<Real> const& si = sir.array(mfi);
                ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    si(i,j,0,cE0) = -si(i,j,0,cE0);
                    si(i,j,0,cE1) = -si(i,j,0,cE1);
                });
            }
        };
        negate_enthalpy();
        advect(IceR::Enth);
        negate_enthalpy();
        {
            const int cA = rcomp(IceR::Aice, linew);
            const int cH = rcomp(IceR::Hice, linew);
            const int cT = rcomp(IceR::Tice, linew);
            const int cE = rcomp(IceR::Enth, linew);
            for (MFIter mfi(sir, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Array4<Real> const& si = sir.array(mfi);
                ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    si(i,j,0,cT) = si(i,j,0,cE)/std::max(si(i,j,0,cH), Real(1.0e-6));
                    if (si(i,j,0,cA) <= min_ai) {
                        si(i,j,0,cE) = zero;
                        si(i,j,0,cT) = Real(-2.0);
                    }
                });
            }
        }
        ice_tibc(lev);

        // Age, through the age-weighted thickness
        {
            const int cHo = rcomp(IceR::Hice, liold), cHn = rcomp(IceR::Hice, linew);
            const int cGo = rcomp(IceR::Iage, liold), cGn = rcomp(IceR::Iage, linew);
            const int cKo = rcomp(IceR::Hage, liold), cKn = rcomp(IceR::Hage, linew);
            for (MFIter mfi(sir, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Array4<Real> const& si = sir.array(mfi);
                ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    si(i,j,0,cKo) = si(i,j,0,cHo)*si(i,j,0,cGo);
                    si(i,j,0,cKn) = si(i,j,0,cHn)*si(i,j,0,cGn);
                    if (si(i,j,0,cHo) <= min_hi) {
                        si(i,j,0,cKo) = zero;
                        si(i,j,0,cGo) = zero;
                    }
                });
            }
            advect(IceR::Hage);
            for (MFIter mfi(sir, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Array4<Real> const& si = sir.array(mfi);
                ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    si(i,j,0,cGn) = si(i,j,0,cKn)/std::max(si(i,j,0,cHn), Real(1.0e-6));
                    if (si(i,j,0,cHn) <= min_hi) {
                        si(i,j,0,cKn) = zero;
                        si(i,j,0,cGn) = zero;
                    }
                });
            }
            ice_bc2d(lev, sir, cGn);
        }
    }

    // The ice age and age-thickness changed at liold too. The linew halos are left to
    // ice_limit, which always follows and clamps cell by cell first.
    ice_exchange(lev, sir, rcomp(IceR::Iage, liold), 1);
    ice_exchange(lev, sir, rcomp(IceR::Hage, liold), 1);
}
