#include <REMORA.H>

using namespace amrex;
using namespace REMORASeaIce;

/**
 * Ice/snow albedo: ROMS ice_albedo.F. Called before the bulk flux, as in
 * main3d. ROMS computes a seawater albedo here too, but its bulk_flux has the
 * line applying it commented out, so it is left out.
 *
 * The default branch is kept with its quirk: it compares the ice surface
 * temperature, which is in Celsius, as if it were in Kelvin, so it always
 * takes the dry albedo.
 *
 * @param[in] lev level to operate on
 */
void
REMORA::ice_albedo (int lev)
{
    BL_PROFILE("REMORA::ice_albedo()");

    const IceAlbedoType albedo_type = seaice_params.albedo_type;
    MultiFab& fi = *vec_ice_fi_r[lev];

    if (albedo_type == IceAlbedoType::none) {
        fi.setVal(zero, IceFR::albedo_ice, 1, fi.nGrowVect());
        return;
    }

    // No PerfectRST restart here, so ROMS's li_stp is always linew.
    const int li_stp = ice_lvl.linew;
    const int cA = rcomp(IceR::Aice, li_stp);
    const int cH = rcomp(IceR::Hice, li_stp);
    const int cS = rcomp(IceR::Hsno, li_stp);
    const int cM = rcomp(IceR::Hmel, li_stp);
    const Real min_ai = seaice_params.min_ai;

    constexpr Real alb_w = Real(0.06);
    constexpr Real alb_i_dry = Real(0.65);
    constexpr Real alb_i_wet = Real(0.60);
    constexpr Real alb_s_dry = Real(0.85);
    constexpr Real alb_s_wet = Real(0.72);

    for (MFIter mfi(fi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<Real const> const& si = vec_ice_si_r[lev]->const_array(mfi);
        Array4<Real      > const& f  = fi.array(mfi);
        const Box gbx = makeSlab(mfi.growntilebox(IntVect(NGROW,NGROW,0)), 2, 0);

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            const Real ai = si(i,j,0,cA);
            const Real cff = one / (ai + Real(0.001));
            const Real ice_thick  = cff * si(i,j,0,cH);
            const Real snow_thick = cff * si(i,j,0,cS);
            const Real sfc_temp = f(i,j,0,IceFR::Isst);
            Real albedo_ice;

            if (albedo_type == IceAlbedoType::ec92) {
                // Ebert and Curry (1992)
                if (ai > min_ai) {
                    const Real Hice = std::max(ice_thick, Real(0.01));
                    const Real Hsnow = snow_thick;
                    Real alb_ice;
                    if (Hice >= Real(2.0)) {
                        alb_ice = Real(0.561632);
                    } else if (Hice >= one) {
                        alb_ice = Real(0.07616) * Hice + Real(0.414492);
                    } else {
                        alb_ice = Real(0.082409) * std::log(Hice) + Real(0.485472);
                    }
                    Real alb_snow;
                    if (si(i,j,0,cM) > zero) {
                        if (Hsnow >= Real(0.1)) {
                            alb_snow = Real(0.701009);
                        } else {
                            alb_snow = alb_ice + (Real(0.701009) - alb_ice) * Hsnow / Real(0.1);
                        }
                    } else {
                        alb_snow = Real(0.83);
                    }
                    if (si(i,j,0,cS) > zero) {
                        albedo_ice = alb_snow;
                    } else if (si(i,j,0,cM) > Real(0.02)) {
                        albedo_ice = Real(0.42);
                    } else {
                        albedo_ice = alb_ice;
                    }
                } else {
                    albedo_ice = alb_w;
                }
            } else {
                // Default branch; sfc_temp is in Celsius but tested as Kelvin, as in ROMS
                const Real cff1 = alb_s_wet - alb_s_dry;
                const Real cff2 = alb_i_wet - alb_i_dry;
                if (ai > min_ai) {
                    if (si(i,j,0,cS) > zero) {
                        albedo_ice = ((sfc_temp - Real(273.16)) > -one)
                                   ? cff1 * (sfc_temp - Real(272.16)) + alb_s_dry : alb_s_dry;
                    } else {
                        albedo_ice = ((sfc_temp - Real(273.16)) > -one)
                                   ? cff2 * (sfc_temp - Real(272.16)) + alb_i_dry : alb_i_dry;
                    }
                } else {
                    albedo_ice = alb_w;
                }
            }
            f(i,j,0,IceFR::albedo_ice) = albedo_ice;
        });
    }
}
