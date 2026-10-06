#include <REMORA.H>

using namespace amrex;
using namespace REMORASeaIce;

/**
 * One ice step per baroclinic step: ROMS seaice.F. Called from setup_step
 * after the bulk flux and the surface flux assembly, which is where ROMS's
 * main3d calls it (after bulk_flux and set_vbc, before the ocean step), so
 * the stresses and fluxes the ice sets are the ones the ocean uses this step.
 *
 * The time-level indices are swapped in ROMS's order, and nothing copies the
 * ice tracers from liold to linew unless advection runs: with advection off,
 * the thermodynamics works on the state of two steps back, as in ROMS.
 *
 * @param[in] lev    level to operate on
 * @param[in] dt_lev ice time step, ROMS dtice = dt
 */
void
REMORA::advance_seaice (int lev, Real dt_lev)
{
    BL_PROFILE("REMORA::advance_seaice()");

    ice_lvl.liold = ice_lvl.linew;
    ice_lvl.linew = 1 - ice_lvl.liold;

    if (seaice_params.copy_old_to_new && !seaice_params.do_advect) {
        MultiFab& sir = *vec_ice_si_r[lev];
        for (int var : {IceR::Aice, IceR::Hice, IceR::Hmel, IceR::Hsno, IceR::Iage,
                        IceR::Tice, IceR::Enth, IceR::Hage}) {
            MultiFab::Copy(sir, sir, rcomp(var, ice_lvl.liold), rcomp(var, ice_lvl.linew),
                           1, sir.nGrowVect());
        }
    }

    ice_spdiw(lev);
    ice_vbc(lev);

    if (seaice_params.do_momentum) {
        amrex::Abort("ice.do_momentum is not implemented yet");
    }

    if (seaice_params.do_advect) {
        amrex::Abort("ice.do_advect is not implemented yet");
    }

    if (seaice_params.do_thermo) {
        ice_thermo(lev, dt_lev);
        ice_limit(lev);
    }
}
