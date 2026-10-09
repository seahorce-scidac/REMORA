#include <REMORA.H>

using namespace amrex;
using namespace REMORASeaIce;

// Lateral boundary conditions of the ice model go through REMORA's FillPatch and its
// physical BCs, as the ocean's do. The ice state takes the ice_bc, ice_u_bc and ice_v_bc
// BC variables, set from ice.bc_type in init_bcs; the fields ROMS fills with its ocean
// bc_r2d, bc_u2d and bc_v2d take foextrap_bc and u2d/v2d_simple_bc.

void
REMORA::ice_exchange (int lev, MultiFab& mf, int scomp, int ncomp)
{
    mf.FillBoundary(scomp, ncomp, geom[lev].periodicity());
}

/**
 * @param[in]    lev     level
 * @param[inout] mfs     the field, one MultiFab per level
 * @param[in]    comp    component to fill
 * @param[in]    bccomp  BC variable
 */
void
REMORA::ice_fill (int lev, Vector<std::unique_ptr<MultiFab>>& mfs, int comp, int bccomp)
{
    // fill_all = false fills only comp, under the one BC variable bccomp
    FillPatch(lev, t_new[lev], *mfs[lev], GetVecOfPtrs(mfs), bccomp, BdyVars::null, comp,
              false, false);
}

void
REMORA::ice_fill_state (int lev)
{
    for (int n = 0; n < vec_ice_si_r[lev]->nComp(); ++n) { ice_fill(lev, vec_ice_si_r, n, ice_bc()); }
    for (int n = 0; n < vec_ice_si_u[lev]->nComp(); ++n) { ice_fill(lev, vec_ice_si_u, n, ice_u_bc()); }
    for (int n = 0; n < vec_ice_si_v[lev]->nComp(); ++n) { ice_fill(lev, vec_ice_si_v, n, ice_v_bc()); }
}
