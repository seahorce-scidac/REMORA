#include <ROMSX.H>
#include <Utils.H>

using namespace amrex;

void
ROMSX::uv3dmix  (const Box& bx,
                 Array4<Real> u  , Array4<Real> v,
                 Array4<Real> rufrc, Array4<Real> rvfrc,
                 Array4<Real> visc3d_r,
                 Array4<Real> Hz,
                 Array4<Real> on_r, Array4<Real> om_r,
                 Array4<Real> pn, Array4<Real> pm,
                 int nrhs, int nnew)
{
    // Need to include uv3dmix
    //
    //-----------------------------------------------------------------------
    //  Add in harmonic viscosity s terms.
    //-----------------------------------------------------------------------
    //

}
