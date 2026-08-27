#include <REMORA.H>

using namespace amrex;

/**
 * @param[in   ] bx             box to apply climatology on
 * @param[in   ] ioff           offset in x-direction
 * @param[in   ] joff           offset in y-direction
 * @param[inout] var            variable to update
 * @param[in   ] var_old        variable to compare against for nudging
 * @param[in   ] var_clim       climatology value to nudge towards
 * @param[in   ] clim_coeff     nudging time scale (1/s)
 * @param[in   ] Hz             vertical cell height
 * @param[in   ] pm             1/dx
 * @param[in   ] pn             1/dy
 * @param[in   ] dt_lev         time step
 */
void
REMORA::apply_clim_nudg (const Box& bx,
                         int ioff, int joff,
                         const Array4<Real      >& var,
                         const Array4<Real const>& var_old,
                         const Array4<Real const>& var_clim,
                         const Array4<Real const>& clim_coeff,
                         const Array4<Real const>& Hz,
                         const Array4<Real const>& pm,
                         const Array4<Real const>& pn,
                         const Real dt_lev)
{
    BL_PROFILE("REMORA::apply_clim_nudg()");
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        Real cff = Real(0.5) * (clim_coeff(i-ioff,j-joff,k) + clim_coeff(i,j,k));
        if (ioff==1 || joff==1) {
            Real om = two / (pm(i-ioff,j-joff,0)+pm(i,j,0));
            Real on = two / (pn(i-ioff,j-joff,0)+pn(i,j,0));
            cff *= Real(0.5) * (Hz(i-ioff,j-joff,k) + Hz(i,j,k)) * om * on;
        } else {
            cff *= dt_lev;
        }
        var(i,j,k) += cff * (var_clim(i,j,k) - var_old(i,j,k));
    });
}

/**
 * @param[in   ] nodal_bx       box from MFIter::nodaltilebox(dir)
 * @param[in   ] dir            direction the box is nodal in (0 for u, 1 for v)
 * @param[in   ] domain         cell-centered problem domain at this level
 * @param[in   ] is_periodic    whether dir is periodic
 */
Box
REMORA::clim_nudg_momentum_box (const Box& nodal_bx,
                                const int dir,
                                const Box& domain,
                                const bool is_periodic)
{
    Box bx(nodal_bx);

    // ROMS sets IstrU=Istr (JstrV=Jstr) under EW_PERIODIC (NS_PERIODIC), so a periodic
    // direction keeps every face.
    if (is_periodic) {
        return bx;
    }

    // MFIter::nodaltilebox calls surroundingNodes(dir) and only strips the high node when the
    // tile stops short of the valid high end, so a box on the domain high edge legitimately
    // ends at domain.bigEnd(dir)+1 -- the same convention as bx_xhi_face.setSmall(0,dom_hi.x+1)
    // in the boundary conditions. Comparing against the cell-centered domain.bigEnd(dir) would
    // never match there, and would wrongly match a box whose last cell is domain hi minus one.
    const int nodal_lo = domain.smallEnd(dir);
    const int nodal_hi = domain.bigEnd(dir) + 1;

    // Independent tests, not else-if: a box spanning the domain in dir sheds both edge faces.
    if (bx.smallEnd(dir) == nodal_lo) {
        bx.growLo(dir,-1);
    }
    if (bx.bigEnd(dir) == nodal_hi) {
        bx.growHi(dir,-1);
    }

    return bx;
}
