#ifndef REMORA_STAGGERED_BILINEAR_H_
#define REMORA_STAGGERED_BILINEAR_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Interpolater.H>
#include <AMReX_Loop.H>
#include <REMORA_NonperiodicGeometry.H>
#include <REMORA_PeriodicImage.H>

#include <cmath>

/**
 * One bilinear value from the four coarse values around it, at weights wx and wy measured from
 * c00 (outside [0, 1] when extrapolating). With a period, the four are first brought within half
 * a period of c00, and the result is then taken to the image nearest the coarse value with the
 * largest weight, so each fine value stays on the branch of its nearest coarse point and a
 * branch cut (the date line) falls between the same points as on the coarse level.
 *
 * @param[in   ] c00     coarse value at (ic, jc)
 * @param[in   ] c10     coarse value at (ic+1, jc)
 * @param[in   ] c01     coarse value at (ic, jc+1)
 * @param[in   ] c11     coarse value at (ic+1, jc+1)
 * @param[in   ] wx      weight of the ic+1 side
 * @param[in   ] wy      weight of the jc+1 side
 * @param[in   ] period  period of the field, or 0 for a field that does not wrap
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
staggered_bilinear_value (amrex::Real c00, amrex::Real c10, amrex::Real c01, amrex::Real c11,
                          amrex::Real wx, amrex::Real wy, amrex::Real period)
{
    const amrex::Real half = amrex::Real(0.5);
    const amrex::Real nearest = (wx < half) ? ((wy < half) ? c00 : c01)
                                            : ((wy < half) ? c10 : c11);
    if (period > amrex::Real(0.0)) {
        c10 = nearest_periodic_image(c10, c00, period);
        c01 = nearest_periodic_image(c01, c00, period);
        c11 = nearest_periodic_image(c11, c00, period);
    }
    const amrex::Real v = (amrex::Real(1.0) - wx) * (amrex::Real(1.0) - wy) * c00
                        + wx * (amrex::Real(1.0) - wy)                      * c10
                        + (amrex::Real(1.0) - wx) * wy                      * c01
                        + wx * wy                                          * c11;
    return (period > amrex::Real(0.0)) ? nearest_periodic_image(v, nearest, period) : v;
}

/**
 * \brief Bilinear coarse-to-fine interpolation of a 2D field at any horizontal staggering
 *        (rho, u, v or psi points), for grid coordinates.
 *
 * The arithmetic reads the coarse field only where a grid file sets it: the domain, plus the
 * one-cell ring outside a physical edge along a cell-centred direction. Points beyond
 * extrapolate linearly from the two nearest values, so a coordinate linear in the index is
 * reproduced exactly. On a NetCDF level 0 the coordinates' ghost cells beyond that ring are
 * never set, which is why AMReX's cell and node bilinear interpolaters, which read whatever the
 * coarse ghost cells hold, don't fit; its face interpolater is also piecewise constant across
 * the face.
 *
 * Used through InterpFromCoarseLevel with a non-periodic coarse Geometry (coordinates do not
 * wrap; see nonperiodic_geometry) and one BCRec per component saying where the ring is: ext_dir
 * at a physical edge, int_dir across a periodic one (see staggered_bilinear_bcrec). Across a
 * periodic edge it reads no coarse ghost cell at all, since those hold the wrapped values from
 * the other side. InterpFromCoarseLevel copies the parent's ghost cells into the coarse patch,
 * so the parent needs ceil(fine ghosts / ratio) + 1 of them.
 *
 * Constructed with a period (360 for a longitude in degrees, 2 pi for an angle in radians),
 * it interpolates a field that wraps (a longitude across the date line) across the cut rather
 * than through the jump, and keeps each fine value on its nearest coarse point's branch (see
 * staggered_bilinear_value).
 */
class StaggeredBilinear
    : public amrex::Interpolater
{
public:
    /**
     * @param[in   ] period  period of the field, or 0 for a field that does not wrap
     */
    explicit StaggeredBilinear (amrex::Real period = amrex::Real(0.0)) : m_period(period) {}

    /**
     * The coarse box a fine box needs, for a uniform refinement ratio
     *
     * @param[in   ] fine   fine box, at the field's staggering
     * @param[in   ] ratio  refinement ratio
     */
    amrex::Box CoarseBox (const amrex::Box& fine, int ratio) override
    {
        return CoarseBox(fine, amrex::IntVect(ratio));
    }

    /**
     * The coarse cells a fine box's stencils can touch: the coarsened box, at the fine box's
     * staggering, plus the far side of each bilinear pair
     *
     * @param[in   ] fine   fine box, at the field's staggering
     * @param[in   ] ratio  refinement ratio
     */
    amrex::Box CoarseBox (const amrex::Box& fine, const amrex::IntVect& ratio) override
    {
        amrex::Box b = amrex::coarsen(fine, ratio);
        b.grow(0, 1);
        b.grow(1, 1);
        return b;
    }

    /**
     * Fill fine_region of fine by bilinear interpolation of crse, extrapolating linearly past
     * the coarse points a grid file sets
     *
     * @param[in   ] crse         coarse patch, as InterpFromCoarseLevel fills it
     * @param[in   ] crse_comp    first coarse component
     * @param[inout] fine         fine FAB
     * @param[in   ] fine_comp    first fine component
     * @param[in   ] ncomp        number of components
     * @param[in   ] fine_region  fine points to fill
     * @param[in   ] ratio        refinement ratio
     * @param[in   ] crse_geom    coarse Geometry; must not be periodic
     * @param[in   ] bcr          per component: ext_dir at a physical edge, int_dir across a
     *                            periodic one
     * @param[in   ] runon        where to run
     */
    void interp (const amrex::FArrayBox& crse, int crse_comp, amrex::FArrayBox& fine,
                 int fine_comp, int ncomp, const amrex::Box& fine_region,
                 const amrex::IntVect& ratio, const amrex::Geometry& crse_geom,
                 const amrex::Geometry& /*fine_geom*/, amrex::Vector<amrex::BCRec> const& bcr,
                 int /*actual_comp*/, int /*actual_state*/, amrex::RunOn runon) override
    {
        // A periodic copy would give a ghost cell beyond a wall two sources at the seam
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!crse_geom.isAnyPeriodic(),
            "StaggeredBilinear needs a non-periodic coarse Geometry (nonperiodic_geometry)");

        const amrex::IntVect stag = fine_region.ixType().toIntVect();
        const amrex::Real sx = amrex::Real(0.5) * amrex::Real(1 - stag[0]);
        const amrex::Real sy = amrex::Real(0.5) * amrex::Real(1 - stag[1]);
        const int rx = ratio[0], ry = ratio[1];

        for (int n = 0; n < ncomp; ++n) {
            amrex::Box filled = amrex::convert(amrex::makeSlab(crse_geom.Domain(), 2, 0), stag);
            for (int dir = 0; dir < 2; ++dir) {
                for (const int side : {bcr[n].lo(dir), bcr[n].hi(dir)}) {
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(side == amrex::BCType::ext_dir ||
                                                     side == amrex::BCType::int_dir,
                        "StaggeredBilinear takes ext_dir or int_dir boundary records only");
                }
                if (stag[dir] != 0) { continue; }
                if (bcr[n].lo(dir) == amrex::BCType::ext_dir) { filled.growLo(dir, 1); }
                if (bcr[n].hi(dir) == amrex::BCType::ext_dir) { filled.growHi(dir, 1); }
            }
            filled &= crse.box();
            AMREX_ALWAYS_ASSERT(filled.length(0) >= 2 && filled.length(1) >= 2);
            const int ilo = filled.smallEnd(0), ihi = filled.bigEnd(0) - 1;
            const int jlo = filled.smallEnd(1), jhi = filled.bigEnd(1) - 1;

            amrex::Array4<amrex::Real const> const& c = crse.const_array(crse_comp + n);
            amrex::Array4<amrex::Real> const& f = fine.array(fine_comp + n);
            const amrex::Real period = m_period;

            // Written out rather than AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG, whose GPU expansion
            // names RunOn and Gpu without the amrex:: this header can't assume
            auto kernel = [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
                // Position in the coarse index space, measured from its first point
                const amrex::Real tx = (amrex::Real(i) + sx) / amrex::Real(rx) - sx;
                const amrex::Real ty = (amrex::Real(j) + sy) / amrex::Real(ry) - sy;
                const int ic = amrex::min(amrex::max(static_cast<int>(std::floor(tx)), ilo), ihi);
                const int jc = amrex::min(amrex::max(static_cast<int>(std::floor(ty)), jlo), jhi);
                const amrex::Real wx = tx - amrex::Real(ic);
                const amrex::Real wy = ty - amrex::Real(jc);
                f(i,j,k) = staggered_bilinear_value(c(ic,jc,k), c(ic+1,jc,k), c(ic,jc+1,k),
                                                    c(ic+1,jc+1,k), wx, wy, period);
            };
            if (runon == amrex::RunOn::Device && amrex::Gpu::inLaunchRegion()) {
                amrex::ParallelFor(fine_region, kernel);
            } else {
                amrex::LoopConcurrentOnCpu(fine_region, kernel);
            }
        }
    }

private:
    amrex::Real m_period;
};

/**
 * The boundary records StaggeredBilinear expects for a coarse level: ext_dir at a physical
 * edge, where a grid file supplies a ring, and int_dir across a periodic one. Build them with
 * this rather than by hand: InterpFromCoarseLevel turns every side of a patch that doesn't reach
 * the domain edge into int_dir before the interpolater sees it, so the interpolater's own check
 * catches a wrong record only on patches at the edge.
 *
 * @param[in   ] crse_geom  the coarse level's Geometry, with its real periodicity
 * @param[in   ] ncomp      number of components
 */
inline amrex::Vector<amrex::BCRec>
staggered_bilinear_bcrec (const amrex::Geometry& crse_geom, int ncomp = 1)
{
    amrex::Vector<amrex::BCRec> bcs(ncomp);
    for (auto& bc : bcs) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            const int type = crse_geom.isPeriodic(dir) ? amrex::BCType::int_dir
                                                       : amrex::BCType::ext_dir;
            bc.setLo(dir, type);
            bc.setHi(dir, type);
        }
    }
    return bcs;
}

#endif
