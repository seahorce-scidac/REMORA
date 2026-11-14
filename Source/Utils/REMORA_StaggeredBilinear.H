#ifndef REMORA_STAGGERED_BILINEAR_H_
#define REMORA_STAGGERED_BILINEAR_H_

#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>

#include <cmath>

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
 * at a physical edge, int_dir across a periodic one. InterpFromCoarseLevel copies the parent's
 * ghost cells into the coarse patch, so the parent needs ceil(fine ghosts / ratio) + 1 of them.
 */
class StaggeredBilinear
    : public amrex::Interpolater
{
public:
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
                if (stag[dir] == 0 && bcr[n].lo(dir) == amrex::BCType::ext_dir) { filled.growLo(dir, 1); }
                if (stag[dir] == 0 && bcr[n].hi(dir) == amrex::BCType::ext_dir) { filled.growHi(dir, 1); }
            }
            filled &= crse.box();
            AMREX_ALWAYS_ASSERT(filled.length(0) >= 2 && filled.length(1) >= 2);
            const int ilo = filled.smallEnd(0), ihi = filled.bigEnd(0) - 1;
            const int jlo = filled.smallEnd(1), jhi = filled.bigEnd(1) - 1;

            amrex::Array4<amrex::Real const> const& c = crse.const_array(crse_comp + n);
            amrex::Array4<amrex::Real> const& f = fine.array(fine_comp + n);

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
                f(i,j,k) = (amrex::Real(1.0) - wx) * (amrex::Real(1.0) - wy) * c(ic  ,jc  ,k)
                         + wx * (amrex::Real(1.0) - wy)                      * c(ic+1,jc  ,k)
                         + (amrex::Real(1.0) - wx) * wy                      * c(ic  ,jc+1,k)
                         + wx * wy                                          * c(ic+1,jc+1,k);
            };
            if (runon == amrex::RunOn::Device && amrex::Gpu::inLaunchRegion()) {
                amrex::ParallelFor(fine_region, kernel);
            } else {
                amrex::LoopConcurrentOnCpu(fine_region, kernel);
            }
        }
    }
};

//! The one instance, used like AMReX's own interpolaters
inline StaggeredBilinear staggered_bilinear_interp;

/**
 * A copy of geom with no periodic direction, for interpolating coordinates from a coarser
 * level. Coordinates do not wrap: across a periodic edge the coarse ghost cells already hold
 * the continued values, and a periodic copy would give a ghost cell that also lies beyond a
 * wall a second source, its periodic image, with the winner undefined.
 *
 * @param[in   ] geom   the coarse Geometry
 */
inline amrex::Geometry
nonperiodic_geometry (const amrex::Geometry& geom)
{
    return amrex::Geometry(geom.Domain(), geom.ProbDomain(), geom.Coord(),
                           amrex::Array<int,AMREX_SPACEDIM>{AMREX_D_DECL(0,0,0)});
}

#endif
