#ifndef _REMORA_ErrorTag_H_
#define _REMORA_ErrorTag_H_

#include <AMReX_ErrorList.H>
#include <AMReX_MultiFab.H>
#include <AMReX_TagBox.H>

/**
 * An amrex::AMRErrorTag that does not let the coastline drive refinement.
 *
 * AMReX writes each of its data-dependent tests twice: a plain form that reads every
 * neighbor, and an embedded-boundary form that first asks the geometry which neighbors are
 * really there. The EB form of GRAD takes a difference only across a face that
 * EBCellFlag::isConnected reports as open, and the EB forms of LESS and GREATER test a cell
 * only when it is not covered.
 *
 * REMORA's land/sea mask plays the part of that connectivity. A face between a water cell
 * and a land cell separates two unrelated states -- the tracers on the land side are zeroed
 * every step by advance_3d_ml -- so differencing across it measures where the coast is, not
 * where the flow has structure. This class is AMReX's EB form of those three tests with
 * mskr in place of the EBCellFlags.
 *
 * REMORA used to suppress the same spurious tags after the fact, by appending two derefine
 * criteria that cleared every land cell and then every cell beside one. Because those ran
 * over the whole TagBoxArray after the user's own criteria, they also cleared tags that a
 * static in_box_lo/in_box_hi region had deliberately set: a box drawn across a coast came
 * back refined only on its water part. Guarding each criterion as it is evaluated leaves
 * every other criterion's tags alone, so such a box is now refined whole.
 *
 * One caution on the inheritance: amrex::AMRErrorTag::operator() is not virtual, so the guard
 * lives only on this class's own overload. Keep criteria in a container of REMORAErrorTag, as
 * REMORA::ref_tags does, and call them directly. Anything that converts one to the base --
 * binding an AMRErrorTag const& or copying into an AMRErrorTag -- would compile and quietly
 * run AMReX's unguarded test. The overload deliberately does not pull the base's seven-argument
 * form in with a using-declaration, so a call that omits the mask fails to compile rather than
 * silently skipping it.
 */
class REMORAErrorTag final : public amrex::AMRErrorTag
{
public:
    //! Every base constructor makes sense here unchanged; this class alters only evaluation.
    using amrex::AMRErrorTag::AMRErrorTag;

    /**
     * Apply this criterion to tba, skipping land.
     *
     * Only GRAD, LESS and GREATER consult the mask. BOX and PARSER tag on position alone,
     * USER hands the whole decision to a functor, and VORT and RELGRAD are unreachable from
     * REMORA's inputs; all of those go to amrex::AMRErrorTag::operator() unchanged.
     *
     * @param[inout] tba       tags to update
     * @param[in]    mf        single-component field this criterion tests, with one ghost
     *                         cell when the test is GRAD
     * @param[in]    mskr3d    land/sea mask on mf's BoxArray and DistributionMap, with enough
     *                         lateral ghost cells for mask_lo/mask_hi below plus one more for
     *                         GRAD (so two rings for vorticity's GRAD, one for a tracer's).
     *                         Pass nullptr to evaluate exactly as amrex::AMRErrorTag would,
     *                         which is what a criterion keyed on the mask itself wants.
     * @param[in]    clearval  value marking an untagged cell
     * @param[in]    tagval    value marking a tagged cell
     * @param[in]    time      current time
     * @param[in]    level     level being tagged
     * @param[in]    geom      geometry of that level
     * @param[in]    mask_lo   lowest offset, in rho-cells relative to a value's own index, of
     *                         the cells whose wet/dry state that value depends on. Defaults to
     *                         zero, which is right for a cell-centered state field: it is
     *                         masked in place, so its own cell is the only one that matters.
     *                         The vertical component must be zero, since the mask is constant
     *                         down a column.
     * @param[in]    mask_hi   highest such offset; see mask_lo.
     *
     *                         Two kinds of field need more. A face velocity is stored at a
     *                         cell index but lives on a face, and vert_mean_3d multiplies it
     *                         by msku = mskr(i-1,j)*mskr(i,j) (or mskv in y), so u at a water
     *                         cell whose i-1 neighbor is land is an exact zero that is purely
     *                         a mask artifact. Vorticity is worse: remora_dervort differences
     *                         the cell-centered velocities at i+-1 and j+-1 without masking
     *                         (see the TODO there), and each of those averages two faces, so
     *                         its value depends on the whole 3x3 block of rho-cells.
     *
     *                         Testing such a value where it is not clean would tag the coast
     *                         just as surely as differencing across the coast would.
     */
    void operator() (amrex::TagBoxArray&    tba,
                     const amrex::MultiFab* mf,
                     const amrex::MultiFab* mskr3d,
                     char                   clearval,
                     char                   tagval,
                     amrex::Real            time,
                     int                    level,
                     const amrex::Geometry& geom,
                     const amrex::IntVect&  mask_lo = amrex::IntVect::TheZeroVector(),
                     const amrex::IntVect&  mask_hi = amrex::IntVect::TheZeroVector()) const;
};

#endif
