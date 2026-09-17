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
     * @param[in]    mskr3d    land/sea mask on mf's BoxArray and DistributionMap, with at
     *                         least one lateral ghost cell. Pass nullptr to evaluate exactly
     *                         as amrex::AMRErrorTag would, which is what a criterion keyed on
     *                         the mask itself wants.
     * @param[in]    clearval  value marking an untagged cell
     * @param[in]    tagval    value marking a tagged cell
     * @param[in]    time      current time
     * @param[in]    level     level being tagged
     * @param[in]    geom      geometry of that level
     * @param[in]    field_radius  how far from its own cell mf's value was computed from.
     *                         Zero for a state field, which is just itself, and the default.
     *                         One for vorticity: remora_dervort centers a difference on
     *                         i+-1 and j+-1 and does not mask it, so a water cell touching
     *                         land carries the land-side velocity inside its own value and
     *                         testing it would tag the coast as surely as differencing
     *                         across it would.
     */
    void operator() (amrex::TagBoxArray&    tba,
                     const amrex::MultiFab* mf,
                     const amrex::MultiFab* mskr3d,
                     char                   clearval,
                     char                   tagval,
                     amrex::Real            time,
                     int                    level,
                     const amrex::Geometry& geom,
                     int                    field_radius = 0) const;
};

#endif
