#ifndef _REMORA_ErrorTag_H_
#define _REMORA_ErrorTag_H_

#include <AMReX_ErrorList.H>
#include <AMReX_MultiFab.H>
#include <AMReX_TagBox.H>

/**
 * An amrex::AMRErrorTag that does not let the coastline drive refinement.
 *
 * AMReX writes each data-dependent test twice: a plain form, and an embedded-boundary form
 * that skips a covered cell and differences only across a face EBCellFlag reports as open.
 * REMORA's land/sea mask plays the part of that connectivity -- tracers on the land side are
 * zeroed every step, so a difference taken across a coast measures the coast, not the flow.
 * This class is that EB form with mskr in place of the flags.
 *
 * REMORA used to suppress the same spurious tags afterwards instead, by appending two
 * derefine criteria that cleared every land cell and every cell beside one. Running after the
 * user's own criteria, they also cleared tags a static in_box_lo/in_box_hi region had
 * deliberately set, so a box drawn across a coast came back refined only on its water part.
 * Guarding each criterion as it is evaluated leaves every other criterion's tags alone.
 *
 * The inheritance is private on purpose: AMRErrorTag::operator() is not virtual, so converting
 * a criterion back to the base would compile and quietly run the unguarded test on masked
 * data. Private inheritance makes that ill-formed while still granting access to the protected
 * m_test, m_value and m_info. Only the constructors and Field() are re-exported; the base's
 * seven-argument operator() deliberately is not, so a call that forgets the mask will not
 * compile.
 */
class REMORAErrorTag final : private amrex::AMRErrorTag
{
public:
    //! Every base constructor makes sense here unchanged; this class alters only evaluation.
    using amrex::AMRErrorTag::AMRErrorTag;

    //! The field name a criterion was built with; ErrorEst switches on it to fill the data.
    using amrex::AMRErrorTag::Field;

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
     *                         lateral ghost cells for mask_lo/mask_hi plus one more for GRAD.
     *                         nullptr evaluates exactly as amrex::AMRErrorTag would, which is
     *                         what a criterion keyed on the mask itself wants, and what a run
     *                         with no mask wants since every value is then clean anyway.
     * @param[in]    clearval  value marking an untagged cell
     * @param[in]    tagval    value marking a tagged cell
     * @param[in]    time      current time
     * @param[in]    level     level being tagged
     * @param[in]    geom      geometry of that level
     * @param[in]    mask_lo   lowest offset, in rho-cells relative to a value's own index, of
     *                         the cells whose wet/dry state that value depends on -- not
     *                         always the cell it is stored in, since a face velocity is masked
     *                         by both cells sharing its face and vorticity by a 3x3 block.
     *                         Testing a value where it is not clean tags the coast just as
     *                         differencing across the coast would. ErrorEst sets these per
     *                         field and says why. Zero, the default, means the value is its
     *                         own cell alone; the vertical component must be zero, the mask
     *                         being column-constant.
     * @param[in]    mask_hi   highest such offset; see mask_lo.
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
