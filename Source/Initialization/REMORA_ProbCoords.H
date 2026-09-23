#ifndef REMORA_PROBCOORDS_H_
#define REMORA_PROBCOORDS_H_

#include <AMReX_MultiFab.H>

/**
 * \brief Coordinates a prob function reads, on the BoxArray and DistributionMapping of the
 *        MultiFab it fills.
 *
 * REMORA::prob_coords(lev) gives a level's own arrays. The hires_init_level path passes
 * full-domain arrays for a level that does not exist yet.
 *
 * In its own header because REMORA.H and REMORA_prob_common.H include each other.
 */
struct ProbCoords
{
    const amrex::MultiFab* x_r = nullptr;
    const amrex::MultiFab* y_r = nullptr;
    const amrex::MultiFab* z_r = nullptr;
    const amrex::MultiFab* z_w = nullptr;
};

#endif
