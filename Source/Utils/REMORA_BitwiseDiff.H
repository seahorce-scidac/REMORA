#ifndef REMORA_BITWISEDIFF_H_
#define REMORA_BITWISEDIFF_H_

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <bit>
#include <cstdint>
#include <type_traits>

/**
 * The number of points, over all ranks, at which a and b differ in bit pattern in any of the
 * components [0, ncomp), over each box grown by ng; with a mask, only where it equals
 * mask_val. a, b and the mask share a layout. The comparison makes no floating-point
 * operation, so the signalling NaNs of cells never written compare equal to themselves
 * without trapping.
 *
 * @param[in] a         first data
 * @param[in] b         second data
 * @param[in] ncomp     number of components compared, from 0
 * @param[in] ng        ghost cells compared
 * @param[in] mask      where to compare, or nullptr for everywhere
 * @param[in] mask_val  the mask value that selects a point
 */
inline amrex::Long
count_bitwise_diff (const amrex::MultiFab& a, const amrex::MultiFab& b, int ncomp,
                    const amrex::IntVect& ng, const amrex::iMultiFab* mask = nullptr,
                    int mask_val = 1)
{
    using Bits = std::conditional_t<sizeof(amrex::Real) == 8, std::uint64_t, std::uint32_t>;
    auto const& x = a.const_arrays();
    auto const& y = b.const_arrays();
    amrex::Long n = 0;
    if (mask) {
        auto const& m = mask->const_arrays();
        n = amrex::ParReduce(amrex::TypeList<amrex::ReduceOpSum>{}, amrex::TypeList<amrex::Long>{},
                             a, ng, ncomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int c)
            -> amrex::GpuTuple<amrex::Long>
        {
            if (m[box_no](i,j,k) != mask_val) { return { amrex::Long(0) }; }
            const bool same = std::bit_cast<Bits>(x[box_no](i,j,k,c))
                           == std::bit_cast<Bits>(y[box_no](i,j,k,c));
            return { same ? amrex::Long(0) : amrex::Long(1) };
        });
    } else {
        n = amrex::ParReduce(amrex::TypeList<amrex::ReduceOpSum>{}, amrex::TypeList<amrex::Long>{},
                             a, ng, ncomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int c)
            -> amrex::GpuTuple<amrex::Long>
        {
            const bool same = std::bit_cast<Bits>(x[box_no](i,j,k,c))
                           == std::bit_cast<Bits>(y[box_no](i,j,k,c));
            return { same ? amrex::Long(0) : amrex::Long(1) };
        });
    }
    amrex::ParallelDescriptor::ReduceLongSum(n);
    return n;
}

#endif
