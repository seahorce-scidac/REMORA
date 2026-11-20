#ifndef REMORA_MASKPOINTS_H_
#define REMORA_MASKPOINTS_H_

#include <AMReX_Dim3.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MFIter.H>
#include <AMReX_iMultiFab.H>

/**
 * The cells of each local box, ghost cells included, at which an integer mask holds a given
 * value: the work list of a pass that writes only those cells, so that it loops over them
 * rather than over the whole box and tests the mask at every cell. The list is built on the
 * host from the mask, which suits masks that depend on the layout alone and are built once
 * per layout, and is kept on the device.
 *
 * Building reads device data on the host, so it must not run inside a BoxLocalLoops scope.
 */
class MaskPoints
{
public:
    MaskPoints () = default;

    /** \brief List the cells of every local box of mask where it equals val */
    void define (const amrex::iMultiFab& mask, int val)
    {
        m_ba = mask.boxArray();
        m_dm = mask.DistributionMap();
        m_ng = mask.nGrowVect();
        const int nlocal = mask.local_size();
        m_offset.assign(nlocal+1, 0);
        amrex::Vector<amrex::Dim3> h_pts;
        amrex::Vector<int> h_mask;
        for (int li = 0; li < nlocal; ++li) {
            const amrex::IArrayBox& fab = mask.atLocalIdx(li);
            const amrex::Box& bx = fab.box();
            h_mask.resize(bx.numPts());
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, fab.dataPtr(), fab.dataPtr() + bx.numPts(),
                             h_mask.begin());
            const auto lo = amrex::lbound(bx);
            const auto hi = amrex::ubound(bx);
            amrex::Long n = 0;
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (h_mask[n++] == val) { h_pts.push_back(amrex::Dim3{i,j,k}); }
            }}}
            m_offset[li+1] = static_cast<int>(h_pts.size());
        }
        m_pts.resize(h_pts.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_pts.begin(), h_pts.end(), m_pts.begin());
    }

    /** \brief Whether the list was built from a mask with fab's layout and ghost cells */
    [[nodiscard]] bool matches (const amrex::FabArrayBase& fab) const
    {
        return !m_offset.empty() && m_ba == fab.boxArray() && m_dm == fab.DistributionMap()
            && m_ng == fab.nGrowVect();
    }

    /** \brief Number of listed cells of the box mfi points at */
    [[nodiscard]] int size (const amrex::MFIter& mfi) const
    {
        const int li = mfi.LocalIndex();
        return m_offset[li+1] - m_offset[li];
    }

    /** \brief The listed cells of the box mfi points at, in device memory */
    [[nodiscard]] const amrex::Dim3* data (const amrex::MFIter& mfi) const
    {
        return m_pts.data() + m_offset[mfi.LocalIndex()];
    }

private:
    amrex::Gpu::DeviceVector<amrex::Dim3> m_pts;
    /** \brief per local box, where its cells start in m_pts; one more entry for the end */
    amrex::Vector<int> m_offset;
    amrex::BoxArray m_ba;
    amrex::DistributionMapping m_dm;
    amrex::IntVect m_ng;
};

#endif
