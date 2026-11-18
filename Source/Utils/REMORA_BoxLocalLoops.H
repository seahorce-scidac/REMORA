#ifndef REMORA_BOX_LOCAL_LOOPS_H_
#define REMORA_BOX_LOCAL_LOOPS_H_

#include <AMReX_FabArrayBase.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuDevice.H>

#include <algorithm>
#include <utility>

/**
 * Scope for consecutive MFIter loops over FabArrays with the layout of fa, in which each box
 * reads and writes only its own data, valid and ghost cells alike.
 *
 * On GPU every MFIter loop waits for its stream before it starts and for every stream it used
 * when it ends, so that the next loop or exchange may read what any box wrote. MFIter sends
 * all the work of one box to the same stream, where it runs in order, so loops that never read
 * another box's data need none of that between them. Inside this scope they skip it: the scope
 * waits for the current stream on entry and for the streams of fa's boxes on exit, as its
 * first and last loop would.
 *
 * Nothing inside may exchange ghost cells or read another box's data in any other way, read
 * device data on the host, or free memory that a kernel may still read: a temporary MultiFab
 * or FArrayBox there must come from The_Async_Arena. On CPU the scope does nothing.
 */
class BoxLocalLoops
{
public:
    explicit BoxLocalLoops (const amrex::FabArrayBase& fa)
        : m_nstreams(std::min(fa.local_size(), amrex::Gpu::numGpuStreams()))
    {
        amrex::Gpu::streamSynchronize();
    }

    ~BoxLocalLoops ()
    {
#ifdef AMREX_USE_GPU
        for (int i = 0; i < m_nstreams; ++i) {
            amrex::Gpu::Device::setStreamIndex(i);
            amrex::Gpu::streamSynchronize();
        }
        amrex::Gpu::Device::resetStreamIndex();
#endif
    }

    BoxLocalLoops (const BoxLocalLoops&) = delete;
    BoxLocalLoops& operator= (const BoxLocalLoops&) = delete;

private:
    amrex::Gpu::NoSyncRegion m_no_sync;
    int m_nstreams;
};

/**
 * Run f, a sequence of ghost-cell exchanges, with one wait for the GPU stream at the end
 * instead of one after each exchange. Every stream must be idle on entry, as after an MFIter
 * loop or a BoxLocalLoops scope; the exchanges then run in order on the current stream. On CPU
 * this just runs f.
 */
template <class F>
void exchange_then_sync (F&& f)
{
    {
        [[maybe_unused]] amrex::Gpu::NoSyncRegion no_sync;
        std::forward<F>(f)();
    }
    amrex::Gpu::streamSynchronize();
}

#endif
