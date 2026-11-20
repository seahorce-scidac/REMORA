#ifndef REMORA_BOX_LOCAL_LOOPS_H_
#define REMORA_BOX_LOCAL_LOOPS_H_

#include <AMReX_GpuControl.H>
#include <AMReX_GpuDevice.H>

#include <utility>

/**
 * Scope for consecutive MFIter loops in which each box reads and writes only its own data,
 * valid and ghost cells alike.
 *
 * On GPU every MFIter loop waits for its stream before it starts and for every stream it used
 * when it ends, so that the next loop or exchange may read what any box wrote. MFIter sends
 * all the work of one box to the same stream, where it runs in order, so loops that never read
 * another box's data need none of that between them. Inside this scope they skip it: the scope
 * waits for the current stream on entry and for every stream on exit.
 *
 * MFIter picks a box's stream from its local index, so every loop inside must iterate
 * FabArrays with the same BoxArray (up to index type) and DistributionMapping, untiled on GPU;
 * otherwise one box's work could land on different streams in different loops.
 *
 * Nothing inside may exchange ghost cells or read another box's data in any other way, read
 * device data on the host, or free memory that a kernel on another stream may still read. A
 * temporary FArrayBox allocated from The_Async_Arena inside a box's own MFIter iteration is
 * safe, since it is allocated and freed in that box's stream order; a temporary MultiFab or
 * FArrayBox made outside the loop is not, since it lives in the order of the current stream
 * while box b runs on stream b. On CPU the scope does nothing.
 */
class BoxLocalLoops
{
public:
    BoxLocalLoops () { amrex::Gpu::streamSynchronize(); }
    ~BoxLocalLoops () { amrex::Gpu::streamSynchronizeAll(); }

    BoxLocalLoops (const BoxLocalLoops&) = delete;
    BoxLocalLoops& operator= (const BoxLocalLoops&) = delete;

private:
    amrex::Gpu::NoSyncRegion m_no_sync;
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
