/**
 * \file remora_subset_compare.cpp
 *
 * Compare a run whose level BoxArrays leave part of the domain out (remora.drop_dry_boxes)
 * with a run on the full layout. amrex_fcompare cannot: it needs the two BoxArrays to cover
 * the same cells. This compares every cell of the first file's boxes, valid cells only, with
 * the same cells of the second, which must cover them.
 *
 *   remora_subset_compare [--abs_tol a] [--rel_tol r] [--chk] <subset> <full>
 *
 * Plotfiles by default: every variable and level the two share. With --chk, two checkpoint
 * directories: every MultiFab under Level_<n>/ that both have (cons, velocities, barotropic
 * fields, free surface, mixing coefficients, ...). The default tolerance is zero, so a pass
 * means bitwise identical. Prints "SUBSET AGREE" and exits 0 on success.
 */
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#include <filesystem>
#include <set>
#include <string>
#include <vector>

using namespace amrex;

namespace {

// Largest |a - b| and |b| over the valid cells of a, b copied onto a's layout
bool compare_mf (const std::string& what, const MultiFab& a, const MultiFab& b_in,
                 Real abs_tol, Real rel_tol)
{
    if (a.nComp() != b_in.nComp() || a.ixType() != b_in.ixType()) {
        Print() << "  " << what << ": component count or index type differs\n";
        return false;
    }
    BoxArray need = a.boxArray();
    if (!b_in.boxArray().contains(need)) {
        Print() << "  " << what << ": the second file does not cover the first's boxes\n";
        return false;
    }
    // Host-accessible copies on a's layout, so the loop below also runs in a GPU build
    MultiFab ah(a.boxArray(), a.DistributionMap(), a.nComp(), 0,
                MFInfo().SetArena(The_Pinned_Arena()));
    MultiFab b (a.boxArray(), a.DistributionMap(), a.nComp(), 0,
                MFInfo().SetArena(The_Pinned_Arena()));
    ah.ParallelCopy(a, 0, 0, a.nComp());
    b.ParallelCopy(b_in, 0, 0, a.nComp());
    Gpu::streamSynchronize();
    bool ok = true;
    for (int n = 0; n < a.nComp(); ++n) {
        Real dmax = 0, bmax = 0;
        Long ndiff = 0, nnan = 0;
        for (MFIter mfi(ah); mfi.isValid(); ++mfi) {
            auto const& x = ah.const_array(mfi);
            auto const& y = b.const_array(mfi);
            const Box bx = mfi.validbox();
            const Dim3 lo = lbound(bx), hi = ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
            for (int i = lo.x; i <= hi.x; ++i) {
                const Real u = x(i,j,k,n), v = y(i,j,k,n);
                // Bitwise, NaN included: two NaNs in the same cell agree
                const bool same = (u == v) || (u != u && v != v);
                if (u != u || v != v) {
                    if (!same) { ++nnan; ++ndiff; }
                    continue;
                }
                const Real d = std::abs(u - v);
                bmax = std::max(bmax, std::abs(v));
                if (!same) {
                    ++ndiff;
                    dmax = std::max(dmax, d);
                }
            }
        }
        ParallelDescriptor::ReduceRealMax(dmax);
        ParallelDescriptor::ReduceRealMax(bmax);
        ParallelDescriptor::ReduceLongSum(ndiff);
        ParallelDescriptor::ReduceLongSum(nnan);
        const bool pass = (nnan == 0) && (dmax <= abs_tol || dmax <= rel_tol * bmax)
                          && (ndiff == 0 || abs_tol > 0 || rel_tol > 0);
        Print() << "  " << what << " comp " << n << ": cells differing " << ndiff
                << ", max |diff| " << dmax << " (max |value| " << bmax << ")"
                << (nnan ? ", NaN mismatch" : "") << (pass ? "" : "  <-- FAIL") << "\n";
        ok = ok && pass;
    }
    return ok;
}

bool compare_plotfiles (const std::string& fa, const std::string& fb, Real abs_tol, Real rel_tol)
{
    PlotFileData pa(fa), pb(fb);
    const int nlev = std::min(pa.finestLevel(), pb.finestLevel()) + 1;
    if (pa.finestLevel() != pb.finestLevel()) {
        Print() << "Finest levels differ: " << pa.finestLevel() << " vs " << pb.finestLevel() << "\n";
        return false;
    }
    const auto& va = pa.varNames();
    const std::set<std::string> vb(pb.varNames().begin(), pb.varNames().end());
    bool ok = true;
    for (int lev = 0; lev < nlev; ++lev) {
        Print() << "level " << lev << ": " << pa.boxArray(lev).size() << " vs "
                << pb.boxArray(lev).size() << " boxes\n";
        for (const auto& v : va) {
            if (vb.count(v) == 0) { continue; }
            const MultiFab a = pa.get(lev, v);
            const MultiFab b = pb.get(lev, v);
            ok = compare_mf(v, a, b, abs_tol, rel_tol) && ok;
        }
    }
    return ok;
}

bool compare_checkpoints (const std::string& fa, const std::string& fb, Real abs_tol, Real rel_tol)
{
    namespace fs = std::filesystem;
    bool ok = true;
    int ncompared = 0;
    for (int lev = 0; ; ++lev) {
        const std::string la = fa + "/Level_" + std::to_string(lev);
        const std::string lb = fb + "/Level_" + std::to_string(lev);
        if (fs::is_directory(la) != fs::is_directory(lb)) {
            Print() << "Level " << lev << " exists in only one of the checkpoints\n";
            return false;
        }
        if (!fs::is_directory(la)) { break; }
        std::vector<std::string> names;
        for (const auto& e : fs::directory_iterator(la)) {
            const std::string f = e.path().filename().string();
            if (f.size() > 2 && f.substr(f.size()-2) == "_H" &&
                fs::exists(lb + "/" + f)) {
                names.push_back(f.substr(0, f.size()-2));
            }
        }
        std::sort(names.begin(), names.end());
        Print() << "level " << lev << ": " << names.size() << " MultiFabs in both\n";
        for (const auto& nm : names) {
            MultiFab a, b;
            VisMF::Read(a, la + "/" + nm);
            VisMF::Read(b, lb + "/" + nm);
            ok = compare_mf(nm, a, b, abs_tol, rel_tol) && ok;
            ++ncompared;
        }
    }
    if (ncompared == 0) {
        Print() << "Nothing to compare in " << fa << " and " << fb << "\n";
        return false;
    }
    return ok;
}
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    int status = 0;
    {
        Real abs_tol = 0, rel_tol = 0;
        bool chk = false;
        std::vector<std::string> files;
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            if (a == "--abs_tol" && i+1 < argc) { abs_tol = std::stod(argv[++i]); }
            else if (a == "--rel_tol" && i+1 < argc) { rel_tol = std::stod(argv[++i]); }
            else if (a == "--chk") { chk = true; }
            else { files.push_back(a); }
        }
        if (files.size() != 2) {
            Print() << "usage: remora_subset_compare [--abs_tol a] [--rel_tol r] [--chk] <subset> <full>\n";
            status = 2;
        } else {
            const bool ok = chk ? compare_checkpoints(files[0], files[1], abs_tol, rel_tol)
                                : compare_plotfiles(files[0], files[1], abs_tol, rel_tol);
            Print() << (ok ? "SUBSET AGREE" : "SUBSET DIFFER") << "\n";
            status = ok ? 0 : 1;
        }
    }
    amrex::Finalize();
    return status;
}
