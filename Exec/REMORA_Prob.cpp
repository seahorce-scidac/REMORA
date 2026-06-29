#include "REMORA_Prob.H"

using namespace amrex;

#ifdef ERF_REMORA_FORCE_PROBINIT_LINK
// Force archive extraction of this TU when ERF is linked as a static library
// inside a parent coupled executable and amrex_probinit is weak.
void erf_probinit_link_anchor_func () noexcept {}
#endif

std::unique_ptr<ProblemBase>
amrex_probinit (const amrex_real* problo, const amrex_real* probhi)
{
    return std::make_unique<Problem>(problo, probhi);
}

Problem::Problem(const amrex::Real* /*problo*/, const amrex::Real* /*probhi*/)
{}
