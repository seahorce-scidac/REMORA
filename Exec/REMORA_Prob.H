#ifndef REMORA_NEWPROB_H_
#define REMORA_NEWPROB_H_

#include <string>

#include "AMReX_REAL.H"
#include "REMORA_ProbCommon.H"

class Problem : public ProblemBase
{
public:
    Problem (const amrex::Real* problo, const amrex::Real* probhi);
};

#endif
