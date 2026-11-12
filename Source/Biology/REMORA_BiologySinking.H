#ifndef REMORA_BIOLOGY_SINKING_H_
#define REMORA_BIOLOGY_SINKING_H_

#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <cmath>

namespace REMORABiology {

/**
 * Vertical sinking of one particulate tracer in column (i,j), by the ROMS
 * piecewise-parabolic reconstruction with WENO-reconciled interfaces and
 * monotone boundary closures, followed by a CFL-free semi-Lagrangian flux.
 * Shared by the Fennel and ECB kernels, whose ROMS sources carry identical
 * copies.
 *
 * \p nlev is the number of layers, k = 0..nlev-1. On entry \p qc holds the
 * concentration to sink (each model restricts it as its ROMS source does).
 * On exit \p bio holds the sunk concentration and FC(i,j,0) the flux through
 * the bottom face, which the caller returns at the seafloor. \p sink_distance
 * is dtdays * |w|. \p qc, \p bR, \p bL, \p WR and \p WL are cell scratch, and
 * \p FC face scratch, over at least the column.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void sink_column (int i, int j, int nlev, amrex::Real sink_distance,
                  amrex::Array4<amrex::Real> const& bio,
                  amrex::Array4<amrex::Real> const& qc,
                  amrex::Array4<amrex::Real> const& bR,
                  amrex::Array4<amrex::Real> const& bL,
                  amrex::Array4<amrex::Real> const& WR,
                  amrex::Array4<amrex::Real> const& WL,
                  amrex::Array4<amrex::Real> const& FC,
                  amrex::Array4<amrex::Real const> const& Hz,
                  amrex::Array4<amrex::Real const> const& z_w,
                  amrex::Array4<amrex::Real> const& inv_hz,
                  amrex::Array4<amrex::Real> const& inv_hz2,
                  amrex::Array4<amrex::Real> const& inv_hz3) noexcept
{
    using amrex::Real;
    constexpr Real zero = Real(0.0);
    constexpr Real one = Real(1.0);
    constexpr Real two = Real(2.0);
    constexpr Real cff_weno = Real(1.0e-14);
    // The boundary closures below touch layers 0, 1, N-1 and N.
    AMREX_ASSERT(nlev >= 2);
    const int N = nlev - 1;   // top layer

    // Reconstruct vertical profiles as parabolic segments within each grid box.
    for (int iface = 1; iface <= N; ++iface) {
        FC(i,j,iface) = (qc(i,j,iface) - qc(i,j,iface-1)) * inv_hz2(i,j,iface-1);
    }
    for (int k = 1; k < N; ++k) {
        Real dltR = Hz(i,j,k) * FC(i,j,k+1);
        Real dltL = Hz(i,j,k) * FC(i,j,k);
        const Real cff = Hz(i,j,k-1) + two * Hz(i,j,k) + Hz(i,j,k+1);
        const Real cffR = cff * FC(i,j,k+1);
        const Real cffL = cff * FC(i,j,k);

        // Apply PPM monotonicity constraint to prevent oscillations within
        // the grid box.
        if ((dltR * dltL) <= zero) {
            dltR = zero;
            dltL = zero;
        } else if (std::abs(dltR) > std::abs(cffL)) {
            dltR = cffL;
        } else if (std::abs(dltL) > std::abs(cffR)) {
            dltL = cffR;
        }

        // Compute right and left side values of parabolic segments; WR and WL
        // are measures of quadratic variations.
        const Real curv = (dltR - dltL) * inv_hz3(i,j,k);
        dltR -= curv * Hz(i,j,k+1);
        dltL += curv * Hz(i,j,k-1);
        bR(i,j,k) = qc(i,j,k) + dltR;
        bL(i,j,k) = qc(i,j,k) - dltL;
        WR(i,j,k) = (two * dltR - dltL) * (two * dltR - dltL);
        WL(i,j,k) = (dltR - two * dltL) * (dltR - two * dltL);
    }
    // Reconcile interfacial values using a WENO procedure so the whole profile
    // remains monotonic.
    for (int k = 1; k < N - 1; ++k) {
        const Real dltL = amrex::max(cff_weno, WL(i,j,k));
        const Real dltR = amrex::max(cff_weno, WR(i,j,k+1));
        bR(i,j,k) = (dltR * bR(i,j,k) + dltL * bL(i,j,k+1)) / (dltR + dltL);
        bL(i,j,k+1) = bR(i,j,k);
    }
    bR(i,j,N) = qc(i,j,N);         // default strictly monotonic conditions
    bL(i,j,N) = qc(i,j,N);
    bR(i,j,N - 1) = qc(i,j,N);
    bL(i,j,1) = qc(i,j,0);         // bottom grid boxes are re-assumed to be
    bR(i,j,0) = qc(i,j,0);         // piecewise constant
    bL(i,j,0) = qc(i,j,0);

    // Apply monotonicity constraint again, since reconciled interfacial values
    // may cause non-monotonic behavior inside the grid box.
    for (int k = 0; k <= N; ++k) {
        Real dltR = bR(i,j,k) - qc(i,j,k);
        Real dltL = qc(i,j,k) - bL(i,j,k);
        const Real cffR = two * dltR;
        const Real cffL = two * dltL;
        if ((dltR * dltL) < zero) {
            dltR = zero;
            dltL = zero;
        } else if (std::abs(dltR) > std::abs(cffL)) {
            dltR = cffL;
        } else if (std::abs(dltL) > std::abs(cffR)) {
            dltL = cffR;
        }
        bR(i,j,k) = qc(i,j,k) + dltR;
        bL(i,j,k) = qc(i,j,k) - dltL;
    }

    // Compute vertical advective fluxes. The algorithm is free of a CFL
    // criterion by allowing semi-Lagrangian integration bounds to use as many
    // upstream grid boxes as necessary. WL is the z-coordinate of the
    // departure point for grid box interface z_w; FC is the finite-volume flux.
    for (int k = 0; k <= N; ++k) {
        WL(i,j,k) = z_w(i,j,k) + sink_distance;
        WR(i,j,k) = Hz(i,j,k) * qc(i,j,k);
    }
    FC(i,j,N + 1) = zero;          // no-flux boundary condition
    for (int k = 0; k <= N; ++k) {
        int ksource = k;
        Real flux = zero;
        for (int ks = k; ks < N; ++ks) {
            if (WL(i,j,k) > z_w(i,j,ks+1)) {
                ksource = ks + 1;
                flux += WR(i,j,ks);
            }
        }
        // Finalize the flux by adding the fractional part from the source box.
        const Real cu = amrex::min(one, (WL(i,j,k) - z_w(i,j,ksource)) * inv_hz(i,j,ksource));
        flux += Hz(i,j,ksource) * cu *
                (bL(i,j,ksource) + cu * (Real(0.5) * (bR(i,j,ksource) - bL(i,j,ksource)) -
                 (Real(1.5) - cu) * (bR(i,j,ksource) + bL(i,j,ksource) -
                  two * qc(i,j,ksource))));
        FC(i,j,k) = flux;
    }
    for (int k = 0; k <= N; ++k) {
        bio(i,j,k) = qc(i,j,k) + (FC(i,j,k+1) - FC(i,j,k)) * inv_hz(i,j,k);
    }
}

}
#endif
