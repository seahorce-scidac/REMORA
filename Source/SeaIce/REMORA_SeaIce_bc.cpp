#include <REMORA.H>

using namespace amrex;
using namespace REMORASeaIce;

// Lateral boundary conditions of the ice model: ROMS ice_bc2d.F, ice_tibc.F,
// ice_uibc.F and ice_vibc.F, for the closed and gradient types.
//
// ROMS applies these per tile on the boundary row, then exchanges halos. Here
// the exchange comes first and the boundary row is set over each box's whole
// ghost extent along the side, from interior data the exchange already filled.
// That reproduces what ROMS's later exchange copies into the boundary rows of
// neighboring tiles. As in ROMS, only the first ghost row (or the boundary
// face of a normal velocity) is set; the corners are set where both sides are
// non-periodic.
//
// Index offset: REMORA index = ROMS index - 1, for rho and for u/v points.

namespace {

struct IceDomain {
    int ilo, ihi, jlo, jhi;
    bool xper, yper;
};

IceDomain ice_domain (const Geometry& geom)
{
    const Box& dom = geom.Domain();
    return {dom.smallEnd(0), dom.bigEnd(0), dom.smallEnd(1), dom.bigEnd(1),
            geom.isPeriodic(0), geom.isPeriodic(1)};
}

/** A one-row box at i = ii (or j = jj) clipped to [lo, hi] and to the grown box `gbx`. */
Box column_box (int ii, int lo, int hi, const Box& gbx)
{
    return Box(IntVect(ii, std::max(lo, gbx.smallEnd(1)), 0),
               IntVect(ii, std::min(hi, gbx.bigEnd(1)), 0));
}

Box row_box (int jj, int lo, int hi, const Box& gbx)
{
    return Box(IntVect(std::max(lo, gbx.smallEnd(0)), jj, 0),
               IntVect(std::min(hi, gbx.bigEnd(0)), jj, 0));
}

bool has_cell (const Box& gbx, int i, int j)
{
    return gbx.contains(IntVect(i, j, 0));
}

} // namespace

void
REMORA::ice_exchange (int lev, MultiFab& mf, int scomp, int ncomp)
{
    mf.FillBoundary(scomp, ncomp, geom[lev].periodicity());
}

/**
 * ROMS ice_bc2d: closed sides copy the interior value into the boundary row,
 * times the rho mask. ice_bc2d has no gradient branch, so a gradient side is
 * left as the exchange found it, as in ROMS.
 */
void
REMORA::ice_bc2d (int lev, MultiFab& mf, int comp)
{
    const IceDomain d = ice_domain(geom[lev]);
    const auto bc = seaice_params.bc_type;

    for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        Array4<Real> const& f = mf.array(mfi, comp);
        Array4<Real const> const& mskr = vec_mskr[lev]->const_array(mfi);
        const Box gbx = amrex::grow(mfi.validbox(), IntVect(1,1,0));

        // Along a side, the interior range; the whole ghost extent if the
        // other direction is periodic (the exchange filled its images).
        const int jl = d.yper ? gbx.smallEnd(1) : d.jlo;
        const int jh = d.yper ? gbx.bigEnd(1)   : d.jhi;
        const int il = d.xper ? gbx.smallEnd(0) : d.ilo;
        const int ih = d.xper ? gbx.bigEnd(0)   : d.ihi;

        if (!d.xper && gbx.smallEnd(0) < d.ilo && bc[0] == IceBCType::closed) {
            ParallelFor(column_box(d.ilo-1, jl, jh, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                f(i,j,0) = f(i+1,j,0) * mskr(i,j,0);
            });
        }
        if (!d.xper && gbx.bigEnd(0) > d.ihi && bc[2] == IceBCType::closed) {
            ParallelFor(column_box(d.ihi+1, jl, jh, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                f(i,j,0) = f(i-1,j,0) * mskr(i,j,0);
            });
        }
        if (!d.yper && gbx.smallEnd(1) < d.jlo && bc[1] == IceBCType::closed) {
            ParallelFor(row_box(d.jlo-1, il, ih, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                f(i,j,0) = f(i,j+1,0) * mskr(i,j,0);
            });
        }
        if (!d.yper && gbx.bigEnd(1) > d.jhi && bc[3] == IceBCType::closed) {
            ParallelFor(row_box(d.jhi+1, il, ih, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                f(i,j,0) = f(i,j-1,0) * mskr(i,j,0);
            });
        }

        if (!d.xper && !d.yper) {
            const int ilo = d.ilo, ihi = d.ihi, jlo = d.jlo, jhi = d.jhi;
            if (has_cell(gbx, ilo-1, jlo-1)) {
                ParallelFor(Box(IntVect(ilo-1,jlo-1,0), IntVect(ilo-1,jlo-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = Real(0.5) * (f(i+1,j,0) + f(i,j+1,0));
                });
            }
            if (has_cell(gbx, ihi+1, jlo-1)) {
                ParallelFor(Box(IntVect(ihi+1,jlo-1,0), IntVect(ihi+1,jlo-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = Real(0.5) * (f(i,j+1,0) + f(i-1,j,0));
                });
            }
            if (has_cell(gbx, ilo-1, jhi+1)) {
                ParallelFor(Box(IntVect(ilo-1,jhi+1,0), IntVect(ilo-1,jhi+1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = Real(0.5) * (f(i,j-1,0) + f(i+1,j,0));
                });
            }
            if (has_cell(gbx, ihi+1, jhi+1)) {
                ParallelFor(Box(IntVect(ihi+1,jhi+1,0), IntVect(ihi+1,jhi+1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = Real(0.5) * (f(i,j-1,0) + f(i-1,j,0));
                });
            }
        }
    }
}

/**
 * ROMS ice_tibc, on level linew: the boundary-row enthalpy is the interior
 * hi*ti. A closed side copies ti; a gradient side divides the enthalpy by the
 * boundary-row thickness. Both zero where that thickness is at most min_hi.
 */
void
REMORA::ice_tibc (int lev)
{
    const IceDomain d = ice_domain(geom[lev]);
    const auto bc = seaice_params.bc_type;
    const int linew = ice_lvl.linew;
    const int cH = rcomp(IceR::Hice, linew);
    const int cT = rcomp(IceR::Tice, linew);
    const int cE = rcomp(IceR::Enth, linew);
    const Real min_hi = seaice_params.min_hi;
    constexpr Real eps = Real(1.0e-6);

    MultiFab& si = *vec_ice_si_r[lev];

    for (MFIter mfi(si, false); mfi.isValid(); ++mfi) {
        Array4<Real const> const& hi  = si.const_array(mfi, cH);
        Array4<Real      > const& ti  = si.array(mfi, cT);
        Array4<Real      > const& enth = si.array(mfi, cE);
        Array4<Real const> const& mskr = vec_mskr[lev]->const_array(mfi);
        const Box gbx = amrex::grow(mfi.validbox(), IntVect(1,1,0));

        const int jl = d.yper ? gbx.smallEnd(1) : d.jlo;
        const int jh = d.yper ? gbx.bigEnd(1)   : d.jhi;
        const int il = d.xper ? gbx.smallEnd(0) : d.ilo;
        const int ih = d.xper ? gbx.bigEnd(0)   : d.ihi;

        // Boundary cell (i,j) takes its value from interior cell (ii,jj).
        auto side = [=] AMREX_GPU_DEVICE (int i, int j, int ii, int jj, bool closed)
        {
            enth(i,j,0) = hi(ii,jj,0) * ti(ii,jj,0);
            if (closed) {
                ti(i,j,0) = ti(ii,jj,0);
                enth(i,j,0) *= mskr(i,j,0);
                ti(i,j,0)   *= mskr(i,j,0);
            } else {
                enth(i,j,0) *= mskr(i,j,0);
                ti(i,j,0) = enth(i,j,0) / std::max(hi(i,j,0), eps);
            }
            if (hi(i,j,0) <= min_hi) {
                enth(i,j,0) = zero;
                ti(i,j,0) = zero;
            }
        };

        if (!d.xper && gbx.smallEnd(0) < d.ilo && bc[0] != IceBCType::periodic) {
            const bool closed = (bc[0] == IceBCType::closed);
            ParallelFor(column_box(d.ilo-1, jl, jh, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                side(i, j, i+1, j, closed);
            });
        }
        if (!d.xper && gbx.bigEnd(0) > d.ihi && bc[2] != IceBCType::periodic) {
            const bool closed = (bc[2] == IceBCType::closed);
            ParallelFor(column_box(d.ihi+1, jl, jh, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                side(i, j, i-1, j, closed);
            });
        }
        if (!d.yper && gbx.smallEnd(1) < d.jlo && bc[1] != IceBCType::periodic) {
            const bool closed = (bc[1] == IceBCType::closed);
            ParallelFor(row_box(d.jlo-1, il, ih, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                side(i, j, i, j+1, closed);
            });
        }
        if (!d.yper && gbx.bigEnd(1) > d.jhi && bc[3] != IceBCType::periodic) {
            const bool closed = (bc[3] == IceBCType::closed);
            ParallelFor(row_box(d.jhi+1, il, ih, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
            {
                side(i, j, i, j-1, closed);
            });
        }

        if (!d.xper && !d.yper) {
            // Corner (i,j) averages its two boundary-row neighbors (ia,ja) and (ib,jb).
            auto corner = [=] AMREX_GPU_DEVICE (int i, int j, int ia, int ja, int ib, int jb)
            {
                enth(i,j,0) = Real(0.5) * (enth(ia,ja,0) + enth(ib,jb,0));
                ti(i,j,0)   = Real(0.5) * (ti(ia,ja,0)   + ti(ib,jb,0));
                enth(i,j,0) *= mskr(i,j,0);
                ti(i,j,0)   *= mskr(i,j,0);
                if (hi(i,j,0) <= min_hi) {
                    enth(i,j,0) = zero;
                    ti(i,j,0) = zero;
                }
            };
            const int ilo = d.ilo, ihi = d.ihi, jlo = d.jlo, jhi = d.jhi;
            if (has_cell(gbx, ilo-1, jlo-1)) {
                ParallelFor(Box(IntVect(ilo-1,jlo-1,0), IntVect(ilo-1,jlo-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    corner(i, j, i+1, j, i, j+1);
                });
            }
            if (has_cell(gbx, ihi+1, jlo-1)) {
                ParallelFor(Box(IntVect(ihi+1,jlo-1,0), IntVect(ihi+1,jlo-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    corner(i, j, i, j+1, i-1, j);
                });
            }
            if (has_cell(gbx, ilo-1, jhi+1)) {
                ParallelFor(Box(IntVect(ilo-1,jhi+1,0), IntVect(ilo-1,jhi+1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    corner(i, j, i, j-1, i+1, j);
                });
            }
            if (has_cell(gbx, ihi+1, jhi+1)) {
                ParallelFor(Box(IntVect(ihi+1,jhi+1,0), IntVect(ihi+1,jhi+1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    corner(i, j, i, j-1, i-1, j);
                });
            }
        }
    }
}

/**
 * ROMS ice_uibc. West/east set the boundary face: closed to zero, gradient to
 * the next face in. South/north set the ghost row: closed to gamma2 times the
 * interior (free or no slip), gradient to the interior.
 */
void
REMORA::ice_uibc (int lev, MultiFab& mf, int comp)
{
    const IceDomain d = ice_domain(geom[lev]);
    const auto bc = seaice_params.bc_type;
    const Real gamma2 = seaice_params.gamma2;

    for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        Array4<Real> const& f = mf.array(mfi, comp);
        Array4<Real const> const& msku = vec_msku[lev]->const_array(mfi);
        // Face-centered: faces ilo..ihi+1 are in the domain.
        const Box gbx = amrex::grow(mfi.validbox(), IntVect(1,1,0));

        const int jl = d.yper ? gbx.smallEnd(1) : d.jlo;
        const int jh = d.yper ? gbx.bigEnd(1)   : d.jhi;

        if (!d.xper && gbx.smallEnd(0) <= d.ilo) {
            if (bc[0] == IceBCType::closed) {
                ParallelFor(column_box(d.ilo, jl, jh, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = zero;
                });
            } else if (bc[0] == IceBCType::gradient) {
                ParallelFor(column_box(d.ilo, jl, jh, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = f(i+1,j,0) * msku(i,j,0);
                });
            }
        }
        if (!d.xper && gbx.bigEnd(0) >= d.ihi+1) {
            if (bc[2] == IceBCType::closed) {
                ParallelFor(column_box(d.ihi+1, jl, jh, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = zero;
                });
            } else if (bc[2] == IceBCType::gradient) {
                ParallelFor(column_box(d.ihi+1, jl, jh, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = f(i-1,j,0) * msku(i,j,0);
                });
            }
        }

        // Tangential sides. Closed: ROMS Imin..Imax, every face including the
        // two on the domain edge unless EW-periodic. Gradient: IstrP..Iend.
        const int il_c = d.xper ? gbx.smallEnd(0) : d.ilo;
        const int ih_c = d.xper ? gbx.bigEnd(0)   : d.ihi+1;
        const int il_g = d.xper ? gbx.smallEnd(0) : d.ilo+1;
        const int ih_g = d.xper ? gbx.bigEnd(0)   : d.ihi;

        if (!d.yper && gbx.smallEnd(1) < d.jlo) {
            if (bc[1] == IceBCType::closed) {
                ParallelFor(row_box(d.jlo-1, il_c, ih_c, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = gamma2 * f(i,j+1,0) * msku(i,j,0);
                });
            } else if (bc[1] == IceBCType::gradient) {
                ParallelFor(row_box(d.jlo-1, il_g, ih_g, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = f(i,j+1,0) * msku(i,j,0);
                });
            }
        }
        if (!d.yper && gbx.bigEnd(1) > d.jhi) {
            if (bc[3] == IceBCType::closed) {
                ParallelFor(row_box(d.jhi+1, il_c, ih_c, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = gamma2 * f(i,j-1,0) * msku(i,j,0);
                });
            } else if (bc[3] == IceBCType::gradient) {
                ParallelFor(row_box(d.jhi+1, il_g, ih_g, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = f(i,j-1,0) * msku(i,j,0);
                });
            }
        }

        if (!d.xper && !d.yper) {
            const int ilo = d.ilo, ihi = d.ihi, jlo = d.jlo, jhi = d.jhi;
            if (has_cell(gbx, ilo, jlo-1)) {
                ParallelFor(Box(IntVect(ilo,jlo-1,0), IntVect(ilo,jlo-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = Real(0.5) * (f(i+1,j,0) + f(i,j+1,0));
                });
            }
            if (has_cell(gbx, ihi+1, jlo-1)) {
                ParallelFor(Box(IntVect(ihi+1,jlo-1,0), IntVect(ihi+1,jlo-1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = Real(0.5) * (f(i-1,j,0) + f(i,j+1,0));
                });
            }
            if (has_cell(gbx, ilo, jhi+1)) {
                ParallelFor(Box(IntVect(ilo,jhi+1,0), IntVect(ilo,jhi+1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = Real(0.5) * (f(i+1,j,0) + f(i,j-1,0));
                });
            }
            if (has_cell(gbx, ihi+1, jhi+1)) {
                ParallelFor(Box(IntVect(ihi+1,jhi+1,0), IntVect(ihi+1,jhi+1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = Real(0.5) * (f(i-1,j,0) + f(i,j-1,0));
                });
            }
        }
    }
}

/**
 * ROMS ice_vibc, the mirror of ice_uibc, with two of its quirks kept: a
 * closed west side multiplies the ghost row by gamma2 and the mask instead of
 * setting it from the interior (ice_vibc.F:197), and a gradient north side
 * masks with the interior row's mask (ice_vibc.F:491).
 */
void
REMORA::ice_vibc (int lev, MultiFab& mf, int comp)
{
    const IceDomain d = ice_domain(geom[lev]);
    const auto bc = seaice_params.bc_type;
    const Real gamma2 = seaice_params.gamma2;

    for (MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        Array4<Real> const& f = mf.array(mfi, comp);
        Array4<Real const> const& mskv = vec_mskv[lev]->const_array(mfi);
        // Face-centered: faces jlo..jhi+1 are in the domain.
        const Box gbx = amrex::grow(mfi.validbox(), IntVect(1,1,0));

        const int il = d.xper ? gbx.smallEnd(0) : d.ilo;
        const int ih = d.xper ? gbx.bigEnd(0)   : d.ihi;

        if (!d.yper && gbx.smallEnd(1) <= d.jlo) {
            if (bc[1] == IceBCType::closed) {
                ParallelFor(row_box(d.jlo, il, ih, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = zero;
                });
            } else if (bc[1] == IceBCType::gradient) {
                ParallelFor(row_box(d.jlo, il, ih, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = f(i,j+1,0) * mskv(i,j,0);
                });
            }
        }
        if (!d.yper && gbx.bigEnd(1) >= d.jhi+1) {
            if (bc[3] == IceBCType::closed) {
                ParallelFor(row_box(d.jhi+1, il, ih, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = zero;
                });
            } else if (bc[3] == IceBCType::gradient) {
                ParallelFor(row_box(d.jhi+1, il, ih, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = f(i,j-1,0) * mskv(i,j-1,0);
                });
            }
        }

        const int jl_c = d.yper ? gbx.smallEnd(1) : d.jlo;
        const int jh_c = d.yper ? gbx.bigEnd(1)   : d.jhi+1;
        const int jl_g = d.yper ? gbx.smallEnd(1) : d.jlo+1;
        const int jh_g = d.yper ? gbx.bigEnd(1)   : d.jhi;

        if (!d.xper && gbx.smallEnd(0) < d.ilo) {
            if (bc[0] == IceBCType::closed) {
                ParallelFor(column_box(d.ilo-1, jl_c, jh_c, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = gamma2 * f(i,j,0) * mskv(i,j,0);
                });
            } else if (bc[0] == IceBCType::gradient) {
                ParallelFor(column_box(d.ilo-1, jl_g, jh_g, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = f(i+1,j,0) * mskv(i,j,0);
                });
            }
        }
        if (!d.xper && gbx.bigEnd(0) > d.ihi) {
            if (bc[2] == IceBCType::closed) {
                ParallelFor(column_box(d.ihi+1, jl_c, jh_c, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = gamma2 * f(i-1,j,0) * mskv(i,j,0);
                });
            } else if (bc[2] == IceBCType::gradient) {
                ParallelFor(column_box(d.ihi+1, jl_g, jh_g, gbx), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = f(i-1,j,0) * mskv(i,j,0);
                });
            }
        }

        if (!d.xper && !d.yper) {
            const int ilo = d.ilo, ihi = d.ihi, jlo = d.jlo, jhi = d.jhi;
            if (has_cell(gbx, ilo-1, jlo)) {
                ParallelFor(Box(IntVect(ilo-1,jlo,0), IntVect(ilo-1,jlo,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = Real(0.5) * (f(i,j+1,0) + f(i+1,j,0));
                });
            }
            if (has_cell(gbx, ihi+1, jlo)) {
                ParallelFor(Box(IntVect(ihi+1,jlo,0), IntVect(ihi+1,jlo,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = Real(0.5) * (f(i-1,j,0) + f(i,j+1,0));
                });
            }
            if (has_cell(gbx, ilo-1, jhi+1)) {
                ParallelFor(Box(IntVect(ilo-1,jhi+1,0), IntVect(ilo-1,jhi+1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = Real(0.5) * (f(i,j-1,0) + f(i+1,j,0));
                });
            }
            if (has_cell(gbx, ihi+1, jhi+1)) {
                ParallelFor(Box(IntVect(ihi+1,jhi+1,0), IntVect(ihi+1,jhi+1,0)), [=] AMREX_GPU_DEVICE (int i, int j, int)
                {
                    f(i,j,0) = Real(0.5) * (f(i,j-1,0) + f(i-1,j,0));
                });
            }
        }
    }
}
