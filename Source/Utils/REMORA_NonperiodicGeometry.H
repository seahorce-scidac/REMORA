#ifndef REMORA_NONPERIODIC_GEOMETRY_H_
#define REMORA_NONPERIODIC_GEOMETRY_H_

#include <AMReX_Geometry.H>

/**
 * A copy of geom with no periodic direction, for interpolating coordinates from a coarser
 * level. Coordinates do not wrap, so a periodic copy of the coarse data is wrong for them: a
 * coarse ghost cell that also lies beyond a wall would get a second source, its periodic image,
 * with the winner undefined. Note that a coarse level's coordinate ghost cells across a periodic
 * edge are still filled by a periodic FillBoundary, so they hold values a domain width away;
 * an interpolater must not read them there.
 *
 * @param[in   ] geom   the coarse Geometry
 */
inline amrex::Geometry
nonperiodic_geometry (const amrex::Geometry& geom)
{
    return amrex::Geometry(geom.Domain(), geom.ProbDomain(), geom.Coord(),
                           amrex::Array<int,AMREX_SPACEDIM>{AMREX_D_DECL(0,0,0)});
}

#endif
