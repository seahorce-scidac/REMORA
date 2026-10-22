#ifndef REMORA_PC_H_
#define REMORA_PC_H_

/**
 * @file REMORA_PC.H
 * @brief Lagrangian particles, advected in the terrain-following coordinates of the grid.
 *
 * @section particle_coords Vertical coordinate
 *
 * Let \f$z_k(x,y)\f$, \f$k = 0,\dots,N_z\f$, be the heights of the w-levels of the column
 * at horizontal position \f$(x,y)\f$: \f$z_0 = -h\f$ is the bed and \f$z_{N_z} = \zeta\f$ the
 * free surface, with \f$N_z\f$ the number of vertical cells. Each level height is the
 * bilinear interpolation of the nodal heights `z_phys_nd` (see level_height()):
 * \f[
 *     z_k(x,y) = \sum_{a,b \in \{0,1\}} \lambda_a(\xi)\, \lambda_b(\eta)\, z^{nd}_{i+a,\,j+b,\,k},
 *     \qquad \lambda_0(t) = 1-t, \quad \lambda_1(t) = t,
 * \f]
 * where \f$i = \lfloor (x - x_{lo})/\Delta x \rfloor\f$, \f$\xi = (x - x_{lo})/\Delta x - i\f$,
 * and likewise \f$j\f$ and \f$\eta\f$ in \f$y\f$.
 *
 * A particle's vertical position is described by a fractional w-level index
 * \f$s = k + f \in [0, N_z]\f$, \f$0 \le f \le 1\f$, which places it the fraction \f$f\f$ of the
 * way from level \f$k\f$ to level \f$k+1\f$:
 * \f[
 *     z = Z(x,y,s) = z_k(x,y) + f\,\bigl(z_{k+1}(x,y) - z_k(x,y)\bigr).
 * \f]
 * level_index_to_height() evaluates \f$Z\f$ and height_to_level_index() inverts it. The
 * particle stores its height \f$z\f$ in `pos(2)` and the cell \f$k\f$ in its integer component,
 * which the assignor uses to locate it.
 *
 * @section particle_eqs Equations of motion
 *
 * In terrain-following coordinates the particle obeys
 * \f[
 *     \frac{dx}{dt} = u, \qquad \frac{dy}{dt} = v, \qquad
 *     \frac{ds}{dt} = \frac{\Omega}{H_z},
 * \f]
 * where \f$u, v\f$ are the horizontal velocities, \f$\Omega\f$ is the velocity across the
 * w-levels (ROMS omega, in m/s; REMORA::compute_omega()), and
 * \f$H_z = z_{k+1}(x,y) - z_k(x,y)\f$ is the thickness of the particle's cell.
 *
 * The physical vertical velocity follows from the chain rule applied to \f$z = Z(x,y,s,t)\f$:
 * \f[
 *     w = \frac{dz}{dt} = \Omega + u \left.\frac{\partial z}{\partial x}\right|_s
 *       + v \left.\frac{\partial z}{\partial y}\right|_s
 *       + \left.\frac{\partial z}{\partial t}\right|_s .
 * \f]
 * Advancing \f$s\f$ instead of \f$z\f$ has two consequences:
 *  - The along-level terms \f$u\,\partial_x z|_s + v\,\partial_y z|_s\f$ never have to be
 *    computed. Moving the particle horizontally at fixed \f$s\f$ and re-evaluating \f$Z\f$ at
 *    the new position raises or lowers it by exactly the change in level height.
 *  - \f$\Omega = 0\f$ at the bed (\f$s = 0\f$) and the free surface (\f$s = N_z\f$), so the
 *    continuous motion cannot leave \f$[0, N_z]\f$. Only time-discretization error can carry a
 *    particle past either boundary, and it is reflected back (see REMORAPC::AdvectWithFlow()).
 *    Integrating \f$z\f$ with \f$w\f$ instead keeps particles in the column only if the
 *    interpolated \f$w\f$ exactly cancels \f$u \cdot \nabla h\f$ at the bed, which it does not.
 *
 * \f$\Omega\f$ is built from the time-averaged mass fluxes that the tracers are advected with,
 * so particles move with the same discrete, volume-conserving flow as the tracers.
 *
 * @section particle_limits Limitation
 *
 * The particle stores \f$z\f$, not \f$s\f$, so \f$s\f$ is recovered from \f$z\f$ on the current
 * grid at the start of every step. While the free surface moves between steps, a particle
 * therefore keeps its height rather than its level, which drops the level-stretching term
 * \f$\partial_t z|_s = \partial_t \zeta \, (z + h)/(\zeta + h)\f$ from its motion. Storing
 * \f$s\f$ as a particle component would remove this, at the cost of a different checkpoint
 * layout.
 */

#ifdef REMORA_USE_PARTICLES

#include <string>
#include <REMORA_Constants.H>
#include <AMReX_Particles.H>

struct REMORAParticlesIntIdxAoS
{
    enum {
        k = 0,
        ncomps
    };
};

struct REMORAParticlesRealIdxAoS
{
    enum {
        ncomps = 0
    };
};

struct REMORAParticlesIntIdxSoA
{
    enum {
        ncomps = 0
    };
};

struct REMORAParticlesRealIdxSoA
{
    enum {
        vx = 0,
        vy,
        vz,
        mass,
        ncomps
    };
};

namespace REMORAParticleInitializations
{
    /* list of particle initializations */
    const std::string init_box_uniform  = "box";
}

namespace REMORAParticleNames
{
    const std::string tracers = "tracer_particles";
    const std::string hydro = "hydro_particles";
}

struct REMORAParticlesAssignor
{
    template <typename P>
    AMREX_GPU_HOST_DEVICE
    amrex::IntVect operator() (  P const& p,
                                 amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                                 amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                                 const amrex::Box& domain ) const noexcept
    {
        amrex::IntVect iv(
            AMREX_D_DECL(   int(amrex::Math::floor((p.pos(0)-plo[0])*dxi[0])),
                            int(amrex::Math::floor((p.pos(1)-plo[1])*dxi[1])),
                            p.idata(REMORAParticlesIntIdxAoS::k) ) );
        iv[0] += domain.smallEnd()[0];
        iv[1] += domain.smallEnd()[1];
        return iv;
    }
};

// The helpers below implement the coordinate described in the file documentation above:
// a fractional w-level index s in [0,nz], with s = 0 the bed and s = nz the free surface.

/**
 * @brief Height \f$z_k(x,y)\f$ of w-level a_k at horizontal position (a_x, a_y).
 *
 * Bilinear in the nodal heights a_z_nd (`z_phys_nd`), whose \f$(i,j,k)\f$ is the low corner
 * of cell \f$(i,j,k)\f$; see @ref particle_coords.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real level_height (  amrex::ParticleReal a_x, amrex::ParticleReal a_y, int a_k,
                            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_plo,
                            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_dxi,
                            amrex::Array4<amrex::Real const> const& a_z_nd )
{
    const auto lx = static_cast<amrex::Real>((a_x-a_plo[0])*a_dxi[0]);
    const auto ly = static_cast<amrex::Real>((a_y-a_plo[1])*a_dxi[1]);
    const int i = static_cast<int>(amrex::Math::floor(lx));
    const int j = static_cast<int>(amrex::Math::floor(ly));
    const amrex::Real sx = lx - static_cast<amrex::Real>(i);
    const amrex::Real sy = ly - static_cast<amrex::Real>(j);
    return a_z_nd(i  ,j  ,a_k) * (one-sx) * (one-sy) +
           a_z_nd(i+1,j  ,a_k) *      sx  * (one-sy) +
           a_z_nd(i  ,j+1,a_k) * (one-sx) *      sy  +
           a_z_nd(i+1,j+1,a_k) *      sx  *      sy;
}

/**
 * @brief Fractional w-level index \f$s = k + f\f$ of height a_z in the column at (a_x, a_y).
 *
 * Searches up or down from cell a_k, the cell the particle was last in, for the \f$k\f$ with
 * \f$z_k \le z \le z_{k+1}\f$, and returns \f$k + f\f$ with
 * \f[
 *     f = \min\!\left(1, \max\!\left(0, \frac{z - z_k}{z_{k+1} - z_k}\right)\right).
 * \f]
 * A height below the bed or above the free surface is clamped onto it.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real height_to_level_index ( amrex::ParticleReal a_x, amrex::ParticleReal a_y,
                                    amrex::ParticleReal a_z, int a_k, int a_nz,
                                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_plo,
                                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_dxi,
                                    amrex::Array4<amrex::Real const> const& a_z_nd )
{
    const auto z = static_cast<amrex::Real>(a_z);
    int k = amrex::max(0, amrex::min(a_k, a_nz-1));
    amrex::Real zlo = level_height(a_x, a_y, k  , a_plo, a_dxi, a_z_nd);
    amrex::Real zhi = level_height(a_x, a_y, k+1, a_plo, a_dxi, a_z_nd);
    while (z > zhi && k < a_nz-1) {
        ++k;
        zlo = zhi;
        zhi = level_height(a_x, a_y, k+1, a_plo, a_dxi, a_z_nd);
    }
    while (z < zlo && k > 0) {
        --k;
        zhi = zlo;
        zlo = level_height(a_x, a_y, k, a_plo, a_dxi, a_z_nd);
    }
    const amrex::Real f = amrex::max(zero, amrex::min(one, (z - zlo) / (zhi - zlo)));
    return static_cast<amrex::Real>(k) + f;
}

/**
 * @brief Cell \f$k = \min(\lfloor s \rfloor, N_z - 1)\f$ holding fractional w-level index a_s.
 *
 * The free surface \f$s = N_z\f$ belongs to the top cell.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int level_index_cell (amrex::Real a_s, int a_nz)
{
    return amrex::max(0, amrex::min(static_cast<int>(amrex::Math::floor(a_s)), a_nz-1));
}

/**
 * @brief Height \f$Z(x,y,s) = z_k + f\,(z_{k+1} - z_k)\f$ of fractional w-level index a_s.
 *
 * The inverse of height_to_level_index(), with \f$k\f$ from level_index_cell() and
 * \f$f = s - k\f$.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real level_index_to_height ( amrex::ParticleReal a_x, amrex::ParticleReal a_y,
                                    amrex::Real a_s, int a_nz,
                                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_plo,
                                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_dxi,
                                    amrex::Array4<amrex::Real const> const& a_z_nd )
{
    const int k = level_index_cell(a_s, a_nz);
    const amrex::Real f = a_s - static_cast<amrex::Real>(k);
    const amrex::Real zlo = level_height(a_x, a_y, k  , a_plo, a_dxi, a_z_nd);
    const amrex::Real zhi = level_height(a_x, a_y, k+1, a_plo, a_dxi, a_z_nd);
    return zlo + f * (zhi - zlo);
}

/**
 * @brief Trilinear interpolation of a_arr to (a_x, a_y) and fractional w-level index a_s.
 *
 * Interpolation is in index space. a_nodal gives a_arr's staggering, and each direction's
 * index coordinate is shifted by \f$\frac{1}{2}\f$ where the data are cell-centered:
 * \f[
 *     l_x = \frac{x - x_{lo}}{\Delta x} - \frac{1}{2} (1 - \nu_x), \qquad
 *     l_y = \frac{y - y_{lo}}{\Delta y} - \frac{1}{2} (1 - \nu_y), \qquad
 *     l_z = \min\!\bigl(k_{max}, \max\bigl(0,\ s - \frac{1}{2} (1 - \nu_z)\bigr)\bigr),
 * \f]
 * with \f$\nu_d \in \{0,1\}\f$ the nodality in direction \f$d\f$, and \f$k_{max} = N_z\f$ for
 * z-nodal data (on the w-levels) or \f$N_z - 1\f$ for cell-centered data. The clamp in \f$l_z\f$
 * holds cell-centered data constant over the half cells at the bed and the surface, so the
 * vertical stencil never leaves the valid levels.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real interpolate_to_level_index ( amrex::ParticleReal a_x, amrex::ParticleReal a_y,
                                         amrex::Real a_s, int a_nz, amrex::IntVect const& a_nodal,
                                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_plo,
                                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_dxi,
                                         amrex::Array4<amrex::Real const> const& a_arr )
{
    const amrex::Real lx = static_cast<amrex::Real>((a_x-a_plo[0])*a_dxi[0])
                         - static_cast<amrex::Real>(!a_nodal[0]) * myhalf;
    const amrex::Real ly = static_cast<amrex::Real>((a_y-a_plo[1])*a_dxi[1])
                         - static_cast<amrex::Real>(!a_nodal[1]) * myhalf;
    const int i0 = static_cast<int>(amrex::Math::floor(lx));
    const int j0 = static_cast<int>(amrex::Math::floor(ly));
    const amrex::Real sx[] = {one - (lx - static_cast<amrex::Real>(i0)), lx - static_cast<amrex::Real>(i0)};
    const amrex::Real sy[] = {one - (ly - static_cast<amrex::Real>(j0)), ly - static_cast<amrex::Real>(j0)};

    const int kmax = a_nodal[2] ? a_nz : a_nz-1;
    const amrex::Real lz = amrex::max(zero, amrex::min(static_cast<amrex::Real>(kmax),
                                      a_s - static_cast<amrex::Real>(!a_nodal[2]) * myhalf));
    const int k0 = amrex::min(static_cast<int>(amrex::Math::floor(lz)), amrex::max(kmax-1, 0));
    const int k1 = amrex::min(k0+1, kmax);
    const amrex::Real sz[] = {one - (lz - static_cast<amrex::Real>(k0)), lz - static_cast<amrex::Real>(k0)};

    amrex::Real val = zero;
    for (int jj = 0; jj < 2; ++jj) {
        for (int ii = 0; ii < 2; ++ii) {
            val += sx[ii] * sy[jj] * ( sz[0] * a_arr(i0+ii, j0+jj, k0)
                                     + sz[1] * a_arr(i0+ii, j0+jj, k1) );
        }
    }
    return val;
}

class REMORAPC : public amrex::ParticleContainer<  REMORAParticlesRealIdxAoS::ncomps,  // AoS real attributes
                                                REMORAParticlesIntIdxAoS::ncomps,   // AoS integer attributes
                                                REMORAParticlesRealIdxSoA::ncomps,  // SoA real attributes
                                                REMORAParticlesIntIdxSoA::ncomps,   // SoA integer attributes
                                                amrex::DefaultAllocator,
                                                REMORAParticlesAssignor >
{
    public:

        /*! Constructor */
        REMORAPC ( amrex::ParGDBBase* a_gdb,
                const std::string& a_name = "particles" )
            : amrex::ParticleContainer< REMORAParticlesRealIdxAoS::ncomps,   // AoS real attributes
                                        REMORAParticlesIntIdxAoS::ncomps,    // AoS integer attributes
                                        REMORAParticlesRealIdxSoA::ncomps,   // SoA real attributes
                                        REMORAParticlesIntIdxSoA::ncomps,    // SoA integer attributes
                                        amrex::DefaultAllocator,
                                        REMORAParticlesAssignor> (a_gdb)
        {
            BL_PROFILE("REMORAPCPC::REMORAPC()");
            m_name = a_name;
            readInputs();
        }

        /*! Constructor */
        REMORAPC ( const amrex::Geometry&            a_geom,
                const amrex::DistributionMapping& a_dmap,
                const amrex::BoxArray&            a_ba,
                const std::string&                a_name = "particles" )
            : amrex::ParticleContainer< REMORAParticlesRealIdxAoS::ncomps,   // AoS real attributes
                                        REMORAParticlesIntIdxAoS::ncomps,    // AoS real attributes
                                        REMORAParticlesRealIdxSoA::ncomps,   // SoA real attributes
                                        REMORAParticlesIntIdxSoA::ncomps,    // SoA integer attributes
                                        amrex::DefaultAllocator,
                                        REMORAParticlesAssignor> ( a_geom, a_dmap, a_ba )
        {
            BL_PROFILE("REMORAPCPC::REMORAPC()");
            m_name = a_name;
            readInputs();
        }

        /*! Initialize particles in domain */
        virtual void InitializeParticles (const std::unique_ptr<amrex::MultiFab>& a_ptr = nullptr);

        /*! Evolve particles for one time step */
        virtual void EvolveParticles (   int a_lev,
                                         amrex::Real a_dt_lev,
                                         amrex::Vector<amrex::MultiFab const*>& a_flow_vel,
                                         const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& a_z_phys_nd );

        /*! Get real-type particle attribute names */
        virtual amrex::Vector<std::string> varNames () const
        {
            BL_PROFILE("REMORAPCPC::varNames()");
            return {AMREX_D_DECL("xvel","yvel","zvel"),"mass"};
        }

        /*! Get real-type particle attribute names */
        virtual amrex::Vector<std::string> meshPlotVarNames () const
        {
            BL_PROFILE("REMORAPCPC::varNames()");
            return {"mass_density"};
        }

        /*! Advances the particles one step with the midpoint method, in terrain-following
         *  coordinates (documented at the definition, REMORA_PC_Evolve.cpp) */
        virtual void AdvectWithFlow (  int a_lev,
                                       amrex::Real a_dt,
                                       amrex::Vector<amrex::MultiFab const*>& a_umac,
                                       const std::unique_ptr<amrex::MultiFab>& a_z_height );

        /*! Compute mass density */
        virtual void massDensity ( amrex::MultiFab&, const int&, const int& a_comp = 0) const;

        /*! Compute mesh variable from particles */
        virtual void computeMeshVar(    const std::string&  a_var_name,
                                        amrex::MultiFab&    a_mf,
                                        const int           a_lev) const
        {
            if (a_var_name == "mass_density") {
                massDensity( a_mf, a_lev );
            } else {
                a_mf.setVal(zero);
            }
        }

        /*! Specify if particles should advect with flow */
        inline void setAdvectWithFlow (bool a_flag)
        {
            BL_PROFILE("REMORAPCPC::setAdvectWithFlow()");
            m_advect_w_flow = a_flag;
        }

        // the following functions should ideally be private or protected, but need to be
        // public due to CUDA extended lambda capture rules

        /*! Default particle initialization */
        void initializeParticlesUniformDistributionInBox (const std::unique_ptr<amrex::MultiFab>& a_ptr,
                                                          const amrex::RealBox& particle_box);

    protected:

        bool m_advect_w_flow;               /*!< advect with flow velocity */

        amrex::RealBox m_particle_box;      /*!< box within which to place particles */

        std::string m_name;                 /*!< name of this particle species */

        std::string m_initialization_type;  /*!< initial particle distribution type */
        int m_ppc_init;                     /*!< initial number of particles per cell */

        /*! read inputs from file */
        virtual void readInputs ();

    private:

        bool place_randomly_in_cells; /*!< place particles at random positions? */
};

#endif
#endif
