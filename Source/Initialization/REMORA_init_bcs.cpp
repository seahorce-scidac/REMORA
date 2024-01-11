#include <AMReX_Vector.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ParmParse.H>

#include <REMORA.H>

using namespace amrex;

void REMORA::init_bcs ()
{
    auto f = [this] (std::string const& bcid, Orientation ori)
    {
        // These are simply defaults for Dirichlet faces -- they should be over-written below
        m_bc_extdir_vals[BCVars::Temp_bc_comp  ][ori] = 1.e19;
        m_bc_extdir_vals[BCVars::Salt_bc_comp  ][ori] = 1.e20;
        m_bc_extdir_vals[BCVars::Scalar_bc_comp][ori] = 1.e21;

        m_bc_extdir_vals[BCVars::xvel_bc][ori] = 0.0; // default
        m_bc_extdir_vals[BCVars::yvel_bc][ori] = 0.0;
        m_bc_extdir_vals[BCVars::zvel_bc][ori] = 0.0;

        ParmParse pp(bcid);
        std::string bc_type_in = "null";
        pp.query("type", bc_type_in);
        std::string bc_type = amrex::toLower(bc_type_in);

        if (bc_type == "symmetry")
        {
            phys_bc_type[ori] = REMORA_BC::symmetry;
            domain_bc_type[ori] = "Symmetry";
        }
        else if (bc_type == "outflow")
        {
            phys_bc_type[ori] = REMORA_BC::outflow;
            domain_bc_type[ori] = "Outflow";
        }
        else if (bc_type == "inflow")
        {
            phys_bc_type[ori] = REMORA_BC::inflow;
            domain_bc_type[ori] = "Inflow";

            std::vector<Real> v;
            pp.getarr("velocity", v, 0, AMREX_SPACEDIM);
            m_bc_extdir_vals[BCVars::xvel_bc][ori] = v[0];
            m_bc_extdir_vals[BCVars::yvel_bc][ori] = v[1];
            m_bc_extdir_vals[BCVars::zvel_bc][ori] = v[2];

            Real scalar_in = 0.;
            if (pp.query("scalar", scalar_in))
            m_bc_extdir_vals[BCVars::Scalar_bc_comp][ori] = scalar_in;
        }
        else if (bc_type == "noslipwall")
        {
            phys_bc_type[ori] = REMORA_BC::no_slip_wall;
            domain_bc_type[ori] = "NoSlipWall";

            std::vector<Real> v;

            // The values of m_bc_extdir_vals default to 0.
            // But if we find "velocity" in the inputs file, use those values instead.
            if (pp.queryarr("velocity", v, 0, AMREX_SPACEDIM))
            {
                v[ori.coordDir()] = 0.0;
                m_bc_extdir_vals[BCVars::xvel_bc][ori] = v[0];
                m_bc_extdir_vals[BCVars::yvel_bc][ori] = v[1];
                m_bc_extdir_vals[BCVars::zvel_bc][ori] = v[2];
            }
        }
        else if (bc_type == "slipwall")
        {
            phys_bc_type[ori] = REMORA_BC::slip_wall;
            domain_bc_type[ori] = "SlipWall";
        }
        else
        {
            phys_bc_type[ori] = REMORA_BC::undefined;
        }

        if (geom[0].isPeriodic(ori.coordDir()))
        {
            domain_bc_type[ori] = "Periodic";
            if (phys_bc_type[ori] == REMORA_BC::undefined)
            {
                phys_bc_type[ori] = REMORA_BC::periodic;
            } else {
                amrex::Abort("Wrong BC type for periodic boundary");
            }
        }

        if (phys_bc_type[ori] == REMORA_BC::undefined)
        {
             amrex::Print() << "BC Type specified for face " << bcid << " is " << bc_type_in << std::endl;
             amrex::Abort("This BC type is unknown");
        }

        if ((bcid == "xlo" || bcid == "xhi" ||
             bcid == "ylo" || bcid == "yhi") &&
            solverChoice.ic_bc_type == IC_BC_Type::Real &&
            phys_bc_type[ori] != REMORA_BC::outflow)
        {
            amrex::Abort("BC type must be outflow in x and y when reading BCs from file");
        }
    };

    f("xlo", Orientation(Direction::x,Orientation::low));
    f("xhi", Orientation(Direction::x,Orientation::high));
    f("ylo", Orientation(Direction::y,Orientation::low));
    f("yhi", Orientation(Direction::y,Orientation::high));
    f("zlo", Orientation(Direction::z,Orientation::low));
    f("zhi", Orientation(Direction::z,Orientation::high));

    // *****************************************************************************
    //
    // Here we translate the physical boundary conditions -- one type per face --
    //     into logical boundary conditions for each velocity component
    //
    // *****************************************************************************
    {
        domain_bcs_type.resize(AMREX_SPACEDIM+NCONS);
        domain_bcs_type_d.resize(AMREX_SPACEDIM+NCONS);

        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            auto const bct = phys_bc_type[ori];
            if ( bct == REMORA_BC::symmetry )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, REMORABCType::reflect_even);
                    domain_bcs_type[BCVars::xvel_bc+dir].setLo(dir, REMORABCType::reflect_odd);
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, REMORABCType::reflect_even);
                    domain_bcs_type[BCVars::xvel_bc+dir].setHi(dir, REMORABCType::reflect_odd);
                }
            }
            else if (bct == REMORA_BC::outflow)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, REMORABCType::foextrap);
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, REMORABCType::foextrap);
                }
            }
            else if (bct == REMORA_BC::inflow)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++) {
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, REMORABCType::ext_dir);
                    }
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++) {
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, REMORABCType::ext_dir);
                    }
                }
            }
            else if (bct == REMORA_BC::no_slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, REMORABCType::ext_dir);
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, REMORABCType::ext_dir);
                }
            }
            else if (bct == REMORA_BC::slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, REMORABCType::foextrap);
                    // Only normal direction has ext_dir
                    domain_bcs_type[BCVars::xvel_bc+dir].setLo(dir, REMORABCType::ext_dir);

                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, REMORABCType::foextrap);
                    // Only normal direction has ext_dir
                    domain_bcs_type[BCVars::xvel_bc+dir].setHi(dir, REMORABCType::ext_dir);
                }
            }
            else if (bct == REMORA_BC::periodic)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, REMORABCType::int_dir);
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, REMORABCType::int_dir);
                }
            }
        }
    }

    // *****************************************************************************
    //
    // Here we translate the physical boundary conditions -- one type per face --
    //     into logical boundary conditions for each cell-centered variable
    //
    // *****************************************************************************
    {
        for (OrientationIter oit; oit; ++oit) {
            Orientation ori = oit();
            int dir = ori.coordDir();
            Orientation::Side side = ori.faceDir();
            auto const bct = phys_bc_type[ori];
            if ( bct == REMORA_BC::symmetry )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, REMORABCType::reflect_even);
                } else {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, REMORABCType::reflect_even);
                }
            }
            else if ( bct == REMORA_BC::outflow )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, REMORABCType::foextrap);
                } else {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, REMORABCType::foextrap);
                }
            }
            else if ( bct == REMORA_BC::no_slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, REMORABCType::foextrap);
                } else {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, REMORABCType::foextrap);
                }
            }
            else if (bct == REMORA_BC::slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, REMORABCType::foextrap);
                } else {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, REMORABCType::foextrap);
                }
            }
            else if (bct == REMORA_BC::inflow)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NCONS; i++) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, REMORABCType::ext_dir);
                    }
                } else {
                    for (int i = 0; i < NCONS; i++) {
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, REMORABCType::ext_dir);
                    }
                }
            }
            else if (bct == REMORA_BC::periodic)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, REMORABCType::int_dir);
                } else {
                    for (int i = 0; i < NCONS; i++)
                       domain_bcs_type[BCVars::cons_bc+i].setHi(dir, REMORABCType::int_dir);
                }
            }
        }
    }

#ifdef AMREX_USE_GPU
    Gpu::htod_memcpy
        (domain_bcs_type_d.data(), domain_bcs_type.data(),
         sizeof(amrex::BCRec)*(NCONS+AMREX_SPACEDIM));
#else
    std::memcpy
        (domain_bcs_type_d.data(), domain_bcs_type.data(),
         sizeof(amrex::BCRec)*(NCONS+AMREX_SPACEDIM));
#endif
}

