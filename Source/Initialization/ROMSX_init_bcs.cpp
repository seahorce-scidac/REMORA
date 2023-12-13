#include <AMReX_Vector.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ParmParse.H>

#include <ROMSX.H>

using namespace amrex;

void ROMSX::init_bcs ()
{
    auto f = [this] (std::string const& bcid, Orientation ori)
    {
        // These are simply defaults for Dirichlet faces -- they should be over-written below
        m_bc_extdir_vals[BCVars::Temp_bc_comp][ori] = -1.0; // It is important to set this negative
                                               // because the sign is tested on below
        m_bc_extdir_vals[BCVars::Scalar_bc_comp][ori] = 0.0;

        m_bc_extdir_vals[BCVars::xvel_bc][ori] = 0.0; // default
        m_bc_extdir_vals[BCVars::yvel_bc][ori] = 0.0;
        m_bc_extdir_vals[BCVars::zvel_bc][ori] = 0.0;

        ParmParse pp(bcid);
        std::string bc_type_in = "null";
        pp.query("type", bc_type_in);
        //if (pp.query("type", bc_type_in))
        //   amrex::Print() << "INPUT BC TYPE " << bcid << " " << bc_type_in << std::endl;
        std::string bc_type = amrex::toLower(bc_type_in);

        if (bc_type == "symmetry")
        {
            // amrex::Print() << bcid << " set to symmetry.\n";

            phys_bc_type[ori] = ROMSX_BC::symmetry;
            domain_bc_type[ori] = "Symmetry";
        }
        else if (bc_type == "outflow")
        {
            // amrex::Print() << bcid << " set to outflow.\n";

            phys_bc_type[ori] = ROMSX_BC::outflow;
            domain_bc_type[ori] = "Outflow";
        }
        else if (bc_type == "inflow")
        {
            // amrex::Print() << bcid << " set to inflow.\n";

            phys_bc_type[ori] = ROMSX_BC::inflow;
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
            // amrex::Print() << bcid <<" set to no-slip wall.\n";

            phys_bc_type[ori] = ROMSX_BC::no_slip_wall;
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
            // amrex::Print() << bcid <<" set to slip wall.\n";

            phys_bc_type[ori] = ROMSX_BC::slip_wall;
            domain_bc_type[ori] = "SlipWall";
        }
        else
        {
            phys_bc_type[ori] = ROMSX_BC::undefined;
        }

        if (geom[0].isPeriodic(ori.coordDir())) {
            domain_bc_type[ori] = "Periodic";
            if (phys_bc_type[ori] == ROMSX_BC::undefined)
            {
                phys_bc_type[ori] = ROMSX_BC::periodic;
            } else {
                amrex::Abort("Wrong BC type for periodic boundary");
            }
        }

        if (phys_bc_type[ori] == ROMSX_BC::undefined)
        {
             amrex::Print() << "BC Type specified for face " << bcid << " is " << bc_type_in << std::endl;
             amrex::Abort("This BC type is unknown");
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
            if ( bct == ROMSX_BC::symmetry )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ROMSXBCType::reflect_even);
                    domain_bcs_type[BCVars::xvel_bc+dir].setLo(dir, ROMSXBCType::reflect_odd);
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ROMSXBCType::reflect_even);
                    domain_bcs_type[BCVars::xvel_bc+dir].setHi(dir, ROMSXBCType::reflect_odd);
                }
            }
            else if (bct == ROMSX_BC::outflow)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ROMSXBCType::foextrap);
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ROMSXBCType::foextrap);
                }
            }
            else if (bct == ROMSX_BC::inflow)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++) {
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ROMSXBCType::ext_dir);
                    }
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++) {
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ROMSXBCType::ext_dir);
                    }
                }
            }
            else if (bct == ROMSX_BC::no_slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ROMSXBCType::ext_dir);
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ROMSXBCType::ext_dir);
                }
            }
            else if (bct == ROMSX_BC::slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ROMSXBCType::foextrap);
                    // Only normal direction has ext_dir
                    domain_bcs_type[BCVars::xvel_bc+dir].setLo(dir, ROMSXBCType::ext_dir);

                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ROMSXBCType::foextrap);
                    // Only normal direction has ext_dir
                    domain_bcs_type[BCVars::xvel_bc+dir].setHi(dir, ROMSXBCType::ext_dir);
                }
            }
            else if (bct == ROMSX_BC::periodic)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setLo(dir, ROMSXBCType::int_dir);
                } else {
                    for (int i = 0; i < AMREX_SPACEDIM; i++)
                        domain_bcs_type[BCVars::xvel_bc+i].setHi(dir, ROMSXBCType::int_dir);
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
            if ( bct == ROMSX_BC::symmetry )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ROMSXBCType::reflect_even);
                } else {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ROMSXBCType::reflect_even);
                }
            }
            else if ( bct == ROMSX_BC::outflow )
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ROMSXBCType::foextrap);
                } else {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ROMSXBCType::foextrap);
                }
            }
            else if ( bct == ROMSX_BC::no_slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ROMSXBCType::foextrap);
                    if (m_bc_extdir_vals[BCVars::Temp_bc_comp][ori] > 0.)
                        domain_bcs_type[BCVars::Temp_bc_comp].setLo(dir, ROMSXBCType::ext_dir);
                } else {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ROMSXBCType::foextrap);
                    if (m_bc_extdir_vals[BCVars::Temp_bc_comp][ori] > 0.)
                        domain_bcs_type[BCVars::Temp_bc_comp].setHi(dir, ROMSXBCType::ext_dir);
                }
            }
            else if (bct == ROMSX_BC::slip_wall)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ROMSXBCType::foextrap);
                    if (m_bc_extdir_vals[BCVars::Temp_bc_comp][ori] > 0.)
                        domain_bcs_type[BCVars::Temp_bc_comp].setLo(dir, ROMSXBCType::ext_dir);
                } else {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ROMSXBCType::foextrap);
                    if (m_bc_extdir_vals[BCVars::Temp_bc_comp][ori] > 0.)
                        domain_bcs_type[BCVars::Temp_bc_comp].setHi(dir, ROMSXBCType::ext_dir);
                }
            }
            else if (bct == ROMSX_BC::inflow)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NCONS; i++) {
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ROMSXBCType::ext_dir);
                    }
                } else {
                    for (int i = 0; i < NCONS; i++) {
                        domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ROMSXBCType::ext_dir);
                    }
                }
            }
            else if (bct == ROMSX_BC::periodic)
            {
                if (side == Orientation::low) {
                    for (int i = 0; i < NCONS; i++)
                        domain_bcs_type[BCVars::cons_bc+i].setLo(dir, ROMSXBCType::int_dir);
                } else {
                    for (int i = 0; i < NCONS; i++)
                       domain_bcs_type[BCVars::cons_bc+i].setHi(dir, ROMSXBCType::int_dir);
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

