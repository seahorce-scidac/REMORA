#include "REMORA_NCFile.H"
#include "AMReX_FArrayBox.H"
#include "REMORA_DataStruct.H"
#include "REMORA_IndexDefines.H"

using namespace amrex;

#ifdef REMORA_USE_NETCDF

namespace REMORABdyTypes {
    enum {
        x_lo,
        x_hi,
        y_lo,
        y_hi
    };
}

AMREX_FORCE_INLINE
std::time_t
getEpochTime (const std::string& dateTime, const std::string& dateTimeFormat)
{
    // Create a stream which we will use to parse the string,
    // which we provide to constructor of stream to fill the buffer.
    std::istringstream ss{ dateTime };

    // Create a tm object to store the parsed date and time.
    std::tm tmTime;
    memset(&tmTime, 0, sizeof(tmTime));

    // Now we read from buffer using get_time manipulator
    // and formatting the input appropriately.
    strptime(dateTime.c_str(), dateTimeFormat.c_str(), &tmTime);

    // Convert the tm structure to time_t value and return.
    // Here we use timegm since the output should be relative to UTC.
    auto epoch = timegm(&tmTime);
    // Print() << "Time Stamp: "<< std::put_time(&tmTime, "%c")
    //         << " , Epoch: " << epoch << std::endl;

    return epoch;
}

Real
read_bdry_from_netcdf (const Box& domain, const std::string& nc_bdry_file,
                       Vector<Vector<FArrayBox>>& bdy_data_xlo,
                       Vector<Vector<FArrayBox>>& bdy_data_xhi,
                       Vector<Vector<FArrayBox>>& bdy_data_ylo,
                       Vector<Vector<FArrayBox>>& bdy_data_yhi,
                       int& width, Real& start_bdy_time)
{
    amrex::Print() << "Loading boundary data from NetCDF file " << nc_bdry_file << std::endl;

    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank

    const auto& lo = domain.loVect();
    const auto& hi = domain.hiVect();

    // *******************************************************************************

    int ntimes;
    Real timeInterval;
    const std::string dateTimeFormat ="%Y-%m-%d_%H:%M:%S";

    Real ocean_times[31];
    // Check units of time stamps. Should be days.
    std::string unit_str;
    unit_str = ReadNetCDFVarAttrStr(nc_bdry_file, "ocean_time", "units"); // works on proc 0
    if (ParallelDescriptor::IOProcessor())
    {
        if (unit_str.find("days") == std::string::npos) {
            amrex::Print() << "Units of ocean_time given as: " << unit_str << std::endl;
            amrex::Abort("Units must be in days.");
        }
    }
    // Read the time stamps
    using RARRAY = NDArray<Real>;
    amrex::Vector<RARRAY> array_ts(1);
    ReadNetCDFFile(nc_bdry_file, {"ocean_time"}, array_ts); // filled only on proc 0
    if (ParallelDescriptor::IOProcessor())
    {
        ntimes = array_ts[0].get_vshape()[0];

        // amrex::Print() << " NTIMES " << ntimes << std::endl;
        for (int nt(0); nt < ntimes; nt++)
        {
            // Convert ocean time from days to seconds
            ocean_times[nt] = (*(array_ts[0].get_data() + nt)) * 60._rt * 60._rt * 24._rt;
            // amrex::Print() << "TIMES " << ocean_times[nt] << std::endl;
        }

        start_bdy_time = ocean_times[0];
        timeInterval = ocean_times[1] - ocean_times[0];

        for (int nt(1); nt < ntimes; nt++)
        {
            AMREX_ALWAYS_ASSERT(ocean_times[nt] - ocean_times[nt-1] == timeInterval);
        }
    }

    ParallelDescriptor::Bcast(&start_bdy_time,1,ioproc);
    ParallelDescriptor::Bcast(&ntimes,1,ioproc);
    ParallelDescriptor::Bcast(&timeInterval,1,ioproc);

    // Even though we may not read in all the variables, we need to make the arrays big enough for them (for now)
    int nvars = BdyVars::NumTypes*4; // 4 = xlo, xhi, ylo, yhi

    // Our outermost loop is time
    bdy_data_xlo.resize(ntimes);
    bdy_data_xhi.resize(ntimes);
    bdy_data_ylo.resize(ntimes);
    bdy_data_yhi.resize(ntimes);

    amrex::IntVect plo(lo);
    amrex::IntVect phi(hi);

    // ******************************************************************
    // Read the netcdf file and fill these FABs
    // NOTE: the order and number of these must match the BdyVars enum!
    // BdyVars:  U, V, R, T, QV, MU, PC
    // ******************************************************************
    Vector<std::string> nc_var_names;
    Vector<std::string> nc_var_prefix = {"u","v","temp","salt","ubar","vbar","zeta"};

    for (int ip = 0; ip < nc_var_prefix.size(); ++ip)
    {
       nc_var_names.push_back(nc_var_prefix[ip] + "_west");
       nc_var_names.push_back(nc_var_prefix[ip] + "_east");
       nc_var_names.push_back(nc_var_prefix[ip] + "_south");
       nc_var_names.push_back(nc_var_prefix[ip] + "_north");
    }

    using RARRAY = NDArray<Real>;
    amrex::Vector<RARRAY> arrays(nc_var_names.size());

    // The width of the boundary region we need to read is 1
    width = 1;

    ReadNetCDFFile(nc_bdry_file, nc_var_names, arrays); // does work on proc 0 only
    if (ParallelDescriptor::IOProcessor())
    {
        // Assert that the data has the same number of time snapshots
        int itimes = static_cast<int>(arrays[0].get_vshape()[0]);
        AMREX_ALWAYS_ASSERT(itimes == ntimes);

        // amrex::Print() << "VSHAPE " << arrays[0].get_vshape()[0] <<  " "
        //                             << arrays[0].get_vshape()[1] <<  " "
        //                             << arrays[0].get_vshape()[2] <<  " "
        //                             << arrays[0].get_vshape()[3] << std::endl;

    }
    ParallelDescriptor::Bcast(&width,1,ioproc);

    // This loops over every variable on every face, so nvars should be 4 * number of "ivartype" below
    for (int iv = 0; iv < nvars; iv++)
    {
        // amrex::Print() << "Building FAB for the NetCDF variable : " << nc_var_names[iv] << std::endl;

        int bdyVarType;

        std::string first1 = nc_var_names[iv].substr(0,1);
        std::string first4 = nc_var_names[iv].substr(0,4);

        if        (first4 == "salt") {
            bdyVarType = BdyVars::s;
        } else if (first4 == "ubar") {
            bdyVarType = BdyVars::ubar;
        } else if (first4 == "vbar") {
            bdyVarType = BdyVars::vbar;
        } else if (first4 == "zeta") {
            bdyVarType = BdyVars::zeta;
        } else if (first1 == "u") {
            bdyVarType = BdyVars::u;
        } else if (first1 == "v") {
            bdyVarType = BdyVars::v;
        } else if (first4 == "temp") {
            bdyVarType = BdyVars::t;
        } else {
            amrex::Print() << "Trying to read " << first1 << " or " << first4 << std::endl;
            amrex::Abort("dont know this variable");
        }

        std::string  last4 = nc_var_names[iv].substr(nc_var_names[iv].size()-4, 4);
        std::string  last5 = nc_var_names[iv].substr(nc_var_names[iv].size()-5, 5);
        int bdyType;

        if        (last4 == "west") {
            bdyType = REMORABdyTypes::x_lo;
        } else if (last4 == "east") {
            bdyType = REMORABdyTypes::x_hi;
        } else if (last5 == "south") {
            bdyType = REMORABdyTypes::y_lo;
        } else if (last5 == "north") {
            bdyType = REMORABdyTypes::y_hi;
        }

        plo[0] = lo[0]; plo[1] = lo[1]; plo[2] = lo[2];
        phi[0] = hi[0]; phi[1] = hi[1]; phi[2] = hi[2];
        const Box pbx(plo, phi);

        Arena* Arena_Used = The_Arena();
#ifdef AMREX_USE_GPU
        Arena_Used = The_Pinned_Arena();
#endif

        if (bdyType == REMORABdyTypes::x_lo) {

                // *******************************************************************************
                // xlo bdy
                // *******************************************************************************
                Box bx_for_u(IntVect(lo[0], lo[1]-1, lo[2]), IntVect(lo[0], hi[1]+1, hi[2]), IntVect(1,0,0));
                // amrex::Print() << "XLO:BX FOR U " << bx_for_u << std::endl;

                Box bx_for_v(IntVect(lo[0]-1, lo[1], lo[2]), IntVect(lo[0]-1, hi[1]+1, hi[2]), IntVect(0,1,0));
                // amrex::Print() << "XLO:BX FOR V " << bx_for_v << std::endl;

                Box bx_for_t(IntVect(lo[0]-1, lo[1]-1, lo[2]), IntVect(lo[0]-1, hi[1]+1, hi[2]));
                // amrex::Print() << "XLO:BX FOR T " << bx_for_t << std::endl;

                if        (bdyVarType == BdyVars::u) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xlo[nt].push_back(FArrayBox(bx_for_u, 1, Arena_Used)); // u
                    }
                } else if (bdyVarType == BdyVars::v) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xlo[nt].push_back(FArrayBox(bx_for_v , 1, Arena_Used)); // v
                    }
                } else if (bdyVarType == BdyVars::t) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xlo[nt].push_back(FArrayBox(bx_for_t, 1, Arena_Used)); // temp
                    }
                } else if (bdyVarType == BdyVars::s) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xlo[nt].push_back(FArrayBox(bx_for_t, 1, Arena_Used)); // salt
                    }
                } else if (bdyVarType == BdyVars::ubar) {
                    Box xlo_ubar(makeSlab(bx_for_u,2,0));
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xlo[nt].push_back(FArrayBox(xlo_ubar, 1, Arena_Used)); // ubar
                    }
                } else if (bdyVarType == BdyVars::vbar) {
                    Box xlo_vbar(makeSlab(bx_for_v,2,0));
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xlo[nt].push_back(FArrayBox(xlo_vbar, 1, Arena_Used)); // ubar
                    }
                } else if (bdyVarType == BdyVars::zeta) {
                    Box xlo_zeta(makeSlab(bx_for_t,2,0));
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xlo[nt].push_back(FArrayBox(xlo_zeta, 1, Arena_Used)); // ubar
                    }
                }

            } else if (bdyType == REMORABdyTypes::x_hi) {

                // *******************************************************************************
                // xhi bdy
                // *******************************************************************************

                Box bx_for_u(IntVect(hi[0]+1, lo[1]-1, lo[2]), IntVect(hi[0]+1, hi[1]+1, hi[2]), IntVect(1,0,0));
                // amrex::Print() << "XHI:BX FOR U " << bx_for_u << std::endl;

                Box bx_for_v(IntVect(hi[0]+1, lo[1], lo[2]), IntVect(hi[0]+1, hi[1]+1, hi[2]), IntVect(0,1,0));
                // amrex::Print() << "XHI:BX FOR V " << bx_for_v << std::endl;

                Box bx_for_t(IntVect(hi[0]+1, lo[1]-1, lo[2]), IntVect(hi[0]+1, hi[1]+1, hi[2]));
                // amrex::Print() << "XHI:BX FOR T " << bx_for_t << std::endl;

                if        (bdyVarType == BdyVars::u) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xhi[nt].push_back(FArrayBox(bx_for_u, 1, Arena_Used)); // u
                    }
                } else if (bdyVarType == BdyVars::v) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xhi[nt].push_back(FArrayBox(bx_for_v , 1, Arena_Used)); // v
                    }
                } else if (bdyVarType == BdyVars::t) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xhi[nt].push_back(FArrayBox(bx_for_t, 1, Arena_Used)); // temp
                    }
                } else if (bdyVarType == BdyVars::s) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xhi[nt].push_back(FArrayBox(bx_for_t, 1, Arena_Used)); // salt
                    }
                } else if (bdyVarType == BdyVars::ubar) {
                    Box xhi_ubar(makeSlab(bx_for_u,2,0));
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xhi[nt].push_back(FArrayBox(xhi_ubar, 1, Arena_Used)); // ubar
                    }
                } else if (bdyVarType == BdyVars::vbar) {
                    Box xhi_vbar(makeSlab(bx_for_v,2,0));
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xhi[nt].push_back(FArrayBox(xhi_vbar, 1, Arena_Used)); // ubar
                    }
                } else if (bdyVarType == BdyVars::zeta) {
                    Box xhi_zeta(makeSlab(bx_for_t,2,0));
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xhi[nt].push_back(FArrayBox(xhi_zeta, 1, Arena_Used)); // ubar
                    }
                }

            } else if (bdyType == REMORABdyTypes::y_lo) {

                // *******************************************************************************
                // ylo bdy
                // *******************************************************************************

                Box bx_for_v(IntVect(lo[0]-1, lo[1], lo[2]), IntVect(hi[0]+1, lo[1], hi[2]), IntVect(0,1,0));
                // amrex::Print() << "YLO:BX FOR V " << bx_for_v << std::endl;

                Box bx_for_u(IntVect(lo[0], lo[1]-1, lo[2]), IntVect(hi[0]+1, lo[1]-1, hi[2]), IntVect(1,0,0));
                // amrex::Print() << "YLO:BX FOR U " << bx_for_u << std::endl;

                Box bx_for_t(IntVect(lo[0]-1, lo[1]-1, lo[2]), IntVect(hi[0]+1, lo[1]-1, hi[2]));
                // amrex::Print() << "YLO:BX FOR T " << bx_for_t << std::endl;

                if        (bdyVarType == BdyVars::u) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_ylo[nt].push_back(FArrayBox(bx_for_u , 1, Arena_Used)); // u
                    }
                } else if (bdyVarType == BdyVars::v) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_ylo[nt].push_back(FArrayBox(bx_for_v, 1, Arena_Used)); // v
                    }
                } else if (bdyVarType == BdyVars::t) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_ylo[nt].push_back(FArrayBox(bx_for_t, 1, Arena_Used)); // temp
                    }
                } else if (bdyVarType == BdyVars::s) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_ylo[nt].push_back(FArrayBox(bx_for_t, 1, Arena_Used)); // salt
                    }
                } else if (bdyVarType == BdyVars::ubar) {
                    Box ylo_ubar(makeSlab(bx_for_u,2,0));
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_ylo[nt].push_back(FArrayBox(ylo_ubar, 1, Arena_Used)); // ubar
                    }
                } else if (bdyVarType == BdyVars::vbar) {
                    Box ylo_vbar(makeSlab(bx_for_v,2,0));
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_ylo[nt].push_back(FArrayBox(ylo_vbar, 1, Arena_Used)); // ubar
                    }
                } else if (bdyVarType == BdyVars::zeta) {
                    Box ylo_zeta(makeSlab(bx_for_t,2,0));
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_ylo[nt].push_back(FArrayBox(ylo_zeta, 1, Arena_Used)); // ubar
                    }
                }

            } else if (bdyType == REMORABdyTypes::y_hi) {

                // *******************************************************************************
                // yhi bdy
                // *******************************************************************************

                Box bx_for_v(IntVect(lo[0]-1, hi[1]+1, lo[2]), IntVect(hi[0]+1, hi[1]+1, hi[2]), IntVect(0,1,0));
                // amrex::Print() << "YHI:BX FOR V " << bx_for_v << std::endl;

                Box bx_for_u(IntVect(lo[0], hi[1]+1, lo[2]), IntVect(hi[0]+1, hi[1]+1, hi[2]), IntVect(1,0,0));
                // amrex::Print() << "YHI:BX FOR U " << bx_for_u << std::endl;

                Box bx_for_t(IntVect(lo[0]-1, hi[1]+1, lo[2]), IntVect(hi[0]+1, hi[1]+1, hi[2]));
                // amrex::Print() << "YHI:BX FOR T " << bx_for_t << std::endl;

                if        (bdyVarType == BdyVars::u) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_yhi[nt].push_back(FArrayBox(bx_for_u , 1, Arena_Used)); // u
                    }
                } else if (bdyVarType == BdyVars::v) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_yhi[nt].push_back(FArrayBox(bx_for_v, 1, Arena_Used)); // v
                    }
                } else if (bdyVarType == BdyVars::t) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_yhi[nt].push_back(FArrayBox(bx_for_t, 1, Arena_Used)); // temp
                    }
                } else if (bdyVarType == BdyVars::s) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_yhi[nt].push_back(FArrayBox(bx_for_t, 1, Arena_Used)); // salt
                    }
                } else if (bdyVarType == BdyVars::ubar) {
                    Box yhi_ubar(makeSlab(bx_for_u,2,0));
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_yhi[nt].push_back(FArrayBox(yhi_ubar, 1, Arena_Used)); // ubar
                    }
                } else if (bdyVarType == BdyVars::vbar) {
                    Box yhi_vbar(makeSlab(bx_for_v,2,0));
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_yhi[nt].push_back(FArrayBox(yhi_vbar, 1, Arena_Used)); // ubar
                    }
                } else if (bdyVarType == BdyVars::zeta) {
                    Box yhi_zeta(makeSlab(bx_for_t,2,0));
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_yhi[nt].push_back(FArrayBox(yhi_zeta, 1, Arena_Used)); // ubar
                    }
                }
        }

        long n_plane;

        // Now fill the data
        if (ParallelDescriptor::IOProcessor())
        {
            // amrex::Print() << "SHAPE0 " << arrays[iv].get_vshape()[0] << std::endl;
            // amrex::Print() << "SHAPE1 " << arrays[iv].get_vshape()[1] << std::endl;
            // amrex::Print() << "SHAPE2 " << arrays[iv].get_vshape()[2] << std::endl;
            // amrex::Print() << "SHAPE3 " << arrays[iv].get_vshape()[3] << std::endl;

            int nx, ny, nz;

            Array4<Real> fab_arr;
            Box my_box;

            if (bdyType == REMORABdyTypes::x_lo) {
                my_box = bdy_data_xlo[0][bdyVarType].box();
            } else if (bdyType == REMORABdyTypes::x_hi) {
                my_box = bdy_data_xhi[0][bdyVarType].box();
            } else if (bdyType == REMORABdyTypes::y_lo) {
                my_box = bdy_data_ylo[0][bdyVarType].box();
            } else if (bdyType == REMORABdyTypes::y_hi) {
                my_box = bdy_data_yhi[0][bdyVarType].box();
            }

            {
                if ( (bdyType == REMORABdyTypes::x_lo) || (bdyType == REMORABdyTypes::x_hi) ) {

                    if (my_box.length()[2] == 1) {
                        nz = 1;
                        ny = arrays[iv].get_vshape()[1];
                    } else {
                        nz = arrays[iv].get_vshape()[1];
                        ny = arrays[iv].get_vshape()[2];
                    }
                    n_plane = ny * nz;

                    AMREX_ALWAYS_ASSERT(my_box.numPts() == n_plane);

                    int i    = my_box.smallEnd()[0];
                    int joff = my_box.smallEnd()[1];

                    for (int nt(0); nt < ntimes; ++nt)
                    {
                        if (bdyType == REMORABdyTypes::x_lo) {
                            fab_arr  = bdy_data_xlo[nt][bdyVarType].array();
                        } else if (bdyType == REMORABdyTypes::x_hi) {
                            fab_arr  = bdy_data_xhi[nt][bdyVarType].array();
                        }
                        int n_off = nt * n_plane;

                        for (int n(0); n < n_plane; ++n) {
                            int k = n / ny;
                            int j = n - (k * ny);
                            fab_arr(i, j+joff, k, 0) = static_cast<Real>(*(arrays[iv].get_data() + n + n_off));
                        }
                    }

                } else if ( (bdyType == REMORABdyTypes::y_lo) || (bdyType == REMORABdyTypes::y_hi) ) {

                    if (my_box.length()[2] == 1) {
                        nz = 1;
                        nx = arrays[iv].get_vshape()[1];
                    } else {
                        nz = arrays[iv].get_vshape()[1];
                        nx = arrays[iv].get_vshape()[2];
                    }
                    n_plane = nx * nz;

                    AMREX_ALWAYS_ASSERT(my_box.numPts() == n_plane);

                    int j    = my_box.smallEnd()[1];
                    int ioff = my_box.smallEnd()[0];

                    for (int nt(0); nt < ntimes; ++nt)
                    {
                        if (bdyType == REMORABdyTypes::y_lo) {
                            fab_arr  = bdy_data_ylo[nt][bdyVarType].array();
                        } else if (bdyType == REMORABdyTypes::y_hi) {
                            fab_arr  = bdy_data_yhi[nt][bdyVarType].array();
                        }
                        int n_off = nt * n_plane;

                        for (int n(0); n < n_plane; ++n) {
                            int k = n / nx;
                            int i = n - (k * nx);

                            fab_arr(i+ioff, j, k, 0) = static_cast<Real>(*(arrays[iv].get_data() + n + n_off));
                        }
                    }
                } // bdyType
            } // bdyVarType
        } // if ParalleDescriptor::IOProcessor()
    } // nc_var_names

    // We put a barrier here so the rest of the processors wait to do anything until they have the data
    amrex::ParallelDescriptor::Barrier();

    // When an FArrayBox is built, space is allocated on every rank.  However, we only
    //    filled the data in these FABs on the IOProcessor.  So here we broadcast
    //    the data to every rank.
    int n_per_time = nc_var_prefix.size();
    for (int nt = 0; nt < ntimes; nt++)
    {
        for (int i = 0; i < n_per_time; i++)
        {
            ParallelDescriptor::Bcast(bdy_data_xlo[nt][i].dataPtr(),bdy_data_xlo[nt][i].box().numPts(),ioproc);
            ParallelDescriptor::Bcast(bdy_data_xhi[nt][i].dataPtr(),bdy_data_xhi[nt][i].box().numPts(),ioproc);
            ParallelDescriptor::Bcast(bdy_data_ylo[nt][i].dataPtr(),bdy_data_ylo[nt][i].box().numPts(),ioproc);
            ParallelDescriptor::Bcast(bdy_data_yhi[nt][i].dataPtr(),bdy_data_yhi[nt][i].box().numPts(),ioproc);
        }
    }

    // Make sure all processors know how timeInterval
    ParallelDescriptor::Bcast(&timeInterval,1,ioproc);

    // Make sure all processors know how many times are stored
    ParallelDescriptor::Bcast(&ntimes,1,ioproc);

    // Return the number of seconds between the boundary plane data
    return timeInterval;
}
#endif // REMORA_USE_NETCDF
