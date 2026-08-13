#include "REMORA_NCFile.H"
#include "AMReX_FArrayBox.H"
#include "REMORA_DataStruct.H"

using namespace amrex;

#ifdef REMORA_USE_NETCDF
/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_temp_fab     container for temperature data
 * @param NC_salt_fab     container for salinity data
 * @param NC_u_fab        container for u velocity data
 * @param NC_v_fab        container for v velocity data
 */
void
read_data_from_netcdf (int /*lev*/,
                       const Box& domain,
                       const std::string& fname,
                       FArrayBox& NC_temp_fab, FArrayBox& NC_salt_fab,
                       FArrayBox& NC_xvel_fab, FArrayBox& NC_yvel_fab)
{
    amrex::Print() << "Loading initial solution data from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_temp_fab); NC_names.push_back("temp");     NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 0
    NC_fabs.push_back(&NC_salt_fab); NC_names.push_back("salt");     NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 1
    NC_fabs.push_back(&NC_xvel_fab); NC_names.push_back("u");        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 2
    NC_fabs.push_back(&NC_yvel_fab); NC_names.push_back("v");        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 3

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}

/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_temp_fab     container for temperature data
 * @param NC_salt_fab     container for salinity data
 * @param NC_u_fab        container for u velocity data
 * @param NC_v_fab        container for v velocity data
 * @param ngrow           number of grow cells to read
 */
void
read_data_full_domain_from_netcdf (int /*lev*/,
                       const Box& domain,
                       const std::string& fname,
                       FArrayBox& NC_temp_fab, FArrayBox& NC_salt_fab,
                       FArrayBox& NC_xvel_fab, FArrayBox& NC_yvel_fab,
                       IntVect ngrow)
{
    amrex::Print() << "Loading initial solution data from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_temp_fab); NC_names.push_back("temp");     NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 0
    NC_fabs.push_back(&NC_salt_fab); NC_names.push_back("salt");     NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 1
    NC_fabs.push_back(&NC_xvel_fab); NC_names.push_back("u");        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 2
    NC_fabs.push_back(&NC_yvel_fab); NC_names.push_back("v");        NC_dim_types.push_back(NC_Data_Dims_Type::Time_BT_SN_WE); // 3

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs, false, 0, ngrow);
}

/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_zeta_fab     container for sea surface height data
 */
void
read_zeta_from_netcdf (int /*lev*/,
                      const Box& domain,
                      const std::string& fname,
                      FArrayBox& NC_zeta_fab)
{
    amrex::Print() << "Loading initial sea surface height from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_zeta_fab )   ; NC_names.push_back("zeta")    ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); // 0

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}

/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_zeta_fab     container for sea surface height data
 * @param ngrow           number of grow cells to read in, if not default
 */
void
read_zeta_full_domain_from_netcdf (int /*lev*/,
                      const Box& domain,
                      const std::string& fname,
                      FArrayBox& NC_zeta_fab,
                      IntVect ngrow)
{
    amrex::Print() << "Loading initial sea surface height from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_zeta_fab )   ; NC_names.push_back("zeta")    ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_SN_WE); // 0

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs, false, 0, ngrow);
}

/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_h_fab        container for bathymetry data
 */
void
read_bathymetry_from_netcdf (int /*lev*/,
                             const Box& domain,
                             const std::string& fname,
                             FArrayBox& NC_h_fab)
{
    amrex::Print() << "Loading initial bathymetry from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_h_fab )   ; NC_names.push_back("h")     ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 0

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}

/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_pm_fab       container for pm data
 * @param NC_pn_fab       container for pn data
 * @param NC_xr_fab       container for x_rho data
 * @param NC_yr_fab       container for y_rho data
 * @param NC_xu_fab       container for x_u data
 * @param NC_yu_fab       container for y_u data
 * @param NC_xv_fab       container for x_v data
 * @param NC_yv_fab       container for y_v data
 * @param NC_xp_fab       container for x_p data
 * @param NC_yp_fab       container for y_p data
 */
void
read_grid_vars_from_netcdf (int /*lev*/,
                             const Box& domain,
                             const std::string& fname,
                             FArrayBox& NC_pm_fab, FArrayBox& NC_pn_fab,
                             FArrayBox& NC_xr_fab, FArrayBox& NC_yr_fab,
                             FArrayBox& NC_xu_fab, FArrayBox& NC_yu_fab,
                             FArrayBox& NC_xv_fab, FArrayBox& NC_yv_fab,
                             FArrayBox& NC_xp_fab, FArrayBox& NC_yp_fab)
{
    amrex::Print() << "Loading grid variables from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_pm_fab)   ; NC_names.push_back("pm")    ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 1
    NC_fabs.push_back(&NC_pn_fab)   ; NC_names.push_back("pn")    ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 2
    NC_fabs.push_back(&NC_xr_fab)   ; NC_names.push_back("x_rho") ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 3
    NC_fabs.push_back(&NC_yr_fab)   ; NC_names.push_back("y_rho") ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 4
    NC_fabs.push_back(&NC_xu_fab)   ; NC_names.push_back("x_u")   ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 5
    NC_fabs.push_back(&NC_yu_fab)   ; NC_names.push_back("y_u")   ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 6
    NC_fabs.push_back(&NC_xv_fab)   ; NC_names.push_back("x_v")   ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 7
    NC_fabs.push_back(&NC_yv_fab)   ; NC_names.push_back("y_v")   ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 8
    NC_fabs.push_back(&NC_xp_fab)   ; NC_names.push_back("x_psi") ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 9
    NC_fabs.push_back(&NC_yp_fab)   ; NC_names.push_back("y_psi") ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 10

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}

/**
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_lonp_fab     container for lon_psi data
 * @param NC_latp_fab     container for lat_psi data
 * @returns whether the file carried spherical psi coordinates and they were read
 *
 * ROMS grid files for spherical grids carry lon_psi/lat_psi alongside
 * x_psi/y_psi, and those degree-valued corner arrays are what a coupled driver
 * needs when the two components do not share a projected Cartesian frame
 * (SCRIP/COAWST selects the same way; see Lib/SCRIP_COAWST/read_roms.f).
 * Idealized grid files carry no spherical coordinates, so this read is optional
 * and reports whether it happened rather than aborting.
 */
bool
read_spherical_grid_vars_from_netcdf (int /*lev*/,
                                      const Box& domain,
                                      const std::string& fname,
                                      FArrayBox& NC_lonp_fab, FArrayBox& NC_latp_fab)
{
    Vector<std::string> NC_names;
    NC_names.push_back("lon_psi");
    NC_names.push_back("lat_psi");

    if (!QueryNetCDFHasVars(fname, NC_names)) {
        amrex::Print() << "Grid file " << fname << " has no lon_psi/lat_psi; "
                       << "spherical psi coordinates unavailable" << std::endl;
        return false;
    }

    amrex::Print() << "Loading spherical psi coordinates from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_lonp_fab) ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE);
    NC_fabs.push_back(&NC_latp_fab) ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE);

    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
    return true;
}

/**
 * @param domain          simulation domain at nc_hires_grid_level
 * @param fname           file name to read from
 * @param NC_h_fab        container for bathymetry data
 * @param ngrow           number of grow cells to read in, if not default
 */
void
read_bathymetry_full_domain_from_netcdf (const Box& domain,
                                         const std::string& fname,
                                         FArrayBox& NC_h_fab, IntVect ngrow)
{
    amrex::Print() << "Loading high resolution bathymetry from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_h_fab )   ; NC_names.push_back("h")     ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 0

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs, false, 0, ngrow);
}

/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_pm_fab       container for pm data
 * @param NC_pn_fab       container for pn data
 * @param ngrow           number of grow cells to read in, if not default
 */
void
read_grid_vars_full_domain_from_netcdf (const Box& domain,
                             const std::string& fname,
                             FArrayBox& NC_pm_fab, FArrayBox& NC_pn_fab,
                             IntVect ngrow)
{
    amrex::Print() << "Loading high resolution grid variables from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_pm_fab)   ; NC_names.push_back("pm")    ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 0
    NC_fabs.push_back(&NC_pn_fab)   ; NC_names.push_back("pn")    ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 1
    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs, false, 0, ngrow);
}

/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_fcor_fab     container for Coriolis parameter data
 */
void
read_coriolis_from_netcdf (int /*lev*/,
                           const Box& domain,
                           const std::string& fname,
                           FArrayBox& NC_fcor_fab)
{
    amrex::Print() << "Loading initial coriolis from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_fcor_fab )   ; NC_names.push_back("f")    ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 0

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}

/**
 * @param lev             level of data to read
 * @param domain          simulation domain
 * @param fname           file name to read from
 * @param NC_mskr_fab     container for rho-point land/sea mask data
 * @param NC_msku_fab     container for u-point land/sea mask data
 * @param NC_mskv_fab     container for v-point land/sea mask data
 */
void
read_masks_from_netcdf (int /*lev*/,
                        const Box& domain,
                        const std::string& fname,
                        FArrayBox& NC_mskr_fab,
                        FArrayBox& NC_msku_fab,
                        FArrayBox& NC_mskv_fab)
{
    amrex::Print() << "Loading masks from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    NC_fabs.push_back(&NC_mskr_fab )   ; NC_names.push_back("mask_rho")  ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 0
    NC_fabs.push_back(&NC_msku_fab )   ; NC_names.push_back("mask_u")    ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 1
    NC_fabs.push_back(&NC_mskv_fab )   ; NC_names.push_back("mask_v")    ; NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE); // 2

    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}

/**
 * @param lev                 level of data to read
 * @param domain              simulation domain
 * @param fname               file name to read from
 * @param do_m2_clim_nudg     whether to do 2d momentum climatology nudging
 * @param do_m3_clim_nudg     whether to do 3d momentum climatology nudging
 * @param do_temp_clim_nudg   whether to do temperature climatology nudging
 * @param do_salt_clim_nudg   whether to do salinity climatology nudging
 * @param NC_M2NC_fab         container for 2d momentum climatology data
 * @param NC_M3NC_fab         container for 3d momentum climatology data
 * @param NC_TempNC_fab       container for temperature climatology data
 * @param NC_SaltNC_fab       container for salinity climatology data
 */
void
read_clim_nudg_coeff_from_netcdf (int /*lev*/,
                        const Box& domain,
                        const std::string& fname,
                        bool do_m2_clim_nudg,
                        bool do_m3_clim_nudg,
                        bool do_temp_clim_nudg,
                        bool do_salt_clim_nudg,
                        FArrayBox& NC_M2NC_fab,
                        FArrayBox& NC_M3NC_fab,
                        FArrayBox& NC_TempNC_fab,
                        FArrayBox& NC_SaltNC_fab)
{
    amrex::Print() << "Loading nudging coefficients from NetCDF file " << fname << std::endl;

    Vector<FArrayBox*> NC_fabs;
    Vector<std::string> NC_names;
    Vector<enum NC_Data_Dims_Type> NC_dim_types;

    if (do_m3_clim_nudg) {
        NC_fabs.push_back(&NC_M3NC_fab ); NC_names.push_back("M3_NudgeCoef"); NC_dim_types.push_back(NC_Data_Dims_Type::BT_SN_WE);
    }
    if (do_m2_clim_nudg) {
        NC_fabs.push_back(&NC_M2NC_fab ); NC_names.push_back("M2_NudgeCoef"); NC_dim_types.push_back(NC_Data_Dims_Type::SN_WE);
    }
    if (do_temp_clim_nudg) {
        NC_fabs.push_back(&NC_TempNC_fab ); NC_names.push_back("temp_NudgeCoef"); NC_dim_types.push_back(NC_Data_Dims_Type::BT_SN_WE);
    }
    if (do_salt_clim_nudg) {
        NC_fabs.push_back(&NC_SaltNC_fab ); NC_names.push_back("salt_NudgeCoef"); NC_dim_types.push_back(NC_Data_Dims_Type::BT_SN_WE);
    }
    // Read the netcdf file and fill these FABs
    BuildFABsFromNetCDFFile<FArrayBox,Real>(domain, fname, NC_names, NC_dim_types, NC_fabs);
}

/**
 * @param lev            level of data to read
 * @param fnames         file name(s) to read from
 * @param field_name     field name to read
 * @param vec_dat        vector to fill data
 */
//template <typename DType>
void read_vec_from_netcdf (int /*lev*/, const amrex::Vector<std::string>& fnames, const std::string& field_name, amrex::Vector<int>& vec_dat)
{
    AMREX_ALWAYS_ASSERT(!fnames.empty());
    const std::string& fname = fnames[0];

    amrex::Print() << "Reading " << field_name << " from NetCDF file" << std::endl;

    // get x-positions and put in array
    using ARRAY = NDArray<int>;
    amrex::Vector<ARRAY> array_dat(1);
    ReadNetCDFFile(fname, {field_name}, array_dat); // filled only on proc 0
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        int n = array_dat[0].get_vshape()[0];
        for (int i(0); i < n; i++)
        {
            vec_dat.push_back((*(array_dat[0].get_data() + i)));
        }
    }
    int nvals = vec_dat.size();
    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&nvals,1,ioproc);
    if (!(amrex::ParallelDescriptor::IOProcessor())) {
        vec_dat.resize(nvals);
    }
    amrex::ParallelDescriptor::Bcast(vec_dat.data(), vec_dat.size(), ioproc);
}

#endif // ROMSX_USE_NETCDF
