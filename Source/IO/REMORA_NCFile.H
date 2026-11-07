#ifndef _REMORA_NCFILE_H_
#define _REMORA_NCFILE_H_

#include <initializer_list>
#include <sstream>
#include <string>
#include <ctime>
#include <memory>

#include "AMReX_FArrayBox.H"
#include "AMReX_IArrayBox.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParallelContext.H"
#include "AMReX_ParallelReduce.H"
#include "AMReX_ParmParse.H"
// this includes pnetcdf.h, and mpi.h
#include "REMORA_NCInterface.H"

using PlaneVector = amrex::Vector<amrex::FArrayBox>;

enum class NC_Data_Dims_Type {
    Time_BT_SN_WE,
    Time_SN_WE,
    BT_SN_WE,
    SN_WE,
    Time_BT_Riv,
    Time_Riv,
    BT_Riv,
};

/**
 * \brief NDArray is the datatype designed to hold any data, including scalars, multidimensional
 * arrays, that read from the NetCDF file.
 *
 * The data read from NetCDF file are stored in a continuous memory, and the data layout is described
 * by using a vector (shape). AMRex Box can be constructed using the data shape information, and MultiFab
 * data array can be setup using the data that stored in the NDArray.
 */
template <typename DataType>
struct NDArray
{
   using DType = typename std::remove_const<DataType>::type;

   /** \brief constructor */
   explicit NDArray (const std::string vname, const std::vector<MPI_Offset>& vshape)
                    : name{vname}, shape{vshape},
                      data{std::shared_ptr<DType[]>(new DType[this->ndim()],
                                                   std::default_delete<DType[]>())} {}

   /** \brief default constructor */
   NDArray () : name{"null"}, data{nullptr} {}

   /** \brief get the data pointer */
   decltype(auto) get_data () {
     return data.get();
   }

   /** \brief get the variable name */
   std::string get_vname () {
     return name;
   }

   /** \brief get the variable data shape */
   std::vector<MPI_Offset> get_vshape () {
     return shape;
   }

   /** \brief return the total number of data */
   size_t ndim () {
     size_t num = 1;
     int isize = static_cast<int>(shape.size());
     for (auto i=0; i<isize; ++i) num *= shape[i];
     return num;
   }

   /** \brief set the data shape information */
   void set_vshape (std::vector<MPI_Offset> vshape) {
     shape = vshape;
   }

 private:
   std::string name;
   std::vector<MPI_Offset> shape;
   std::shared_ptr<DType[]> data;
};

/** \brief Read in data from netcdf file and save to data arrays
 *
 * Only the I/O processor opens the file and fills \p arrays; the other ranks
 * return at once, and the caller broadcasts what it needs.
 *
 * @param[in   ] fname       file name to read from
 * @param[in   ] names       field names to read
 * @param[inout] arrays      vector of data arrays to fill
 * @param[in   ] one_time    whether to read in a single time step of data
 * @param[in   ] fill_time   if reading a single time, what time step to read from the file?
 */
template<typename DType>
void ReadNetCDFFile (const std::string& fname, amrex::Vector<std::string> names,
                     amrex::Vector<NDArray<DType> >& arrays, bool one_time=false, int fill_time=0)
{
    AMREX_ASSERT(arrays.size() == names.size());

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        auto ncf = ncutils::NCFile::open_local(fname);

        /*
        // get the dimension information
        int Time               = static_cast<int>(ncf.dim("Time").len());
        int DateStrLen         = static_cast<int>(ncf.dim("DateStrLen").len());
        int west_east          = static_cast<int>(ncf.dim("west_east").len());
        int south_north        = static_cast<int>(ncf.dim("south_north").len());
        int bottom_top         = static_cast<int>(ncf.dim("bottom_top").len());
        int bottom_top_stag    = static_cast<int>(ncf.dim("bottom_top_stag").len());
        int west_east_stag     = static_cast<int>(ncf.dim("west_east_stag").len());
        int south_north_stag   = static_cast<int>(ncf.dim("south_north_stag").len());
        */

//        amrex::Print() << "Reading the dimensions from the netcdf file " << "\n";

        for (auto n=0; n<arrays.size(); ++n) {
            // read the data from NetCDF file
            std::string vname_to_write = names[n];
            std::string vname_to_read  = names[n];

            if (!ncf.has_var(vname_to_read)) {
               amrex::Print() << "Variable name " << vname_to_read << " not found!" << std::endl;
            }

            // amrex::AllPrint() << "About to read " << vname_to_read
            //                   << " while filling the array for " << vname_to_write << std::endl;

            std::vector<MPI_Offset> shape = ncf.var(vname_to_read).shape();
            std::vector<MPI_Offset> start(shape.size(), 0);
            // If only reading one time step, adjust start and shape accordingly
            if (one_time) {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!shape.empty(),
                    "Cannot read one time step from a scalar NetCDF variable");
                start[0] = fill_time;
                shape[0] = 1;
            }

            arrays[n]                 = NDArray<DType>(vname_to_read,shape);
            DType* dataPtr            = arrays[n].get_data();

            // auto numPts               = arrays[n].ndim();
            // amrex::Print() << "NetCDF Variable name = " << vname_to_read << std::endl;
            // amrex::Print() << "numPts read from NetCDF file/var = " << numPts << std::endl;
            // amrex::Print()  << "Dims of the variable = (";
            // for (auto &dim:shape)
            //     amrex::Print() << dim << ", " ;
            // amrex::Print()  << ")" << std::endl;

            ncf.var(vname_to_read).get_all(dataPtr, start, shape);
        }
        ncf.close();
    }
}

/**
 * \brief Helper function for reading a single variable attribute
 */
std::string ReadNetCDFVarAttrStr (const std::string& fname,
                                  const std::string& var_name,
                                  const std::string& attr_name);

/**
 * \brief Helper function for testing for the presence of a single variable attribute
 */
bool QueryNetCDFVarAttrStr (const std::string& fname,
                            const std::string& var_name,
                            const std::string& attr_name);

/**
 * \brief Offset, in days, that puts a file's time variable onto the model clock.
 *
 * Reads the `units` attribute of \p time_name, which must count days, and its
 * `calendar` attribute if it has one. A `days since <date>` reference date is
 * converted to days since the date `remora.time_ref` names, so a file whose
 * epoch differs from the run's is still read at the right time; the offset is
 * zero when the two dates agree. A time variable with no reference date is
 * taken to be on the model clock already. A calendar that disagrees with
 * \p time_ref is an error, since the two would count days differently.
 *
 * Collective: the attributes are read on the I/O processor and the result is
 * broadcast. The mapping applied is printed, so a reader of the log can see it.
 *
 * @param[in] fname           NetCDF file
 * @param[in] time_name       name of the time variable
 * @param[in] time_ref        remora.time_ref
 * @param[in] units_required  abort when the time variable has no `units`
 *                            attribute; otherwise warn and assume days on the
 *                            model clock
 * @param[in] announce        print the mapping; false for a series that is set
 *                            up but will never be read
 */
double ReadNetCDFTimeOffsetDays (const std::string& fname,
                                 const std::string& time_name,
                                 double time_ref,
                                 bool units_required,
                                 bool announce);

/**
 * \brief Helper function for testing whether a file carries every named variable.
 *
 * Call this before reading optional variables: ReadNetCDFFile() dereferences
 * ncf.var(name) unconditionally, so a missing name is a hard failure.
 */
bool QueryNetCDFHasVars (const std::string& fname,
                         const amrex::Vector<std::string>& var_names);

/**
 * \brief Where fill_fab_from_arrays puts a variable's data, before the shift by
 * the domain's lower corner.
 *
 * @param[in ] var_name      variable name, which sets the staggering
 * @param[in ] NC_dim_type   dimension type of the variable in the file
 * @param[in ] shape         the variable's shape in the file (any time count)
 * @param[in ] custom_ngrow  grow cells of the FAB, as in BuildFABsFromNetCDFFile
 * @param[out] nsx, nsy, nsz extents of the file data in x, y and z
 * @returns the FAB's box. The data fill the nsx x nsy x nsz corner of it that
 *          starts at its lower corner.
 */
inline amrex::Box
nc_fab_box (const std::string& var_name, NC_Data_Dims_Type NC_dim_type,
            const std::vector<MPI_Offset>& shape, amrex::IntVect custom_ngrow,
            int& nsx, int& nsy, int& nsz)
{
    if (NC_dim_type == NC_Data_Dims_Type::SN_WE) {
        nsz = 1;
        nsy = shape[0];
        nsx = shape[1];
    } else if (NC_dim_type == NC_Data_Dims_Type::BT_SN_WE) {
        nsz = shape[0];
        nsy = shape[1];
        nsx = shape[2];
    } else if (NC_dim_type == NC_Data_Dims_Type::Time_SN_WE) {
        nsz = 1;
        nsy = shape[1];
        nsx = shape[2];
    } else if (NC_dim_type == NC_Data_Dims_Type::Time_BT_SN_WE) {
        nsz = shape[1];
        nsy = shape[2];
        nsx = shape[3];
    } else if (NC_dim_type == NC_Data_Dims_Type::Time_BT_Riv) {
        nsz = shape[1]; // number of z-coords
        nsy = 1; // dummy
        nsx = shape[2]; // number of rivers
    } else if (NC_dim_type == NC_Data_Dims_Type::Time_Riv) {
        nsz = 1; // number of z-coords
        nsy = 1; // dummy
        nsx = shape[1]; // number of rivers
    } else if (NC_dim_type == NC_Data_Dims_Type::BT_Riv) {
        nsz = shape[0]; // number of z-coords
        nsy = 1; // dummy
        nsx = shape[1]; // number of rivers
    } else {
        amrex::Abort("Dont know this NC_Data_Dims_Type");
    }

    amrex::Box my_box;

    if (var_name == "u" || var_name == "ubar" || var_name == "mask_u" || var_name == "x_u" ||
          var_name == "y_u" || var_name == "sustr")
    {
        my_box.setSmall(amrex::IntVect(0,-1,0));
        my_box.setBig(amrex::IntVect(nsx-1,nsy-2,nsz-1));
        my_box.setType(amrex::IndexType(amrex::IntVect(1,0,0)));
    }
    else if (var_name == "v" || var_name == "vbar" || var_name == "mask_v" || var_name == "x_v" ||
          var_name == "y_v" || var_name == "svstr")
    {
        my_box.setSmall(amrex::IntVect(-1,0,0));
        my_box.setBig(amrex::IntVect(nsx-2,nsy-1,nsz-1));
        my_box.setType(amrex::IndexType(amrex::IntVect(0,1,0)));
    }
    else if (var_name == "mask_psi" || var_name == "x_psi" || var_name == "y_psi" ||
          var_name == "lon_psi" || var_name == "lat_psi")
    {
        my_box.setSmall(amrex::IntVect(0,0,0));
        my_box.setBig(amrex::IntVect(nsx-1,nsy-1,nsz-1));
        my_box.setType(amrex::IndexType(amrex::IntVect(1,1,0)));
    // River data has no ghost cells
    } else if (NC_dim_type == NC_Data_Dims_Type::Time_BT_Riv ||
               NC_dim_type == NC_Data_Dims_Type::Time_Riv ||
               NC_dim_type == NC_Data_Dims_Type::BT_Riv) {
        my_box.setSmall(amrex::IntVect(0,0,0));
        my_box.setBig(amrex::IntVect(nsx-1,nsy-1,nsz-1));
    }
    else  // everything cell-centered -- these have one ghost cell
    {
        my_box.setSmall(amrex::IntVect(-1,-1,0));
        my_box.setBig(amrex::IntVect(nsx-2,nsy-2,nsz-1));
    }
    // Adjust box size if a non-standard number of grow cells
    my_box.grow(custom_ngrow - amrex::IntVect(1,1,0));
    return my_box;
}

/**
 * \brief Helper function for reading data from NetCDF file into a
 * provided FAB.
 *
 * @param iv Index for which variable we are going to fill
 * @param nc_arrays Arrays of data from NetCDF file
 * @param var_name Variable name
 * @param NC_dim_type Dimension type for the variable as stored in the NetCDF file
 * @param temp FAB where we store the variable data from the NetCDF Arrays
 * @param custom_ngrow Custom number of grow cells for the FAB
 */
template<class FAB,typename DType>
void
fill_fab_from_arrays (int iv,
                      amrex::Vector<NDArray<amrex::Real>>& nc_arrays,
                      const std::string& var_name,
                      NC_Data_Dims_Type& NC_dim_type,
                      FAB& temp, amrex::IntVect custom_ngrow)
{
    int nsx, nsy, nsz;
    const amrex::Box my_box = nc_fab_box(var_name, NC_dim_type, nc_arrays[iv].get_vshape(),
                                         custom_ngrow, nsx, nsy, nsz);

#ifdef AMREX_USE_GPU
    // Make sure temp lives on CPU since nc_arrays lives on CPU only
    temp.resize(my_box,1,amrex::The_Pinned_Arena());
#else
    temp.resize(my_box,1);
#endif
    amrex::Array4<DType> fab_arr = temp.array();

    int ioff = temp.box().smallEnd()[0];
    int joff = temp.box().smallEnd()[1];

    auto num_pts_in_box = my_box.numPts();
    auto num_pts_in_src = nsx*nsy*nsz;
    // amrex::Print() <<" NUM PTS IN BOX " << num_pts_in_box << std::endl;
    // amrex::Print() <<" NUM PTS IN SRC " << num_pts_in_src << std::endl;
    AMREX_ALWAYS_ASSERT(num_pts_in_box >= num_pts_in_src);

    for (int n(0); n < num_pts_in_src; ++n) {
        int nplane = nsx*nsy;
        int k  =  n / nplane;
        int j  = (n - k*nplane) / nsx + joff;
        int i  =  n - k*nplane - (j-joff) * nsx + ioff;
        fab_arr(i,j,k,0) = static_cast<DType>(*(nc_arrays[iv].get_data()+n));
    }
}

/**
 * \brief Function to read NetCDF variables and fill the corresponding Array4's
 *
 * @param fname Name of the NetCDF file to be read
 * @param nc_var_names Variable names in the NetCDF file
 * @param NC_dim_types NetCDF data dimension types
 * @param fab_vars Fab data we are to fill
 * @param one_time Whether to read in a single time step from file
 * @param fill_time What time step to read from file, if only reading one step
 * @param custom_ngrow Number of grow cells to read in the fab; default is based off of one cell in each of x and y directions
 */
template<class FAB,typename DType>
void
BuildFABsFromNetCDFFile (const amrex::Box& domain,
                         const std::string &fname,
                         amrex::Vector<std::string> nc_var_names,
                         amrex::Vector<enum NC_Data_Dims_Type> NC_dim_types,
                         amrex::Vector<FAB*> fab_vars,
                         bool one_time=false, int fill_time=0,
                         amrex::IntVect custom_ngrow=amrex::IntVect(1,1,0))
{
    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();  // I/O rank

    amrex::Vector<NDArray<amrex::Real>> nc_arrays(nc_var_names.size());


    ReadNetCDFFile(fname, nc_var_names, nc_arrays, one_time, fill_time); // this is filled only on proc 0

    for (int iv = 0; iv < nc_var_names.size(); iv++)
    {
        FAB tmp;
        if (amrex::ParallelDescriptor::IOProcessor()) {
            fill_fab_from_arrays<FAB,DType>(iv, nc_arrays, nc_var_names[iv], NC_dim_types[iv], tmp, custom_ngrow);
        }

        int ncomp = tmp.nComp();
        amrex::Box box = tmp.box();

        amrex::ParallelDescriptor::Bcast(&box,   1, ioproc);
        amrex::ParallelDescriptor::Bcast(&ncomp, 1, ioproc);

        if (!amrex::ParallelDescriptor::IOProcessor()) {
#ifdef AMREX_USE_GPU
            tmp.resize(box,ncomp,amrex::The_Pinned_Arena());
#else
            tmp.resize(box,ncomp);
#endif
        }

        amrex::ParallelDescriptor::Bcast(tmp.dataPtr(), tmp.size(), ioproc);

        // Shift box by the domain lower corner
        amrex::Box  fab_bx = tmp.box();
        amrex::Dim3 dom_lb = lbound(domain);
        fab_bx += amrex::IntVect(dom_lb.x,dom_lb.y,dom_lb.z);
        // fab_vars points to data on device
        fab_vars[iv]->resize(fab_bx,1);
#ifdef AMREX_USE_GPU
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         tmp.dataPtr(), tmp.dataPtr() + tmp.size(),
                         fab_vars[iv]->dataPtr());
#else
        // Provided by BaseFab inheritance through FArrayBox
        fab_vars[iv]->copy(tmp,tmp.box(),0,fab_bx,0,1);
#endif
    }
}

/**
 * \brief The cell-centered bounding box of every FAB, ghost cells included,
 * that this rank owns in each of \p mfs: the region BuildLocalFABsFromNetCDFFile
 * has to read for a copy into them. Empty on a rank that owns none. Null
 * entries are skipped.
 *
 * One box per rank keeps each variable to one collective read per rank. Its
 * cost is the bounding box: a rank whose FABs lie far apart reads, and holds,
 * everything between them. With remora.v >= 1 that cost is printed: the
 * largest ratio over the ranks of cells in the region to cells the rank owns.
 *
 * Collective over amrex::ParallelContext::CommunicatorSub() when remora.v >= 1.
 */
inline amrex::Box
nc_local_region (const amrex::Vector<const amrex::MultiFab*>& mfs)
{
    amrex::Box region;
    amrex::BoxList owned_boxes;
    for (const amrex::MultiFab* mf : mfs) {
        if (mf == nullptr) { continue; }
        for (amrex::MFIter mfi(*mf, false); mfi.isValid(); ++mfi) {
            const amrex::Box b = amrex::convert(mfi.fabbox(), amrex::IndexType::TheCellType());
            region = region.ok() ? region.minBox(b) : b;
            owned_boxes.push_back(b);
        }
    }

    int verbose = 0;
    amrex::ParmParse("remora").query("v", verbose);
    if (verbose >= 1) {
        // Cells this rank's FABs cover, each counted once although several
        // MultiFabs, or a FAB's ghost cells and its neighbor's, may overlap.
        long long owned = 0;
        if (!owned_boxes.isEmpty()) {
            amrex::BoxArray ba(owned_boxes);
            ba.removeOverlap();
            owned = ba.numPts();
        }
        const long long read = region.ok() ? region.numPts() : 0;
        struct { double ratio; int rank; } mine, worst;
        mine.ratio = (owned > 0) ? double(read) / double(owned) : 1.0;
        mine.rank = amrex::ParallelContext::MyProcSub();
        MPI_Comm comm = amrex::ParallelContext::CommunicatorSub();
        MPI_Allreduce(&mine, &worst, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
        double mean = mine.ratio;
        amrex::ParallelAllReduce::Sum(mean, comm);
        mean /= amrex::ParallelContext::NProcsSub();
        long long worst_cells[2] = {read, owned};
        MPI_Bcast(worst_cells, 2, MPI_LONG_LONG, worst.rank, comm);
        amrex::Print() << "NetCDF read region: up to " << worst.ratio
                       << " times the cells owned, on rank " << worst.rank << " ("
                       << worst_cells[0] << " read, " << worst_cells[1]
                       << " owned); mean " << mean << " over ranks" << std::endl;
    }
    return region;
}

//! nc_local_region for a fixed list of MultiFabs
inline amrex::Box
nc_local_region (std::initializer_list<const amrex::MultiFab*> mfs)
{
    return nc_local_region(amrex::Vector<const amrex::MultiFab*>(mfs));
}

//! nc_local_region for one MultiFab
inline amrex::Box
nc_local_region (const amrex::MultiFab& mf)
{
    return nc_local_region({&mf});
}

/**
 * \brief Read the part of each variable that this rank needs, collectively.
 *
 * Fills the same FABs as BuildFABsFromNetCDFFile, in the same index space,
 * but each rank's FABs cover only the file data inside \p region (converted to
 * each variable's staggering), and each rank reads only that part. A copy on
 * intersection into a MultiFab whose FABs this rank owns lie inside \p region
 * gives the same values as from BuildFABsFromNetCDFFile; nc_local_region gives
 * that region for one MultiFab.
 *
 * Collective over amrex::ParallelContext::CommunicatorSub(): every rank calls
 * it, including ranks with an empty \p region, which read nothing. Each
 * variable is one collective read of one hyperslab per rank, which works for
 * every format PnetCDF opens, compressed NetCDF-4 included. River data, which
 * is not laid out on the model grid, is not supported here.
 *
 * @param domain       domain whose lower corner the file data starts at
 * @param region       cell-centered region of the domain this rank needs
 * @param fname        file to read
 * @param nc_var_names variables to read
 * @param NC_dim_types dimension type of each variable in the file
 * @param fab_vars     FABs to fill, one per variable
 * @param one_time     whether to read a single time step
 * @param fill_time    time index to read when \p one_time
 * @param custom_ngrow grow cells, as in BuildFABsFromNetCDFFile
 */
template<class FAB,typename DType>
void
BuildLocalFABsFromNetCDFFile (const amrex::Box& domain,
                              const amrex::Box& region,
                              const std::string &fname,
                              amrex::Vector<std::string> nc_var_names,
                              amrex::Vector<enum NC_Data_Dims_Type> NC_dim_types,
                              amrex::Vector<FAB*> fab_vars,
                              bool one_time=false, int fill_time=0,
                              amrex::IntVect custom_ngrow=amrex::IntVect(1,1,0))
{
    auto ncf = ncutils::NCFile::open(fname, NC_NOWRITE, amrex::ParallelContext::CommunicatorSub());

    for (int iv = 0; iv < nc_var_names.size(); iv++)
    {
        const NC_Data_Dims_Type dim_type = NC_dim_types[iv];
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dim_type == NC_Data_Dims_Type::Time_BT_SN_WE ||
                                         dim_type == NC_Data_Dims_Type::Time_SN_WE ||
                                         dim_type == NC_Data_Dims_Type::BT_SN_WE ||
                                         dim_type == NC_Data_Dims_Type::SN_WE,
                                         "BuildLocalFABsFromNetCDFFile reads gridded data only");
        const bool has_time = (dim_type == NC_Data_Dims_Type::Time_BT_SN_WE ||
                               dim_type == NC_Data_Dims_Type::Time_SN_WE);
        const bool has_z    = (dim_type == NC_Data_Dims_Type::Time_BT_SN_WE ||
                               dim_type == NC_Data_Dims_Type::BT_SN_WE);

        const auto var = ncf.var(nc_var_names[iv]);
        const std::vector<MPI_Offset> shape = var.shape();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            static_cast<int>(shape.size()) == int(has_time) + int(has_z) + 2,
            "NetCDF variable " + nc_var_names[iv] + " in " + fname +
            " does not have the dimensions expected of it");

        // The same layout as BuildFABsFromNetCDFFile: the file data fill the
        // corner of the FAB box that starts at its lower corner.
        int nsx, nsy, nsz;
        amrex::Box fab_box = nc_fab_box(nc_var_names[iv], dim_type, shape, custom_ngrow,
                                        nsx, nsy, nsz);
        amrex::Box data_box(fab_box.smallEnd(),
                            fab_box.smallEnd() + amrex::IntVect(nsx-1, nsy-1, nsz-1),
                            fab_box.ixType());
        AMREX_ALWAYS_ASSERT(fab_box.contains(data_box));
        data_box.shift(domain.smallEnd());

        const amrex::Box local = region.ok()
            ? (data_box & amrex::convert(region, data_box.ixType())) : amrex::Box();

        // Hyperslab of the file holding `local`; empty on a rank that needs nothing.
        std::vector<MPI_Offset> start(shape.size(), 0), count(shape.size(), 0);
        int d = 0;
        if (has_time) {
            start[d] = one_time ? fill_time : 0;
            count[d] = local.ok() ? 1 : 0;
            ++d;
        }
        if (local.ok()) {
            if (has_z) {
                start[d] = local.smallEnd(2) - data_box.smallEnd(2);
                count[d] = local.length(2);
            }
            start[d + int(has_z)]     = local.smallEnd(1) - data_box.smallEnd(1);
            count[d + int(has_z)]     = local.length(1);
            start[d + int(has_z) + 1] = local.smallEnd(0) - data_box.smallEnd(0);
            count[d + int(has_z) + 1] = local.length(0);
        }

        // The file stores x fastest, then y, then z: the FAB's own order, so
        // the hyperslab reads straight into it. A rank with nothing to read
        // still gets a valid one-cell FAB.
        const amrex::Box fab_bx = local.ok() ? local
                                             : amrex::Box(data_box.smallEnd(), data_box.smallEnd(),
                                                          data_box.ixType());
        fab_vars[iv]->resize(fab_bx, 1);
#ifdef AMREX_USE_GPU
        // fab_vars live on the device; read on the host and copy over.
        FAB tmp(fab_bx, 1, amrex::The_Pinned_Arena());
        if (!local.ok()) { tmp.template setVal<amrex::RunOn::Host>(0.0); }
        var.get_all(tmp.dataPtr(), start, count);
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         tmp.dataPtr(), tmp.dataPtr() + tmp.size(),
                         fab_vars[iv]->dataPtr());
#else
        if (!local.ok()) { fab_vars[iv]->template setVal<amrex::RunOn::Host>(0.0); }
        var.get_all(fab_vars[iv]->dataPtr(), start, count);
#endif
    }
    ncf.close();
}

#endif
