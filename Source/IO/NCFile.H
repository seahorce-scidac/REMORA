#ifndef _NCFILE_H_
#define _NCFILE_H_

#include <sstream>
#include <string>
#include <ctime>
#include <atomic>

#include "AMReX_FArrayBox.H"
#include "AMReX_IArrayBox.H"
#include "NCInterface.H"

using PlaneVector = amrex::Vector<amrex::FArrayBox>;

enum class NC_Data_Dims_Type {
    Time_BT_SN_WE,
    Time_SN_WE,
    SN_WE,
};

//
// NDArray is the datatype designed to hold any data, including scalars, multidimensional
// arrays, that read from the NetCDF file.
//
// The data read from NetCDF file are stored in a continuous memory, and the data layout is descriped
// by using a vector (shape). AMRex Box can be constructed using the data shape information, and MultiFab
// data array can be setup using the data that stored in the NDArray.
//
template <typename DataType>
struct NDArray
{
   using DType = typename std::remove_const<DataType>::type;

   // constructor
   explicit NDArray (const std::string vname, const std::vector<size_t>& vshape)
                    : name{vname}, shape{vshape}, ref_counted{1}, owned{true} {
      data = new DType [this->ndim()];
   }

   // default constructor
   NDArray () : name{"null"}, ref_counted{1}, owned{false}, data{nullptr} {}

   // copy constructor
   NDArray (const NDArray& array) {
     name  = array.name;
     shape = array.shape;
     data  = array.data;
     owned = false;
     ref_counted.fetch_add(1, std::memory_order_relaxed);
   }

   // copy assignment
   NDArray& operator=(const NDArray& array) {
      name  = array.name;
      shape = array.shape;
      data  = array.data;
      owned = false;
      ref_counted.fetch_add(1, std::memory_order_relaxed);
      return *this;
    }

   // destructor
   ~NDArray () {
     ref_counted.fetch_sub(1, std::memory_order_acq_rel);
     if (ref_counted == 1 && owned) delete[] data;
   }

   // get the data pointer
   decltype(auto) get_data () {
     ref_counted.fetch_add(1, std::memory_order_relaxed);
     return data;
   }

   // get the variable name
   std::string get_vname () {
     return name;
   }

   // get the variable data shape
   std::vector<size_t> get_vshape () {
     return shape;
   }

   // return the total number of data
   size_t ndim () {
     size_t num = 1;
     int isize = static_cast<int>(shape.size());
     for (auto i=0; i<isize; ++i) num *= shape[i];
     return num;
   }

   // set the data shape information
   void set_vshape (std::vector<size_t> vshape) {
     shape = vshape;
    }

 private:
   std::string name;
   std::vector<size_t> shape;
   mutable std::atomic<size_t> ref_counted;
   bool owned;
   DType* data;
};

template<typename DType>
void ReadNetCDFFile (const std::string& fname, amrex::Vector<std::string> names,
                     amrex::Vector<NDArray<DType> >& arrays)
{
    AMREX_ASSERT(arrays.size() == names.size());

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        auto ncf = ncutils::NCFile::open(fname, NC_CLOBBER | NC_NETCDF4);

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
        int bdy_width          = static_cast<int>(ncf.dim("bdy_width").len());
        */

        amrex::Print() << "Reading the dimensions from the netcdf file " << "\n";

        for (auto n=0; n<arrays.size(); ++n) {
            // read the data from NetCDF file
            std::string vname_to_write = names[n];
            std::string vname_to_read  = names[n];

            // amrex::AllPrint() << "About to read " << vname_to_read
            //                   << " while filling the array for " << vname_to_write << std::endl;

            std::vector<size_t> shape = ncf.var(vname_to_read).shape();
            arrays[n]                 = NDArray<DType>(vname_to_read,shape);
            DType* dataPtr            = arrays[n].get_data();

            std::vector<size_t> start(shape.size(), 0);

            // auto numPts               = arrays[n].ndim();
            // amrex::Print() << "NetCDF Variable name = " << vname_to_read << std::endl;
            // amrex::Print() << "numPts read from NetCDF file/var = " << numPts << std::endl;
            // amrex::Print()  << "Dims of the variable = (";
            // for (auto &dim:shape)
            //     amrex::Print() << dim << ", " ;
            // amrex::Print()  << ")" << std::endl;

            ncf.var(vname_to_read).get(dataPtr, start, shape);
        }
        ncf.close();
    }
}

/**
 * Helper function for reading a single variable attribute
 *
 * @param fname Name of NetCDF file
 * @param var_name Name of variable
 * @param attr_name Name of attribute to read
 * @returns attribute value
 */
std::string ReadNetCDFVarAttrStr (const std::string& fname,
                                  const std::string& var_name,
                                  const std::string& attr_name);

/**
 * Helper function for reading data from NetCDF file into a
 * provided FAB.
 *
 * @param iv Index for which variable we are going to fill
 * @param nc_arrays Arrays of data from NetCDF file
 * @param var_name Variable name
 * @param NC_dim_type Dimension type for the variable as stored in the NetCDF file
 * @param temp FAB where we store the variable data from the NetCDF Arrays
 */
template<class FAB,typename DType>
void
fill_fab_from_arrays (int iv,
                      amrex::Vector<NDArray<float>>& nc_arrays,
                      const std::string& var_name,
                      NC_Data_Dims_Type& NC_dim_type,
                      FAB& temp)
{
    int nsx, nsy, nsz;
    if (NC_dim_type == NC_Data_Dims_Type::SN_WE) {
        nsz = 1;
        nsy = nc_arrays[iv].get_vshape()[0];
        nsx = nc_arrays[iv].get_vshape()[1];
        // amrex::Print() << "TYPE SN WE " << nsy << " " << nsx << std::endl;
    } else if (NC_dim_type == NC_Data_Dims_Type::Time_SN_WE) {
        nsz = 1;
        nsy = nc_arrays[iv].get_vshape()[1];
        nsx = nc_arrays[iv].get_vshape()[2];
        // amrex::Print() << "TYPE TIME SN WE " << nsy << " " << nsx << std::endl;
    } else if (NC_dim_type == NC_Data_Dims_Type::Time_BT_SN_WE) {
        nsz = nc_arrays[iv].get_vshape()[1];
        nsy = nc_arrays[iv].get_vshape()[2];
        nsx = nc_arrays[iv].get_vshape()[3];
        // amrex::Print() << "TYPE TIME BT SN WE " << nsx << " " << nsy << " " << nsz << std::endl;
    } else {
        amrex::Abort("Dont know this NC_Data_Dims_Type");
    }

    amrex::Box my_box;

    if (var_name == "u" || var_name == "ubar" || var_name == "mask_u")
    {
        my_box.setSmall(amrex::IntVect(0,-1,0));
        my_box.setBig(amrex::IntVect(nsx-1,nsy-2,nsz-1));
        my_box.setType(amrex::IndexType(amrex::IntVect(1,0,0)));
    }
    else if (var_name == "v" || var_name == "vbar" || var_name == "mask_v")
    {
        my_box.setSmall(amrex::IntVect(-1,0,0));
        my_box.setBig(amrex::IntVect(nsx-2,nsy-1,nsz-1));
        my_box.setType(amrex::IndexType(amrex::IntVect(0,1,0)));
    }
    else  // everything cell-centered -- these have one ghost cell
    {
        my_box.setSmall(amrex::IntVect(-1,-1,0));
        my_box.setBig(amrex::IntVect(nsx-2,nsy-2,nsz-1));
    }

    // amrex::Print() <<" NAME and BOX for this VAR " << var_name << " " << my_box << std::endl;

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
 * Function to read NetCDF variables and fill the corresponding Array4's
 *
 * @param fname Name of the NetCDF file to be read
 * @param nc_var_names Variable names in the NetCDF file
 * @param NC_dim_types NetCDF data dimension types
 * @param fab_vars Fab data we are to fill
 */
template<class FAB,typename DType>
void
BuildFABsFromNetCDFFile (const amrex::Box& domain,
                         const std::string &fname,
                         amrex::Vector<std::string> nc_var_names,
                         amrex::Vector<enum NC_Data_Dims_Type> NC_dim_types,
                         amrex::Vector<FAB*> fab_vars)
{
    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();  // I/O rank

    amrex::Vector<NDArray<float>> nc_arrays(nc_var_names.size());

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        ReadNetCDFFile(fname, nc_var_names, nc_arrays);
    }

    for (int iv = 0; iv < nc_var_names.size(); iv++)
    {
        FAB tmp;
        if (amrex::ParallelDescriptor::IOProcessor()) {
            fill_fab_from_arrays<FAB,DType>(iv, nc_arrays, nc_var_names[iv], NC_dim_types[iv], tmp);
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

#endif
