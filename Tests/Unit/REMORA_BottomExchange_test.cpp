// Full solver column budgets, using an independent uniform initial state and
// analytic volume / salt / dye / horizontal momentum balances as the oracle.
#include <REMORA.H>
#include <REMORA_DepthStretchTransform.H>
#include <REMORA_prob_common.H>
#include <REMORA_BottomExchange.H>
#include <AMReX.H>
#include <AMReX_Reduce.H>
#include <algorithm>
#include <cmath>
#include <string>
#include <limits>

using namespace amrex;
std::string inputs_name;

namespace {
class Column : public ProblemBase {
public:
    void init_analytic_bathymetry (int, const Geometry&, const SolverChoice&,
        const REMORA&, MultiFab& h) override { h.setVal(10); }
    void init_analytic_zeta (int, const Geometry&, const SolverChoice&,
        const REMORA&, const ProbCoords&, MultiFab& z) override { z.setVal(0); }
    void init_analytic_masks (int, const Geometry&, const SolverChoice&,
        const REMORA&, MultiFab& mask) override { mask.setVal(1); }
    void init_analytic_coriolis (int, const Geometry&, const SolverChoice&,
        const REMORA&, MultiFab& f) override { f.setVal(0); }
    void init_analytic_grid_scale (int, const Geometry& g, const SolverChoice&,
        const REMORA&, MultiFab& pm, MultiFab& pn) override {
        pm.setVal(Real(1)/g.CellSize(0)); pn.setVal(Real(1)/g.CellSize(1));
    }
    void init_analytic_vmix (int, const Geometry&, const SolverChoice&,
        const REMORA&, MultiFab& km, MultiFab& kh) override { km.setVal(0); kh.setVal(0); }
    void init_analytic_hmix (int, const Geometry&, const SolverChoice&, const REMORA&,
        MultiFab& p, MultiFab& r, MultiFab& d) override { p.setVal(0); r.setVal(0); d.setVal(0); }
    void init_analytic_smflux (int, const Geometry&, const SolverChoice&, const REMORA&,
        MultiFab& x, MultiFab& y) override { x.setVal(0); y.setVal(0); }
    void init_analytic_prob (int, const Geometry&, const SolverChoice&, const REMORA&,
        const ProbCoords&, MultiFab& c, MultiFab& u, MultiFab& v) override {
        c.setVal(10,0,1); c.setVal(35,1,1); c.setVal(2,2,1);
        Real u0=Real(0.05); ParmParse("remora.prob").query("u_0",u0);
        u.setVal(u0); v.setVal(0);
    }
    void bottom_exchange (int, Real time, const Geometry&, const REMORA&, MultiFab& raw) override {
        ParmParse pp("remora"); Real q=0, salt=0, dye=7;
        pp.get("bottom_q",q); pp.query("bottom_salt",salt); pp.query("bottom_tracer",dye);
        bool ramp=false; ParmParse("test").query("ramp",ramp);
        if (ramp) { q *= Real(1)+time/Real(8); }
        raw.setVal(q,0,1); raw.setVal(10,1,1); raw.setVal(salt,2,1); raw.setVal(dye,3,1);
        bool spatial=false; ParmParse("test").query("spatial",spatial);
        if (spatial) {
            for (MFIter mfi(raw); mfi.isValid(); ++mfi) {
                const auto a=raw.array(mfi);
                ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i,int j,int) {
                    a(i,j,0,0)=q*(i<4 ? Real(0.5):Real(1.5));
                    a(i,j,0,3)=i<4 ? 2*dye : Real(0);
                });
            }
        }
    }
};

void require_close (const std::string& name, double actual, double expected,
                    double tolerance=2.e-10, double scale=1)
{
    // Allow for rounding in Real state while keeping sensitivity to the transfer.
    tolerance = std::max(tolerance, 16.0*std::numeric_limits<Real>::epsilon());
    Print() << name << ": actual=" << actual << " expected=" << expected << '\n';
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(actual) &&
        std::abs(actual-expected) <= tolerance*std::max(scale,std::abs(expected)), name);
}

Vector<double> budget (REMORA& ocean)
{
    Vector<double> total(5,0); // volume, temperature, salt, dye, u momentum
    for (int lev=0; lev<=ocean.finestLevel(); ++lev) {
        const MultiFab* fine = lev<ocean.finestLevel() ? &ocean.build_fine_mask(lev+1) : nullptr;
        ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum,ReduceOpSum> op;
        ReduceData<double,double,double,double,double> data(op);
        for (MFIter mfi(*ocean.cons_new[lev]); mfi.isValid(); ++mfi) {
            const auto c=ocean.cons_new[lev]->const_array(mfi);
            const auto u=ocean.xvel_new[lev]->const_array(mfi);
            const auto hz=ocean.vec_Hz[lev]->const_array(mfi);
            const auto pm=ocean.vec_pm[lev]->const_array(mfi), pn=ocean.vec_pn[lev]->const_array(mfi);
            const auto mask=fine ? fine->const_array(mfi) : Array4<Real const>();
            // Accumulate the oracle in double even when the evolved state is float.
            op.eval(mfi.validbox(),data,[=] AMREX_GPU_DEVICE (int i,int j,int k)
                -> GpuTuple<double,double,double,double,double> {
                const double volume=double(hz(i,j,k))/(double(pm(i,j,0))*double(pn(i,j,0))) *
                                    (mask ? double(mask(i,j,k)) : 1.0);
                return {volume,volume*double(c(i,j,k,0)),volume*double(c(i,j,k,1)),
                        volume*double(c(i,j,k,2)),
                        volume*0.5*(double(u(i,j,k))+double(u(i+1,j,k)))};
            });
        }
        const auto result=data.value(op);
        total[0]+=get<0>(result); total[1]+=get<1>(result); total[2]+=get<2>(result);
        total[3]+=get<3>(result); total[4]+=get<4>(result);
    }
    ParallelDescriptor::ReduceRealSum(total.data(),total.size());
    return total;
}

#ifdef REMORA_USE_NETCDF
void make_file (int kind, Real q, Real salt, const std::string& path="bottom.nc", int record=0)
{
    int nc, eta, xi, tim=-1, vars[4];
    auto check=[] (int status) { AMREX_ALWAYS_ASSERT_WITH_MESSAGE(status==NC_NOERR,ncmpi_strerror(status)); };
    check(ncmpi_create(ParallelDescriptor::Communicator(), path.c_str(), NC_CLOBBER, MPI_INFO_NULL, &nc));
    check(ncmpi_def_dim(nc,"eta_rho",10,&eta)); check(ncmpi_def_dim(nc,"xi_rho",10,&xi));
    const bool timed=kind!=1;
    if (timed) { check(ncmpi_def_dim(nc,"bottom_time",kind==2 ? 1 : 2,&tim)); }
    int td=-1;
    if (timed) {
        check(ncmpi_def_var(nc,"bottom_time",NC_DOUBLE,1,&tim,&td));
        const std::string units="days since 2000-01-01 00:00:00";
        check(ncmpi_put_att_text(nc,td,"units",units.size(),units.c_str()));
    }
    int dims[3]={tim,eta,xi};
    const char* names[4]={"bottom_q","bottom_temp","bottom_salt","bottom_tracer"};
    for (int n=0;n<4;++n) { check(ncmpi_def_var(nc,names[n],NC_DOUBLE,timed?3:2,timed?dims:dims+1,&vars[n])); }
    check(ncmpi_enddef(nc)); check(ncmpi_begin_indep_data(nc));
    if (ParallelDescriptor::IOProcessor()) {
        const double times[2]={record*16.0/86400.0,16.0/86400.0};
        if (timed) { check(ncmpi_put_var_double(nc,td,times)); }
        const int records=kind==3 ? 2 : 1;
        const double vals[4]={double(q),10,double(salt),7};
        for (int n=0;n<4;++n) {
            Vector<double> values(100*records,vals[n]);
            bool ring=false; ParmParse("test").query("ring",ring);
            if (ring && n==0) {
                for (int i=0;i<10;++i) { values[i]=2*vals[n]; values[90+i]=2*vals[n]; }
            }
            if (kind==3 && n==0) { for (int i=100;i<200;++i) { values[i]=3*vals[n]; } }
            check(ncmpi_put_var_double(nc,vars[n],values.data()));
        }
    }
    check(ncmpi_end_indep_data(nc)); check(ncmpi_close(nc));
}
#endif
}

std::unique_ptr<ProblemBase> amrex_probinit (const amrex_real*, const amrex_real*)
{
    return std::make_unique<Column>();
}

int main (int argc, char** argv)
{
    amrex::Initialize(argc,argv,true,MPI_COMM_WORLD,[] {
        ParmParse pp("remora"), amr("amr");
        Vector<int> cells,periods; Vector<Real> lo,hi;
        pp.getarr("n_cell",cells); pp.getarr("is_periodic",periods);
        pp.getarr("prob_lo",lo); pp.getarr("prob_hi",hi);
        amr.addarr("n_cell",cells); amr.add("refine_grid_layout_z",0); amr.add("max_grid_size_z",cells[2]);
        int max_grid=8; amr.queryAdd("max_grid_size",max_grid); amr.add("blocking_factor",1); amr.add("n_proper",2);
        ParmParse geom("geometry"); geom.addarr("is_periodic",periods);
        geom.addarr("prob_lo",lo); geom.addarr("prob_hi",hi);
    });
    {
        ParmParse pp("remora"), test("test");
        bool omit_salt=false; test.query("omit_salt",omit_salt);
        if (omit_salt) { pp.remove("bottom_salt"); }
        Real q=0, incoming_salt=0, incoming_dye=7; pp.query("bottom_q",q); pp.query("bottom_salt",incoming_salt); pp.query("bottom_tracer",incoming_dye);
        bool resume=false, reverse=false, unstable=false, spatial=false;
        test.query("resume",resume); test.query("reverse",reverse); test.query("unstable",unstable); test.query("spatial",spatial);
        std::string type; pp.get("bottom_exchange_type",type);
        int file_kind=0; test.query("file_kind",file_kind);
#ifdef REMORA_USE_NETCDF
        if (file_kind==4) {
            make_file(2,q,incoming_salt,"bottom0.nc",0);
            make_file(2,3*q,incoming_salt,"bottom1.nc",1);
        } else if (file_kind) { make_file(file_kind,q,incoming_salt); }
#endif
        REMORA ocean; ocean.InitData();
        if (type=="constant" || type=="analytic" || type=="netcdf") {
            for (int lev=0; lev<=ocean.finestLevel(); ++lev) {
                if (!spatial) {
                    require_close("initial vertical velocity", ocean.zvel_new[lev]->min(0), q,Real(2.e-10),Real(0.001));
                    require_close("initial old vertical velocity", ocean.zvel_old[lev]->min(0), q,Real(2.e-10),Real(0.001));
                }
            }
        }
        if (type=="coupler") {
            BoxArray ba; DistributionMapping dm; ocean.GetBottomExchangeLayout(ba,dm);
            MultiFab water(ba,dm,1,0), c(ba,dm,3,0);
            water.setVal(q); c.setVal(10,0,1); c.setVal(incoming_salt,1,1); c.setVal(incoming_dye,2,1);
            Real valid_end=100; test.query("valid_end",valid_end);
            ocean.SetBottomExchangeData(water,c,0,valid_end);
            // Mutating the producer's buffers must not mutate REMORA's owned snapshot.
            water.setVal(123); c.setVal(123);
        }
        bool ghosts=false; test.query("ghosts",ghosts);
        if (ghosts) {
            for (int lev=0; lev<=ocean.finestLevel(); ++lev) {
                const auto* transport=ocean.GetBottomExchangeTransport(lev);
                ReduceOps<ReduceOpMax> op;
                ReduceData<int> data(op);
                for (MFIter mfi(*transport); mfi.isValid(); ++mfi) {
                    const auto a=transport->const_array(mfi);
                    op.eval(mfi.fabbox(),data,[=] AMREX_GPU_DEVICE (int i,int j,int k) -> GpuTuple<int> {
                        return {std::isfinite(a(i,j,k)) ? 0 : 1};
                    });
                }
                int bad=get<0>(data.value(op));
                ParallelDescriptor::ReduceIntMax(bad);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!bad,"bottom transport ghosts must be finite");
            }
        }
        bool ring=false; test.query("ring",ring);
        if (ring) {
            const auto* transport=ocean.GetBottomExchangeTransport(0);
            for (MFIter mfi(*transport); mfi.isValid(); ++mfi) {
                const Box edge=makeSlab(mfi.fabbox(),1,-1) & mfi.fabbox();
                if (edge.ok()) {
                    const Real value=(*transport)[mfi].min<RunOn::Device>(edge,0);
                    require_close("NetCDF physical rho ring min",value,2*q*Real(100*100));
                    require_close("NetCDF physical rho ring max",
                        (*transport)[mfi].max<RunOn::Device>(edge,0),2*q*Real(100*100));
                }
            }
        }
        const auto initial=budget(ocean);
        const Real start_time = ocean.get_t_new(0);
        Real time=start_time, last_sample=start_time;
        bool ramp=false; test.query("ramp",ramp);
        for (int step=0;step<(resume ? 2 : 4);++step) {
            if (reverse && step==2) {
                BoxArray ba; DistributionMapping dm; ocean.GetBottomExchangeLayout(ba,dm);
                MultiFab water(ba,dm,1,0), c(ba,dm,3,0);
                water.setVal(-q); c.setVal(10,0,1); c.setVal(35,1,1); c.setVal(2,2,1);
                ocean.SetBottomExchangeData(water,c,time,100);
            }
            const Real advanced=ocean.EvolveOneStep(time,Real(2));
            if (advanced > 0) { last_sample=time+Real(0.5)*advanced; }
            time += advanced;

        }
        const auto final=budget(ocean);
        const Real duration = (file_kind>=3 || ramp) ? time-start_time+(time*time-start_time*start_time)/Real(16) : time-start_time;
        const Real added=(type=="none" || reverse ? Real(0):q)*Real(800*800)*duration;
        require_close("volume",final[0],initial[0]+added);
        require_close("temperature",final[1],Real(10)*final[0]);
        require_close("salt",final[2],q>=0 ? initial[2]+added*incoming_salt : Real(35)*final[0]);
        require_close("dye",final[3],q>=0 ? initial[3]+added*incoming_dye*(spatial ? Real(0.5) : Real(1)) : Real(2)*final[0]);
        require_close("u momentum",final[4],q>=0 ? initial[4] : Real(0.05)*final[0],Real(2.e-8));
        for (int lev=0;lev<=ocean.finestLevel();++lev) {
            ocean.compute_wvel(lev, *ocean.zvel_new[lev]);
            MultiFab bottom(ocean.vec_mskr[lev]->boxArray(),ocean.cons_new[lev]->DistributionMap(),3,0);
            ocean.PackBottomState(lev,bottom);
            require_close("bottom temperature",bottom.min(0),10);
            if (!spatial) {
                const Real imposed=type=="none" ? Real(0) : reverse ? -q :
                    (file_kind>=3 || ramp) ? q*(1+last_sample/Real(8)) : q;
                require_close("vertical velocity min",ocean.zvel_new[lev]->min(0),imposed,Real(1.e-9),Real(0.001));
                require_close("vertical velocity max",ocean.zvel_new[lev]->max(0),imposed,Real(1.e-9),Real(0.001));
            }
            AMREX_ALWAYS_ASSERT(bottom.min(1)>=0);
            if (incoming_salt==35 || q<0 || type=="none") {
                require_close("uniform salinity min",bottom.min(1),35);
                require_close("uniform salinity max",bottom.max(1),35);
            }
        }
        if (unstable) {
            AMREX_ALWAYS_ASSERT(ocean.vec_bvf[0]->min(0) < 0);
            AMREX_ALWAYS_ASSERT(std::isfinite(ocean.vec_tke[0]->max(0)));
        }
        // These checks verify the documented scaling against a neutral log-layer oracle.
        const Real friction=Real(0.01), speed=Real(1.e-5), d=Real(1);
        const Real k=friction*friction/Real(0.3), km=Real(0.4)*friction*d;
        const Real diss=friction*friction*friction/(Real(0.4)*d);
        const auto ratios=BottomExchange::scales(speed,d,km,km,km/Real(1.3),k,diss,
            friction*friction,Real(0.1));
        require_close("GLS transport scaling",ratios.transport,Real(0.00325),Real(2.e-10),Real(0));
        require_close("GLS adjustment scaling",ratios.adjustment,Real(1)/Real(750),Real(2.e-10),Real(0));
        Print() << "Bottom exchange column passed\n";
    }
    amrex::Finalize();
}
