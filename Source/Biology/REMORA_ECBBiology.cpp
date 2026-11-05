/*
Copyright (c) 2002-2026 Regional Ocean Modeling System (ROMS)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
// Native port of ROMS ECB; see Docs/sphinx_doc/ECB.rst for pinned
// source/options.
#include <REMORA_ECBBiology.H>
#include <REMORA_BiologyCarbonate.H>
#include <REMORA_DataStruct.H>
#include <REMORA_Constants.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParmParse.H>
#include <cmath>
#include <iterator>
using namespace amrex;
namespace {
// Local scratch inventory follows ROMS idbio, after temperature and salinity.
enum ECBComponent {
    sc_temp,
    sc_salt,
    sc_no3,
    sc_nh4,
    sc_chlo,
    sc_phyt,
    sc_zoop,
    sc_lden,
    sc_sden,
    sc_don,
    sc_po4,
    sc_ldec,
    sc_sdec,
    sc_tic,
    sc_talk,
    sc_doc,
    sc_oxyg,
    sc_rdon,
    sc_rdoc,
    sc_calc
};
constexpr int NumScratch = sc_calc + 1;
constexpr int NumECBTracers = NumScratch - sc_no3;
// State tracer name of each scratch biology component, sc_no3 to sc_calc. The
// kernel finds each one's state component by name, so the declared order is
// free to change.
const char* const ecb_tracer_names[NumECBTracers] = {
  "NO3",        "NH4",        "chlorophyll", "phytoplankton", "zooplankton",
  "LdetritusN", "SdetritusN", "DON",         "PO4",           "LdetritusC",
  "SdetritusC", "TIC",        "alkalinity",  "DOC",           "oxygen",
  "Rdon",       "Rdoc",       "calcium"};
} // namespace
namespace REMORABiology {
void
ECBBiology::configure(const std::string& remora_prefix)
{
    ParmParse pp(remora_prefix + ".ecb");
    for (const char* option :
         {"carbon", "oxygen", "bio_sediment", "denitrification", "po4",
          "talk_nonconserv", "co2_schmidt", "oxygen_schmidt", "pco2air_type",
          "pco2_rz", "sediment_proxy", "net_sulfate_reduction", "wet_dry"}) {
        if (pp.contains(option))
            Abort(
              "ECB implements a fixed ROMS option set; unsupported option: " +
              remora_prefix + ".ecb." + option);
    }
    pp.query("BioIter", m_parameters.BioIter);
    pp.get("AttSW", m_parameters.AttSW);
    pp.get("AttChl", m_parameters.AttChl);
    pp.get("Chl2C_m", m_parameters.Chl2C_m);
    pp.get("ChlMin", m_parameters.ChlMin);
    pp.get("CoagR", m_parameters.CoagR);
    pp.get("D_p5NH4", m_parameters.D_p5NH4);
    pp.get("I_thNH4", m_parameters.I_thNH4);
    pp.get("K_NH4", m_parameters.K_NH4);
    pp.get("K_NO3", m_parameters.K_NO3);
    pp.get("K_PO4", m_parameters.K_PO4);
    pp.get("K_Phy", m_parameters.K_Phy);
    pp.get("LDeRRN", m_parameters.LDeRRN);
    pp.get("LDeRRC", m_parameters.LDeRRC);
    pp.get("NitriR", m_parameters.NitriR);
    pp.get("PARfrac", m_parameters.PARfrac);
    pp.get("PhyCN", m_parameters.PhyCN);
    pp.get("R_P2N", m_parameters.R_P2N);
    pp.get("PhyIS", m_parameters.PhyIS);
    pp.get("PhyMin", m_parameters.PhyMin);
    pp.get("PhyMR", m_parameters.PhyMR);
    pp.get("SDeRRN", m_parameters.SDeRRN);
    pp.get("SDeRRC", m_parameters.SDeRRC);
    pp.get("RDeRRN", m_parameters.RDeRRN);
    pp.get("RDeRRC", m_parameters.RDeRRC);
    pp.get("Vp0", m_parameters.Vp0);
    pp.get("wLDet", m_parameters.wLDet);
    pp.get("wPhy", m_parameters.wPhy);
    pp.get("wSDet", m_parameters.wSDet);
    pp.get("ZooAE_N", m_parameters.ZooAE_N);
    pp.get("ZooBM", m_parameters.ZooBM);
    pp.get("ZooCN", m_parameters.ZooCN);
    pp.get("ZooER", m_parameters.ZooER);
    pp.get("ZooGR", m_parameters.ZooGR);
    pp.get("ZooMin", m_parameters.ZooMin);
    pp.get("ZooMR", m_parameters.ZooMR);
    pp.get("ElDON", m_parameters.ElDON);
    pp.get("EsDON", m_parameters.EsDON);
    pp.get("deltN", m_parameters.deltN);
    pp.get("deltC", m_parameters.deltC);
    pp.get("gammaC", m_parameters.gammaC);
    pp.get("slCexc", m_parameters.slCexc);
    pp.get("pCO2air", m_parameters.pCO2air);
    // Reject values that divide by zero in the kernel, where they would
    // otherwise surface as NaN in the tracers with no diagnostic.
    auto require =
      [&](bool ok, const char* name, const char* condition, Real value) {
          if (!ok) {
              Abort(
                remora_prefix + ".ecb." + name + " must be " + condition +
                "; got " + std::to_string(value));
          }
      };
    const auto& prm = m_parameters;
    // Light attenuation divides by (AttSW + AttChl*chl + DOC term) * dz.
    require(prm.AttSW > Real(0.0), "AttSW", "positive", prm.AttSW);
    require(prm.AttChl >= Real(0.0), "AttChl", "non-negative", prm.AttChl);
    // Carbon excess uptake divides by slCexc.
    require(prm.slCexc > Real(0.0), "slCexc", "positive", prm.slCexc);
    // Grazing divides by K_Phy + phyt^2, which is K_Phy when phyt is zero.
    require(prm.K_Phy > Real(0.0), "K_Phy", "positive", prm.K_Phy);
    m_declaration.tracer_names.assign(
      std::begin(ecb_tracer_names), std::end(ecb_tracer_names));
    m_declaration.plot_aliases = {"ecb"};
}
void
ECBBiology::advance(const BiologyAdvanceContext& context) const
{
    if (m_parameters.BioIter <= 0)
        return;
    const bool bulk = context.solver_choice.bulk_fluxes;
    if (
      !context.hz || !context.z_w || !context.mask || !context.shortwave ||
      !context.bottom_stress_u || !context.bottom_stress_v ||
      (bulk && (!context.wind_u || !context.wind_v)) ||
      (!bulk && (!context.stress_u || !context.stress_v)))
        Abort("ECB requires Hz, z_w, mask, shortwave, bottom stresses and "
              "surface wind/stress");
    const int nz = context.top_layer + 1;
    if (nz < 2)
        Abort("ECB sinking requires at least two vertical layers");
    const int base = context.component_start;
    // State component of each scratch component: temperature and salinity,
    // then each biology tracer looked up by name in the declaration.
    GpuArray<int, NumScratch> state_comp;
    state_comp[sc_temp] = Temp_comp;
    state_comp[sc_salt] = Salt_comp;
    for (int n = 0; n < NumECBTracers; ++n) {
        state_comp[sc_no3 + n] =
          m_declaration.component_index(ecb_tracer_names[n], base);
        AMREX_ALWAYS_ASSERT(state_comp[sc_no3 + n] >= 0);
    }
    const Real dt = context.dt;
    if (dt <= 0)
        Abort("ECB requires a positive timestep when BioIter > 0");
    const Real rho0 = context.solver_choice.rho0;
    const auto p = m_parameters;
    for (MFIter mfi(context.new_state, false); mfi.isValid(); ++mfi) {
        Box column = mfi.tilebox();
        if (column.smallEnd(2) != 0 || column.bigEnd(2) != nz - 1)
            Abort("ECB requires each grid box to contain the entire vertical "
                  "column");
        Box horizontal = column;
        horizontal.makeSlab(2, 0);
        Box work = column;
        work.setBig(2, nz);
        // Bio and Bio_old use zero-based physical layers; the column work
        // arrays retain ROMS's indices 0..nz for an auditable sinking stencil.
        FArrayBox storage(work, 49, The_Async_Arena());
        IArrayBox source(work, 1, The_Async_Arena());
        const auto bio = storage.array(0);
        const auto bio_old = storage.array(20);
        const auto fc = storage.array(40);
        const auto hz_inv = storage.array(41);
        const auto hz_inv2 = storage.array(42);
        const auto hz_inv3 = storage.array(43);
        const auto wl = storage.array(44);
        const auto wr = storage.array(45);
        const auto bl = storage.array(46);
        const auto br = storage.array(47);
        const auto qc = storage.array(48);
        const auto ksource = source.array();
        const auto old = context.old_state.const_array(mfi);
        const auto next = context.new_state.array(mfi);
        const auto hz = context.hz->const_array(mfi);
        const auto z_w = context.z_w->const_array(mfi);
        const auto mask = context.mask->const_array(mfi);
        const auto srflx = context.shortwave->const_array(mfi);
        const auto bustr = context.bottom_stress_u->const_array(mfi);
        const auto bvstr = context.bottom_stress_v->const_array(mfi);
        Array4<const Real> sustr, svstr, uwind, vwind;
        if (bulk) {
            uwind = context.wind_u->const_array(mfi);
            vwind = context.wind_v->const_array(mfi);
        } else {
            sustr = context.stress_u->const_array(mfi);
            svstr = context.stress_v->const_array(mfi);
        }
        ParallelFor(
          horizontal, [=] AMREX_GPU_DEVICE(int i, int j, int) noexcept {
              if (mask(i, j, 0) == Real(0.0))
                  return;
              constexpr int nsink = 6;
              int iter, ibio, isink, itrc, k, ks;
              constexpr Real eps = Real(1.0e-20);
              Real u10squ;
              constexpr Real a_o2 = Real(1920.4);
              constexpr Real b_o2 = Real(135.6);
              constexpr Real c_o2 = Real(5.2122);
              constexpr Real d_o2 = Real(0.10939);
              constexpr Real e_o2 = Real(0.00093777);
              constexpr Real oa0 = Real(2.00907);
              constexpr Real oa1 = Real(3.22014);
              constexpr Real oa2 = Real(4.05010);
              constexpr Real oa3 = Real(4.94457);
              constexpr Real oa4 = -Real(0.256847);
              constexpr Real oa5 = Real(3.88767);
              constexpr Real ob0 = -Real(0.00624523);
              constexpr Real ob1 = -Real(0.00737614);
              constexpr Real ob2 = -Real(0.0103410);
              constexpr Real ob3 = -Real(0.00817083);
              constexpr Real oc0 = -Real(0.000000488682);
              constexpr Real roxno3 = Real(8.625);
              constexpr Real roxnh4 = Real(6.625);
              constexpr Real k_ntr = Real(1.0);
              constexpr Real k_dnf = k_ntr;
              constexpr Real k_wno3 = Real(0.5);
              constexpr Real k_bo2 = Real(26.5);
              Real l2mol = Real(1000.0) / Real(22.3916);
              constexpr Real a_co2 = Real(2116.8);
              constexpr Real b_co2 = Real(136.25);
              constexpr Real c_co2 = Real(4.7353);
              constexpr Real d_co2 = Real(0.092307);
              constexpr Real e_co2 = Real(0.0007555);
              constexpr Real a1 = -Real(60.2409);
              constexpr Real a2 = Real(93.4517);
              constexpr Real a3 = Real(23.3585);
              constexpr Real b1 = Real(0.023517);
              constexpr Real b2 = -Real(0.023656);
              constexpr Real b3 = Real(0.0047036);
              Real att, attfac, expatt, itop, par;
              Real epp, l_nh4, l_no3, ltot, vp;
              constexpr Real minval = Real(1.0e-6);
              Real l_po4, lmin, mu;
              Real chl2c, dtdays, t_ppmax, inhnh4;
              Real cff, cff1, cff2, cff3, cff4, cff5, cff6;
              Real cff7, cff8;
              Real fac, fac1, fac2, fac3;
              Real cffl, cffr, cu, dltl, dltr;
              Real schmidtn_ox, o2satu, o2_flux;
              Real ts, aa;
              Real c_flux_reminel, c_flux_remines, c_flux_reminer;
              Real co2_flux, co2_sol, schmidtn, tempk;
              Real f_ntr, f_dnf;
              Real lbo2, o2satb;
              Real n_flux_assim;
              Real n_flux_coagd, n_flux_coagp;
              Real n_flux_egest;
              Real n_flux_newprod, n_flux_regprod;
              Real n_flux_nitrifi;
              Real n_flux_pmortal, n_flux_zmortal;
              Real n_flux_reminel, n_flux_remines, n_flux_reminer;
              Real n_flux_zexcret, n_flux_zmetabo;
              // cbe is always set before use; initialize it for
              // -Wmaybe-uninitialized.
              Real resuspr, cbe = Real(0.0), ustarb;
              Real c_exc_up, qdoc, rdoc, cnbur;
              Real g_gmax, qdon, rdon;
              Real n_flux_sloppy_ldon, n_flux_sloppy_sldon;
              int idsink[6];
              Real wbio[6];
              Real parsur, pco2;
              // Subcycle the complete ordered ECB reaction/column calculation.
              dtdays = dt * Real(1.0 / 86400.0) / Real(p.BioIter);
              idsink[0] = sc_phyt;
              idsink[1] = sc_chlo;
              idsink[2] = sc_sden;
              idsink[3] = sc_lden;
              idsink[4] = sc_sdec;
              idsink[5] = sc_ldec;
              wbio[1 - 1] = p.wPhy;
              wbio[2 - 1] = p.wPhy;
              wbio[3 - 1] = p.wSDet;
              wbio[4 - 1] = p.wLDet;
              wbio[5 - 1] = p.wSDet;
              wbio[6 - 1] = p.wLDet;
              for (k = 1; k <= nz; k += 1) {
                  hz_inv(i, j, k) = Real(1.0) / hz(i, j, k - 1);
              }
              for (k = 1; k <= nz - 1; k += 1) {
                  hz_inv2(i, j, k) =
                    Real(1.0) / (hz(i, j, k - 1) + hz(i, j, k + 1 - 1));
              }
              for (k = 2; k <= nz - 1; k += 1) {
                  hz_inv3(i, j, k) =
                    Real(1.0) / (hz(i, j, k - 1 - 1) + hz(i, j, k - 1) +
                                 hz(i, j, k + 1 - 1));
              }
              // Reference extraction applies positivity and TIC bounds to
              // Bio_old.
              for (itrc = 1; itrc <= 18; itrc += 1) {
                  ibio = (itrc + 1);
                  for (k = 1; k <= nz; k += 1) {
                      bio_old(i, j, k - 1, ibio) = amrex::max(
                        Real(0.0),
                        old(i, j, k - 1, state_comp[ibio]));
                      bio(i, j, k - 1, ibio) = bio_old(i, j, k - 1, ibio);
                  }
              }
              for (k = 1; k <= nz; k += 1) {
                  bio_old(i, j, k - 1, sc_tic) =
                    amrex::min(bio_old(i, j, k - 1, sc_tic), Real(3000.0));
                  bio_old(i, j, k - 1, sc_tic) =
                    amrex::max(bio_old(i, j, k - 1, sc_tic), Real(1.0));
                  bio(i, j, k - 1, sc_tic) = bio_old(i, j, k - 1, sc_tic);
              }
              for (k = 1; k <= nz; k += 1) {
                  bio(i, j, k - 1, sc_temp) =
                    amrex::min(old(i, j, k - 1, 0), Real(35.0));
                  bio(i, j, k - 1, sc_salt) =
                    amrex::max(old(i, j, k - 1, 1), Real(0.0));
              }
              // REMORA stores srflx in W/m2, so only PARfrac applies. ROMS
              // ecb.h multiplies by rho0*Cp because its srflx is kinematic
              // (degC m/s).
              parsur = p.PARfrac * srflx(i, j, 0);
              // Primary production, exudation and nitrification, from surface
              // down.
              for (iter = 1; iter <= p.BioIter; iter += 1) {
                  par = parsur;
                  if (parsur > Real(0.0)) {
                      for (k = nz; k >= 1; k += -1) {
                          attfac = amrex::max(
                                     bio(i, j, k - 1, sc_doc) +
                                       bio(i, j, k - 1, sc_rdoc) - Real(70.819),
                                     Real(0.0)) *
                                   Real(3.786e-3);
                          att =
                            (p.AttSW + p.AttChl * bio(i, j, k - 1, sc_chlo) +
                             attfac) *
                            (z_w(i, j, k) - z_w(i, j, k - 1));
                          expatt = std::exp(-att);
                          itop = par;
                          par = itop * (Real(1.0) - expatt) / att;
                          cff = p.PhyCN * Real(12.0);
                          chl2c = amrex::min(
                            bio(i, j, k - 1, sc_chlo) /
                              (bio(i, j, k - 1, sc_phyt) * cff + eps),
                            p.Chl2C_m);
                          vp = amrex::max(
                            p.Vp0 *
                              std::exp(
                                Real(0.08065) * bio(i, j, k - 1, sc_temp)),
                            Real(4.00));
                          fac1 = par * p.PhyIS;
                          epp = vp / std::sqrt(vp * vp + fac1 * fac1);
                          t_ppmax = epp * fac1;
                          cff1 = bio(i, j, k - 1, sc_nh4) * p.K_NH4;
                          cff2 = bio(i, j, k - 1, sc_no3) * p.K_NO3;
                          inhnh4 = Real(1.0) / (Real(1.0) + cff1);
                          l_nh4 = cff1 / (Real(1.0) + cff1);
                          l_no3 = cff2 * inhnh4 / (Real(1.0) + cff2);
                          ltot = l_no3 + l_nh4;
                          cff3 = bio(i, j, k - 1, sc_po4) * p.K_PO4;
                          l_po4 = cff3 / (Real(1.0) + cff3);
                          lmin = amrex::min(ltot, l_po4);
                          mu = dtdays * t_ppmax * lmin;
                          cff4 = mu * bio(i, j, k - 1, sc_phyt) * l_no3 /
                                 amrex::max(minval, ltot) /
                                 amrex::max(minval, bio(i, j, k - 1, sc_no3));
                          cff5 = mu * bio(i, j, k - 1, sc_phyt) * l_nh4 /
                                 amrex::max(minval, ltot) /
                                 amrex::max(minval, bio(i, j, k - 1, sc_nh4));
                          cff6 = p.R_P2N * mu * bio(i, j, k - 1, sc_phyt) /
                                 amrex::max(minval, bio(i, j, k - 1, sc_po4));
                          bio(i, j, k - 1, sc_no3) =
                            bio(i, j, k - 1, sc_no3) / (Real(1.0) + cff4);
                          bio(i, j, k - 1, sc_nh4) =
                            bio(i, j, k - 1, sc_nh4) / (Real(1.0) + cff5);
                          bio(i, j, k - 1, sc_po4) =
                            bio(i, j, k - 1, sc_po4) / (Real(1.0) + cff6);
                          n_flux_newprod = bio(i, j, k - 1, sc_no3) * cff4;
                          n_flux_regprod = bio(i, j, k - 1, sc_nh4) * cff5;
                          f_ntr = bio(i, j, k - 1, sc_oxyg) /
                                  (bio(i, j, k - 1, sc_oxyg) + k_ntr);
                          f_dnf = k_dnf / (bio(i, j, k - 1, sc_oxyg) + k_dnf) *
                                  bio(i, j, k - 1, sc_no3) /
                                  (bio(i, j, k - 1, sc_no3) + k_wno3);
                          bio(i, j, k - 1, sc_phyt) =
                            bio(i, j, k - 1, sc_phyt) +
                            (Real(1.0) - p.EsDON - (f_ntr + f_dnf) * p.ElDON) *
                              (n_flux_newprod + n_flux_regprod);
                          bio(i, j, k - 1, sc_chlo) =
                            bio(i, j, k - 1, sc_chlo) +
                            (dtdays * t_ppmax * t_ppmax * lmin * lmin *
                             p.Chl2C_m * bio(i, j, k - 1, sc_chlo)) /
                              (p.PhyIS * amrex::max(chl2c, eps) * par + eps);
                          bio(i, j, k - 1, sc_don) =
                            bio(i, j, k - 1, sc_don) +
                            p.EsDON * (n_flux_newprod + n_flux_regprod);
                          bio(i, j, k - 1, sc_nh4) =
                            bio(i, j, k - 1, sc_nh4) +
                            (f_ntr + f_dnf) * p.ElDON *
                              (n_flux_newprod + n_flux_regprod);
                          bio(i, j, k - 1, sc_talk) =
                            bio(i, j, k - 1, sc_talk) +
                            (f_ntr + f_dnf) * p.ElDON *
                              (n_flux_newprod + n_flux_regprod);
                          bio(i, j, k - 1, sc_oxyg) =
                            bio(i, j, k - 1, sc_oxyg) +
                            roxno3 * n_flux_newprod + roxnh4 * n_flux_regprod -
                            roxnh4 * f_ntr * p.ElDON *
                              (n_flux_newprod + n_flux_regprod);
                          bio(i, j, k - 1, sc_oxyg) =
                            amrex::max(bio(i, j, k - 1, sc_oxyg), Real(0.0));
                          cff1 = p.PhyCN * (n_flux_newprod + n_flux_regprod);
                          cff1 =
                            amrex::min(cff1, bio(i, j, k - 1, sc_tic) - eps);
                          bio(i, j, k - 1, sc_tic) =
                            bio(i, j, k - 1, sc_tic) - cff1;
                          c_exc_up = p.gammaC * p.PhyCN * dtdays * t_ppmax *
                                     bio(i, j, k - 1, sc_phyt) *
                                     (Real(1.0) - ltot);
                          c_exc_up = amrex::min(
                            c_exc_up,
                            amrex::max(bio(i, j, k - 1, sc_tic), Real(0.0)) /
                              p.slCexc);
                          bio(i, j, k - 1, sc_tic) =
                            bio(i, j, k - 1, sc_tic) - c_exc_up * p.slCexc;
                          bio(i, j, k - 1, sc_doc) =
                            bio(i, j, k - 1, sc_doc) + c_exc_up * p.slCexc;
                          bio(i, j, k - 1, sc_oxyg) =
                            bio(i, j, k - 1, sc_oxyg) + c_exc_up;
                          bio(i, j, k - 1, sc_doc) =
                            bio(i, j, k - 1, sc_doc) +
                            p.PhyCN * p.EsDON *
                              (n_flux_newprod + n_flux_regprod);
                          bio(i, j, k - 1, sc_tic) =
                            bio(i, j, k - 1, sc_tic) +
                            p.PhyCN * p.ElDON * (f_ntr + f_dnf) *
                              (n_flux_newprod + n_flux_regprod);
                          bio(i, j, k - 1, sc_talk) =
                            bio(i, j, k - 1, sc_talk) + n_flux_newprod -
                            n_flux_regprod;
                          fac2 =
                            amrex::max(bio(i, j, k - 1, sc_oxyg), Real(0.0));
                          fac3 = amrex::max(fac2 / (k_ntr + fac2), Real(0.0));
                          fac1 = dtdays * p.NitriR * fac3;
                          cff1 = (par - p.I_thNH4) /
                                 (p.D_p5NH4 + par - Real(2.0) * p.I_thNH4);
                          cff2 = Real(1.0) - amrex::max(Real(0.0), cff1);
                          cff3 = fac1 * cff2;
                          bio(i, j, k - 1, sc_nh4) =
                            bio(i, j, k - 1, sc_nh4) / (Real(1.0) + cff3);
                          n_flux_nitrifi = bio(i, j, k - 1, sc_nh4) * cff3;
                          bio(i, j, k - 1, sc_no3) =
                            bio(i, j, k - 1, sc_no3) + n_flux_nitrifi;
                          bio(i, j, k - 1, sc_oxyg) =
                            bio(i, j, k - 1, sc_oxyg) -
                            Real(2.0) * n_flux_nitrifi;
                          bio(i, j, k - 1, sc_talk) =
                            bio(i, j, k - 1, sc_talk) -
                            Real(2.0) * n_flux_nitrifi;
                          par = itop * expatt;
                      }
                  } else {
                      for (k = nz; k >= 1; k += -1) {
                          f_ntr = bio(i, j, k - 1, sc_oxyg) /
                                  (bio(i, j, k - 1, sc_oxyg) + k_ntr);
                          cff3 = dtdays * p.NitriR * f_ntr;
                          bio(i, j, k - 1, sc_nh4) =
                            bio(i, j, k - 1, sc_nh4) / (Real(1.0) + cff3);
                          n_flux_nitrifi = bio(i, j, k - 1, sc_nh4) * cff3;
                          bio(i, j, k - 1, sc_no3) =
                            bio(i, j, k - 1, sc_no3) + n_flux_nitrifi;
                          bio(i, j, k - 1, sc_oxyg) =
                            bio(i, j, k - 1, sc_oxyg) -
                            Real(2.0) * n_flux_nitrifi;
                          bio(i, j, k - 1, sc_oxyg) =
                            amrex::max(Real(0.0), bio(i, j, k - 1, sc_oxyg));
                          bio(i, j, k - 1, sc_talk) =
                            bio(i, j, k - 1, sc_talk) -
                            Real(2.0) * n_flux_nitrifi;
                      }
                  }
                  // Grazing and phytoplankton mortality.
                  cff2 = dtdays * p.PhyMR;
                  for (k = 1; k <= nz; k += 1) {
                      fac1 = dtdays * p.ZooGR *
                             std::exp(Real(0.0875) * bio(i, j, k - 1, sc_temp));
                      cff1 = fac1 * bio(i, j, k - 1, sc_zoop) *
                             bio(i, j, k - 1, sc_phyt) /
                             (p.K_Phy + bio(i, j, k - 1, sc_phyt) *
                                          bio(i, j, k - 1, sc_phyt));
                      cff3 = Real(1.0) / (Real(1.0) + cff1);
                      bio(i, j, k - 1, sc_phyt) =
                        cff3 * bio(i, j, k - 1, sc_phyt);
                      bio(i, j, k - 1, sc_chlo) =
                        cff3 * bio(i, j, k - 1, sc_chlo);
                      n_flux_assim =
                        cff1 * bio(i, j, k - 1, sc_phyt) * p.ZooAE_N;
                      n_flux_egest = bio(i, j, k - 1, sc_phyt) * cff1 *
                                     (Real(1.0) - p.ZooAE_N);
                      bio(i, j, k - 1, sc_zoop) =
                        bio(i, j, k - 1, sc_zoop) + n_flux_assim;
                      g_gmax = bio(i, j, k - 1, sc_phyt) *
                               bio(i, j, k - 1, sc_phyt) /
                               (p.K_Phy + bio(i, j, k - 1, sc_phyt) *
                                            bio(i, j, k - 1, sc_phyt));
                      qdon = Real(0.71) * g_gmax;
                      bio(i, j, k - 1, sc_lden) =
                        bio(i, j, k - 1, sc_lden) +
                        n_flux_egest * (Real(1.0) - qdon);
                      n_flux_sloppy_ldon =
                        n_flux_egest * qdon * (Real(1.0) - p.deltN);
                      n_flux_sloppy_sldon = n_flux_egest * qdon * p.deltN;
                      bio(i, j, k - 1, sc_don) =
                        bio(i, j, k - 1, sc_don) + n_flux_sloppy_sldon;
                      bio(i, j, k - 1, sc_nh4) =
                        bio(i, j, k - 1, sc_nh4) + n_flux_sloppy_ldon;
                      bio(i, j, k - 1, sc_talk) =
                        bio(i, j, k - 1, sc_talk) + n_flux_sloppy_ldon;
                      n_flux_pmortal =
                        cff2 *
                        amrex::max(
                          bio(i, j, k - 1, sc_phyt) - p.PhyMin, Real(0.0));
                      bio(i, j, k - 1, sc_phyt) =
                        bio(i, j, k - 1, sc_phyt) - n_flux_pmortal;
                      bio(i, j, k - 1, sc_chlo) =
                        bio(i, j, k - 1, sc_chlo) -
                        cff2 *
                          amrex::max(
                            bio(i, j, k - 1, sc_chlo) - p.ChlMin, Real(0.0));
                      bio(i, j, k - 1, sc_sden) =
                        bio(i, j, k - 1, sc_sden) + n_flux_pmortal;
                      qdoc = qdon;
                      bio(i, j, k - 1, sc_ldec) =
                        bio(i, j, k - 1, sc_ldec) +
                        p.PhyCN * n_flux_egest * (Real(1.0) - qdoc);
                      bio(i, j, k - 1, sc_tic) =
                        bio(i, j, k - 1, sc_tic) +
                        p.PhyCN * n_flux_egest * qdoc * (Real(1.0) - p.deltC);
                      bio(i, j, k - 1, sc_doc) =
                        bio(i, j, k - 1, sc_doc) +
                        p.PhyCN * n_flux_egest * qdoc * p.deltC;
                      bio(i, j, k - 1, sc_sdec) =
                        bio(i, j, k - 1, sc_sdec) + p.PhyCN * n_flux_pmortal;
                      bio(i, j, k - 1, sc_tic) =
                        bio(i, j, k - 1, sc_tic) +
                        (p.PhyCN - p.ZooCN) * n_flux_assim;
                      bio(i, j, k - 1, sc_oxyg) =
                        bio(i, j, k - 1, sc_oxyg) - roxnh4 * n_flux_sloppy_ldon;
                      bio(i, j, k - 1, sc_oxyg) =
                        amrex::max(Real(0.0), bio(i, j, k - 1, sc_oxyg));
                  }
                  // Zooplankton mortality, excretion and basal metabolism.
                  cff1 = dtdays * p.ZooBM;
                  fac2 = dtdays * p.ZooMR;
                  fac3 = dtdays * p.ZooER;
                  for (k = 1; k <= nz; k += 1) {
                      fac1 = fac3 * bio(i, j, k - 1, sc_phyt) *
                             bio(i, j, k - 1, sc_phyt) /
                             (p.K_Phy + bio(i, j, k - 1, sc_phyt) *
                                          bio(i, j, k - 1, sc_phyt));
                      cff2 = fac2 * bio(i, j, k - 1, sc_zoop);
                      cff3 = fac1 * p.ZooAE_N;
                      bio(i, j, k - 1, sc_zoop) =
                        bio(i, j, k - 1, sc_zoop) / (Real(1.0) + cff2 + cff3);
                      n_flux_zmortal = cff2 * bio(i, j, k - 1, sc_zoop);
                      n_flux_zexcret = cff3 * bio(i, j, k - 1, sc_zoop);
                      bio(i, j, k - 1, sc_nh4) =
                        bio(i, j, k - 1, sc_nh4) + n_flux_zexcret;
                      bio(i, j, k - 1, sc_po4) =
                        bio(i, j, k - 1, sc_po4) + p.R_P2N * n_flux_zexcret;
                      bio(i, j, k - 1, sc_sden) =
                        bio(i, j, k - 1, sc_sden) + n_flux_zmortal;
                      n_flux_zmetabo =
                        cff1 *
                        amrex::max(
                          bio(i, j, k - 1, sc_zoop) - p.ZooMin, Real(0.0));
                      bio(i, j, k - 1, sc_zoop) =
                        bio(i, j, k - 1, sc_zoop) - n_flux_zmetabo;
                      bio(i, j, k - 1, sc_nh4) =
                        bio(i, j, k - 1, sc_nh4) + n_flux_zmetabo;
                      bio(i, j, k - 1, sc_po4) =
                        bio(i, j, k - 1, sc_po4) + p.R_P2N * n_flux_zmetabo;
                      bio(i, j, k - 1, sc_oxyg) =
                        bio(i, j, k - 1, sc_oxyg) -
                        roxnh4 * (n_flux_zmetabo + n_flux_zexcret);
                      bio(i, j, k - 1, sc_oxyg) =
                        amrex::max(Real(0.), bio(i, j, k - 1, sc_oxyg));
                      bio(i, j, k - 1, sc_sdec) =
                        bio(i, j, k - 1, sc_sdec) + p.ZooCN * n_flux_zmortal;
                      bio(i, j, k - 1, sc_tic) =
                        bio(i, j, k - 1, sc_tic) +
                        p.ZooCN * (n_flux_zmetabo + n_flux_zexcret);
                      bio(i, j, k - 1, sc_talk) = bio(i, j, k - 1, sc_talk) +
                                                  n_flux_zmetabo +
                                                  n_flux_zexcret;
                  }
                  // Coagulation into large detritus.
                  fac1 = dtdays * p.CoagR;
                  for (k = 1; k <= nz; k += 1) {
                      cff1 = fac1 * (bio(i, j, k - 1, sc_sden) +
                                     bio(i, j, k - 1, sc_phyt));
                      cff2 = Real(1.0) / (Real(1.0) + cff1);
                      bio(i, j, k - 1, sc_phyt) =
                        bio(i, j, k - 1, sc_phyt) * cff2;
                      bio(i, j, k - 1, sc_chlo) =
                        bio(i, j, k - 1, sc_chlo) * cff2;
                      bio(i, j, k - 1, sc_sden) =
                        bio(i, j, k - 1, sc_sden) * cff2;
                      bio(i, j, k - 1, sc_sdec) =
                        bio(i, j, k - 1, sc_sdec) * cff2;
                      n_flux_coagp = bio(i, j, k - 1, sc_phyt) * cff1;
                      n_flux_coagd = bio(i, j, k - 1, sc_sden) * cff1;
                      bio(i, j, k - 1, sc_lden) =
                        bio(i, j, k - 1, sc_lden) + n_flux_coagp + n_flux_coagd;
                      bio(i, j, k - 1, sc_ldec) =
                        bio(i, j, k - 1, sc_ldec) + p.PhyCN * n_flux_coagp +
                        bio(i, j, k - 1, sc_sdec) * cff1;
                  }
                  for (k = 1; k <= nz; k += 1) {
                      f_ntr = bio(i, j, k - 1, sc_oxyg) /
                              (bio(i, j, k - 1, sc_oxyg) + k_ntr);
                      f_dnf = k_dnf / (bio(i, j, k - 1, sc_oxyg) + k_dnf) *
                              bio(i, j, k - 1, sc_no3) /
                              (bio(i, j, k - 1, sc_no3) + k_wno3);
                      cff = dtdays *
                            std::exp(Real(0.0875) * bio(i, j, k - 1, sc_temp));
                      cff1 = cff * p.SDeRRN;
                      cff2 = cff1 * p.deltN;
                      cff3 = cff1 * (Real(1.0) - p.deltN);
                      cff4 = cff * p.LDeRRN;
                      cff5 = cff4 * p.deltN;
                      cff6 = cff4 * (Real(1.0) - p.deltN);
                      rdon = cff * p.RDeRRN;
                      fac = (f_ntr + f_dnf);
                      bio(i, j, k - 1, sc_sden) =
                        bio(i, j, k - 1, sc_sden) /
                        (Real(1.0) + cff2 + cff3 * fac);
                      bio(i, j, k - 1, sc_lden) =
                        bio(i, j, k - 1, sc_lden) /
                        (Real(1.0) + cff5 + cff6 * fac);
                      bio(i, j, k - 1, sc_don) =
                        bio(i, j, k - 1, sc_don) / (Real(1.0) + rdon * fac);
                      n_flux_remines = bio(i, j, k - 1, sc_sden) * cff3 * fac;
                      n_flux_reminel = bio(i, j, k - 1, sc_lden) * cff6 * fac;
                      n_flux_reminer = bio(i, j, k - 1, sc_don) * rdon * fac;
                      bio(i, j, k - 1, sc_nh4) =
                        bio(i, j, k - 1, sc_nh4) + n_flux_remines +
                        n_flux_reminel + n_flux_reminer;
                      bio(i, j, k - 1, sc_po4) =
                        bio(i, j, k - 1, sc_po4) +
                        p.R_P2N *
                          (n_flux_remines + n_flux_reminel + n_flux_reminer);
                      cff7 = (bio(i, j, k - 1, sc_don) * rdon +
                              bio(i, j, k - 1, sc_sden) * cff3 +
                              bio(i, j, k - 1, sc_lden) * cff6) *
                             roxnh4 / (bio(i, j, k - 1, sc_oxyg) + k_ntr);
                      bio(i, j, k - 1, sc_oxyg) =
                        bio(i, j, k - 1, sc_oxyg) / (Real(1.0) + cff7);
                      bio(i, j, k - 1, sc_talk) =
                        bio(i, j, k - 1, sc_talk) + n_flux_remines +
                        n_flux_reminel + n_flux_reminer;
                      cff8 = (bio(i, j, k - 1, sc_sden) * cff3 +
                              bio(i, j, k - 1, sc_lden) * cff6 +
                              bio(i, j, k - 1, sc_don) * rdon) *
                             Real(5.3) / (bio(i, j, k - 1, sc_oxyg) + k_dnf) *
                             k_dnf / (bio(i, j, k - 1, sc_no3) + k_wno3);
                      bio(i, j, k - 1, sc_no3) =
                        bio(i, j, k - 1, sc_no3) / (Real(1.0) + cff8);
                      bio(i, j, k - 1, sc_talk) =
                        bio(i, j, k - 1, sc_talk) +
                        bio(i, j, k - 1, sc_no3) * cff8;
                      bio(i, j, k - 1, sc_don) =
                        bio(i, j, k - 1, sc_don) +
                        bio(i, j, k - 1, sc_sden) * cff2 +
                        bio(i, j, k - 1, sc_lden) * cff5;
                  }
                  // Air-sea exchange uses the reference Schmidt/rate
                  // coefficient pair.
                  cff1 = rho0 * Real(550.0);
                  cff2 = dtdays * Real(0.251) * Real(24.0) / Real(100.0);
                  k = nz;
                  u10squ =
                    bulk ? uwind(i, j, 0) * uwind(i, j, 0) +
                             vwind(i, j, 0) * vwind(i, j, 0)
                         : cff1 * std::sqrt(
                                    std::pow(
                                      (Real(0.5) *
                                       (sustr(i, j, 0) + sustr(i + 1, j, 0))),
                                      2) +
                                    std::pow(
                                      (Real(0.5) *
                                       (svstr(i, j, 0) + svstr(i, j + 1, 0))),
                                      2));
                  schmidtn_ox =
                    a_o2 -
                    bio(i, j, k - 1, sc_temp) *
                      (b_o2 -
                       bio(i, j, k - 1, sc_temp) *
                         (c_o2 - bio(i, j, k - 1, sc_temp) *
                                   (d_o2 - bio(i, j, k - 1, sc_temp) * e_o2)));
                  cff3 = cff2 * u10squ * std::sqrt(Real(660.0) / schmidtn_ox);
                  ts = std::log(
                    (Real(298.15) - bio(i, j, k - 1, sc_temp)) /
                    (Real(273.15) + bio(i, j, k - 1, sc_temp)));
                  aa =
                    oa0 +
                    ts *
                      (oa1 + ts * (oa2 + ts * (oa3 + ts * (oa4 + ts * oa5)))) +
                    bio(i, j, k - 1, sc_salt) *
                      (ob0 + ts * (ob1 + ts * (ob2 + ts * ob3))) +
                    oc0 * bio(i, j, k - 1, sc_salt) * bio(i, j, k - 1, sc_salt);
                  o2satu = l2mol * std::exp(aa);
                  o2_flux = cff3 * (o2satu - bio(i, j, k - 1, sc_oxyg));
                  bio(i, j, k - 1, sc_oxyg) =
                    bio(i, j, k - 1, sc_oxyg) + o2_flux * hz_inv(i, j, k);
                  for (k = 1; k <= nz; k += 1) {
                      f_ntr = bio(i, j, k - 1, sc_oxyg) /
                              (bio(i, j, k - 1, sc_oxyg) + k_ntr);
                      f_dnf = k_dnf / (bio(i, j, k - 1, sc_oxyg) + k_dnf) *
                              bio(i, j, k - 1, sc_no3) /
                              (bio(i, j, k - 1, sc_no3) + k_wno3);
                      cff = dtdays *
                            std::exp(Real(0.0875) * bio(i, j, k - 1, sc_temp));
                      cff1 = cff * p.SDeRRC;
                      cff2 = cff1 * p.deltC;
                      cff3 = cff1 * (Real(1.0) - p.deltC);
                      cff4 = cff * p.LDeRRC;
                      cff5 = cff4 * p.deltC;
                      cff6 = cff4 * (Real(1.0) - p.deltC);
                      rdoc = cff * p.RDeRRC;
                      fac = (f_ntr + f_dnf);
                      bio(i, j, k - 1, sc_sdec) =
                        bio(i, j, k - 1, sc_sdec) /
                        (Real(1.0) + cff2 + cff3 * fac);
                      bio(i, j, k - 1, sc_ldec) =
                        bio(i, j, k - 1, sc_ldec) /
                        (Real(1.0) + cff5 + cff6 * fac);
                      bio(i, j, k - 1, sc_doc) =
                        bio(i, j, k - 1, sc_doc) / (Real(1.0) + rdoc * fac);
                      c_flux_remines = bio(i, j, k - 1, sc_sdec) * cff3 * fac;
                      c_flux_reminel = bio(i, j, k - 1, sc_ldec) * cff6 * fac;
                      c_flux_reminer = bio(i, j, k - 1, sc_doc) * rdoc * fac;
                      bio(i, j, k - 1, sc_tic) =
                        bio(i, j, k - 1, sc_tic) + c_flux_remines +
                        c_flux_reminel + c_flux_reminer;
                      bio(i, j, k - 1, sc_doc) =
                        bio(i, j, k - 1, sc_doc) +
                        cff2 * bio(i, j, k - 1, sc_sdec) +
                        cff5 * bio(i, j, k - 1, sc_ldec);
                  }
                  k = nz;
                  pco2 = REMORABiology::pco2_water(
                    bio(i, j, k - 1, sc_temp), bio(i, j, k - 1, sc_salt),
                    bio(i, j, k - 1, sc_tic), bio(i, j, k - 1, sc_talk));
                  // The remineralization loop reuses cff2, so restore the rate
                  // coefficient; u10squ from the O2 exchange is still current.
                  cff2 = dtdays * Real(0.251) * Real(24.0) / Real(100.0);
                  schmidtn =
                    a_co2 -
                    bio(i, j, k - 1, sc_temp) *
                      (b_co2 - bio(i, j, k - 1, sc_temp) *
                                 (c_co2 - bio(i, j, k - 1, sc_temp) *
                                            (d_co2 - bio(i, j, k - 1, sc_temp) *
                                                       e_co2)));
                  cff3 = cff2 * u10squ * std::sqrt(Real(660.0) / schmidtn);
                  tempk =
                    Real(0.01) * (bio(i, j, k - 1, sc_temp) + Real(273.15));
                  co2_sol = std::exp(
                    a1 + a2 / tempk + a3 * std::log(tempk) +
                    bio(i, j, k - 1, sc_salt) *
                      (b1 + tempk * (b2 + b3 * tempk)));
                  co2_flux = cff3 * co2_sol * (p.pCO2air - pco2);
                  bio(i, j, k - 1, sc_tic) =
                    bio(i, j, k - 1, sc_tic) + co2_flux * hz_inv(i, j, k);
                  // PPM/WENO sinking with the reference monotone boundary
                  // closures.
                  for (isink = 1; isink <= nsink; isink += 1) {
                      ibio = idsink[isink - 1];
                      for (k = 1; k <= nz; k += 1) {
                          qc(i, j, k) =
                            amrex::max(bio(i, j, k - 1, ibio), Real(0.));
                      }
                      for (k = nz - 1; k >= 1; k += -1) {
                          fc(i, j, k) =
                            (qc(i, j, k + 1) - qc(i, j, k)) * hz_inv2(i, j, k);
                      }
                      for (k = 2; k <= nz - 1; k += 1) {
                          dltr = hz(i, j, k - 1) * fc(i, j, k);
                          dltl = hz(i, j, k - 1) * fc(i, j, k - 1);
                          cff = hz(i, j, k - 1 - 1) +
                                Real(2.0) * hz(i, j, k - 1) +
                                hz(i, j, k + 1 - 1);
                          cffr = cff * fc(i, j, k);
                          cffl = cff * fc(i, j, k - 1);
                          if ((dltr * dltl) <= Real(0.0)) {
                              dltr = Real(0.0);
                              dltl = Real(0.0);
                          } else if (std::abs(dltr) > std::abs(cffl)) {
                              dltr = cffl;
                          } else if (std::abs(dltl) > std::abs(cffr)) {
                              dltl = cffr;
                          }
                          cff = (dltr - dltl) * hz_inv3(i, j, k);
                          dltr = dltr - cff * hz(i, j, k + 1 - 1);
                          dltl = dltl + cff * hz(i, j, k - 1 - 1);
                          br(i, j, k) = qc(i, j, k) + dltr;
                          bl(i, j, k) = qc(i, j, k) - dltl;
                          wr(i, j, k) = std::pow((Real(2.0) * dltr - dltl), 2);
                          wl(i, j, k) = std::pow((dltr - Real(2.0) * dltl), 2);
                      }
                      cff = Real(1.0e-14);
                      for (k = 2; k <= nz - 2; k += 1) {
                          dltl = amrex::max(cff, wl(i, j, k));
                          dltr = amrex::max(cff, wr(i, j, k + 1));
                          br(i, j, k) =
                            (dltr * br(i, j, k) + dltl * bl(i, j, k + 1)) /
                            (dltr + dltl);
                          bl(i, j, k + 1) = br(i, j, k);
                      }
                      fc(i, j, nz) = Real(0.0);
                      br(i, j, nz) = qc(i, j, nz);
                      bl(i, j, nz) = qc(i, j, nz);
                      br(i, j, nz - 1) = qc(i, j, nz);
                      bl(i, j, 2) = qc(i, j, 1);
                      br(i, j, 1) = qc(i, j, 1);
                      bl(i, j, 1) = qc(i, j, 1);
                      for (k = 1; k <= nz; k += 1) {
                          dltr = br(i, j, k) - qc(i, j, k);
                          dltl = qc(i, j, k) - bl(i, j, k);
                          cffr = Real(2.0) * dltr;
                          cffl = Real(2.0) * dltl;
                          if ((dltr * dltl) < Real(0.0)) {
                              dltr = Real(0.0);
                              dltl = Real(0.0);
                          } else if (std::abs(dltr) > std::abs(cffl)) {
                              dltr = cffl;
                          } else if (std::abs(dltl) > std::abs(cffr)) {
                              dltl = cffr;
                          }
                          br(i, j, k) = qc(i, j, k) + dltr;
                          bl(i, j, k) = qc(i, j, k) - dltl;
                      }
                      cff = dtdays * std::abs(wbio[isink - 1]);
                      for (k = 1; k <= nz; k += 1) {
                          fc(i, j, k - 1) = Real(0.0);
                          wl(i, j, k) = z_w(i, j, k - 1) + cff;
                          wr(i, j, k) = hz(i, j, k - 1) * qc(i, j, k);
                          ksource(i, j, k) = k;
                      }
                      for (k = 1; k <= nz; k += 1) {
                          for (ks = k; ks <= nz - 1; ks += 1) {
                              if (wl(i, j, k) > z_w(i, j, ks)) {
                                  ksource(i, j, k) = ks + 1;
                                  fc(i, j, k - 1) =
                                    fc(i, j, k - 1) + wr(i, j, ks);
                              }
                          }
                      }
                      for (k = 1; k <= nz; k += 1) {
                          ks = ksource(i, j, k);
                          cu = amrex::min(
                            Real(1.0), (wl(i, j, k) - z_w(i, j, ks - 1)) *
                                         hz_inv(i, j, ks));
                          fc(i, j, k - 1) =
                            fc(i, j, k - 1) +
                            hz(i, j, ks - 1) * cu *
                              (bl(i, j, ks) +
                               cu * (Real(0.5) * (br(i, j, ks) - bl(i, j, ks)) -
                                     (Real(1.5) - cu) *
                                       (br(i, j, ks) + bl(i, j, ks) -
                                        Real(2.0) * qc(i, j, ks))));
                      }
                      for (k = 1; k <= nz; k += 1) {
                          bio(i, j, k - 1, ibio) =
                            qc(i, j, k) +
                            (fc(i, j, k) - fc(i, j, k - 1)) * hz_inv(i, j, k);
                      }
                      // Stress-dependent resuspension and instantaneous
                      // sediment recycling/burial.
                      cff2 = Real(4.0) / Real(16.0);
                      cff3 = Real(115.0) / Real(16.0);
                      cff4 = Real(106.0) / Real(16.0);
                      if (
                        (ibio == sc_phyt) || (ibio == sc_sdec) ||
                        (ibio == sc_ldec) || (ibio == sc_sden) ||
                        (ibio == sc_lden)) {
                          ustarb = std::sqrt(std::sqrt(
                            std::pow(
                              (Real(0.5) *
                               (bustr(i, j, 0) + bustr(i + 1, j, 0))),
                              2) +
                            std::pow(
                              (Real(0.5) *
                               (bvstr(i, j, 0) + bvstr(i, j + 1, 0))),
                              2)));
                          resuspr = amrex::min(
                            std::pow(ustarb, 2) / Real(9.61e-6), Real(1.0));
                          cff1 = fc(i, j, 0) * hz_inv(i, j, 1) *
                                 (Real(1.0) - resuspr);
                          cff5 = fc(i, j, 0) * hz_inv(i, j, 1) * resuspr;
                          if (
                            (ibio == sc_phyt) || (ibio == sc_sden) ||
                            (ibio == sc_lden)) {
                              bio(i, j, 0, sc_sden) =
                                bio(i, j, 0, sc_sden) + cff5;
                          }
                          if ((ibio == sc_sdec) || (ibio == sc_ldec)) {
                              bio(i, j, 0, sc_sdec) =
                                bio(i, j, 0, sc_sdec) + cff5;
                          }
                          if (ibio == sc_phyt) {
                              bio(i, j, 0, sc_sdec) =
                                bio(i, j, 0, sc_sdec) + cff5 * p.PhyCN;
                          }
                      }
                      if (
                        (ibio == sc_phyt) || (ibio == sc_sden) ||
                        (ibio == sc_lden)) {
                          if (ibio == sc_phyt) {
                              cnbur = p.PhyCN;
                          } else {
                              cnbur = Real(9.3);
                          }
                          cff5 = Real(0.01);
                          fac1 = cnbur * Real(12.0) / Real(1000.0) * cff1 *
                                 hz(i, j, 1 - 1) * (Real(365.0) / dtdays);
                          cbe = amrex::min(
                            Real(0.50),
                            Real(0.023) * std::pow(fac1, Real(0.5797)));
                          cff1 = (Real(1.0) - cbe) * cff1;
                          ts = std::log(
                            (Real(298.15) - bio(i, j, 0, sc_temp)) /
                            (Real(273.15) + bio(i, j, 0, sc_temp)));
                          aa =
                            oa0 +
                            ts * (oa1 +
                                  ts * (oa2 +
                                        ts * (oa3 + ts * (oa4 + ts * oa5)))) +
                            bio(i, j, 0, sc_salt) *
                              (ob0 + ts * (ob1 + ts * (ob2 + ts * ob3))) +
                            oc0 * bio(i, j, 0, sc_salt) * bio(i, j, 0, sc_salt);
                          o2satb = l2mol * std::exp(aa);
                          lbo2 = (k_bo2 * (o2satb - bio(i, j, 0, sc_oxyg))) /
                                 (o2satb * (bio(i, j, 0, sc_oxyg) + k_bo2));
                          fac = cff1 * (Real(1.0) + Real(3.0) * lbo2);
                          bio(i, j, 0, sc_nh4) =
                            bio(i, j, 0, sc_nh4) + cff2 * fac;
                          bio(i, j, 0, sc_talk) =
                            bio(i, j, 0, sc_talk) + cff2 * fac;
                          bio(i, j, 0, sc_don) =
                            bio(i, j, 0, sc_don) + cff5 * fac;
                          bio(i, j, 0, sc_oxyg) =
                            bio(i, j, 0, sc_oxyg) -
                            cff1 * cff3 * (Real(1.0) - lbo2);
                      }
                      if (
                        (ibio == sc_phyt) || (ibio == sc_sdec) ||
                        (ibio == sc_ldec)) {
                          if (ibio == sc_phyt) {
                              cff1 = cff1 / (Real(1.0) - cbe);
                          }
                          fac1 = Real(12.0) / Real(1000.0) * cff1 *
                                 hz(i, j, 1 - 1) * (Real(365.0) / dtdays);
                          if (ibio == sc_phyt) {
                              fac1 = fac1 * p.PhyCN;
                          }
                          cbe = amrex::min(
                            Real(0.50),
                            Real(0.023) * std::pow(fac1, Real(0.5797)));
                          if (ibio == sc_phyt) {
                              fac2 = p.PhyCN;
                          } else {
                              fac2 = Real(1.0);
                          }
                          cff1 = (Real(1.0) - cbe) * cff1;
                          cff5 = Real(0.01);
                          fac = fac2 * cff1;
                          bio(i, j, 0, sc_doc) =
                            bio(i, j, 0, sc_doc) + fac * cff5;
                          bio(i, j, 0, sc_tic) =
                            bio(i, j, 0, sc_tic) + fac * (Real(1.0) - cff5);
                      }
                  }
              }
              for (k = 1; k <= nz; k += 1) {
                  bio(i, j, k - 1, sc_tic) =
                    amrex::min(bio(i, j, k - 1, sc_tic), Real(3000.0));
                  bio(i, j, k - 1, sc_tic) =
                    amrex::max(bio(i, j, k - 1, sc_tic), Real(1.0));
              }
              for (itrc = 1; itrc <= 18; itrc += 1) {
                  ibio = (itrc + 1);
                  for (k = 1; k <= nz; k += 1) {
                      cff = bio(i, j, k - 1, ibio) - bio_old(i, j, k - 1, ibio);
                      cff = cff * mask(i, j, 0);
                      next(i, j, k - 1, state_comp[ibio]) =
                        next(i, j, k - 1, state_comp[ibio]) +
                        cff * hz(i, j, k - 1);
                  }
              }
          });
    }
}
} // namespace REMORABiology
