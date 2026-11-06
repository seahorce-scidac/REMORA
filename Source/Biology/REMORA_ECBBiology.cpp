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
// Native port of ROMS ECB (ROMS/Nonlinear/Biology/ecb.h); see
// Docs/sphinx_doc/ECB.rst for the pinned source and option set.
#include <REMORA_ECBBiology.H>
#include <REMORA_BiologyCarbonate.H>
#include <REMORA_BiologyGasExchange.H>
#include <REMORA_DataStruct.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParmParse.H>
#include <cmath>
#include <iterator>

using namespace amrex;

namespace {
// State tracer names, in the order REMORA lays them out. The kernel looks each
// one up by name, so this order is only the plotfile and checkpoint layout.
const char* const ecb_tracer_names[] = {
    "NO3", "NH4", "chlorophyll", "phytoplankton", "zooplankton",
    "LdetritusN", "SdetritusN", "DON", "PO4", "LdetritusC",
    "SdetritusC", "TIC", "alkalinity", "DOC", "oxygen",
    "Rdon", "Rdoc", "calcium"};
}

void
REMORABiology::ECBParameters::init_params (const std::string& remora_prefix)
{
    ParmParse pp(remora_prefix + ".ecb");

    // ECB is ported with one ROMS option set; reject the Fennel-style switches
    // rather than ignore them.
    for (const char* option : {"carbon", "oxygen", "bio_sediment", "denitrification", "po4",
                               "talk_nonconserv", "co2_schmidt", "oxygen_schmidt",
                               "pco2air_type", "pco2_rz", "sediment_proxy",
                               "net_sulfate_reduction", "wet_dry"}) {
        if (pp.contains(option)) {
            amrex::Abort("ECB implements a fixed ROMS option set; unsupported option: " +
                         remora_prefix + ".ecb." + option);
        }
    }

    // The ROMS Fortran kernel and the parity diagnostics exist only for Fennel,
    // so asking for either with ECB is an error rather than a silent C++ run.
    {
        ParmParse pp_remora(remora_prefix);
        int use_cpp = 1;
        int debug = 0;
        pp_remora.query("use_biology_cpp_answer", use_cpp);
        pp_remora.query("biology_debug", debug);
        if (use_cpp != 1) {
            amrex::Abort(remora_prefix + ".use_biology_cpp_answer = 0 selects the ROMS Fortran "
                         "kernel, which exists only for Fennel, not ECB.");
        }
        if (debug > 0) {
            amrex::Abort(remora_prefix + ".biology_debug > 0 requests the parity diagnostics, "
                         "which exist only for Fennel, not ECB.");
        }
    }

    pp.queryAdd("BioIter", BioIter);
    pp.queryAdd("AttSW", AttSW);
    pp.queryAdd("AttChl", AttChl);
    pp.queryAdd("PARfrac", PARfrac);
    pp.queryAdd("Vp0", Vp0);
    pp.queryAdd("I_thNH4", I_thNH4);
    pp.queryAdd("D_p5NH4", D_p5NH4);
    pp.queryAdd("NitriR", NitriR);
    pp.queryAdd("K_NO3", K_NO3);
    pp.queryAdd("K_NH4", K_NH4);
    pp.queryAdd("K_PO4", K_PO4);
    pp.queryAdd("K_Phy", K_Phy);
    pp.queryAdd("Chl2C_m", Chl2C_m);
    pp.queryAdd("ChlMin", ChlMin);
    pp.queryAdd("PhyCN", PhyCN);
    pp.queryAdd("R_P2N", R_P2N);
    pp.queryAdd("PhyIS", PhyIS);
    pp.queryAdd("PhyMin", PhyMin);
    pp.queryAdd("PhyMR", PhyMR);
    pp.queryAdd("ElDON", ElDON);
    pp.queryAdd("EsDON", EsDON);
    pp.queryAdd("deltN", deltN);
    pp.queryAdd("deltC", deltC);
    pp.queryAdd("gammaC", gammaC);
    pp.queryAdd("slCexc", slCexc);
    pp.queryAdd("ZooAE_N", ZooAE_N);
    pp.queryAdd("ZooBM", ZooBM);
    pp.queryAdd("ZooCN", ZooCN);
    pp.queryAdd("ZooER", ZooER);
    pp.queryAdd("ZooGR", ZooGR);
    pp.queryAdd("ZooMin", ZooMin);
    pp.queryAdd("ZooMR", ZooMR);
    pp.queryAdd("LDeRRN", LDeRRN);
    pp.queryAdd("LDeRRC", LDeRRC);
    pp.queryAdd("CoagR", CoagR);
    pp.queryAdd("SDeRRN", SDeRRN);
    pp.queryAdd("SDeRRC", SDeRRC);
    pp.queryAdd("RDeRRN", RDeRRN);
    pp.queryAdd("RDeRRC", RDeRRC);
    pp.queryAdd("wPhy", wPhy);
    pp.queryAdd("wLDet", wLDet);
    pp.queryAdd("wSDet", wSDet);
    pp.queryAdd("pCO2air", pCO2air);

    // Reject values that divide by zero in the kernel, where they would
    // otherwise surface as NaN in the tracers with no diagnostic.
    auto require = [&] (bool ok, const char* name, const char* condition, Real value) {
        if (!ok) {
            amrex::Abort(remora_prefix + ".ecb." + name + " must be " + condition +
                         "; got " + std::to_string(value));
        }
    };
    // Light attenuation (AttSW + AttChl*chl + DOC term) * dz may be zero, which the
    // kernel handles, but a negative coefficient would amplify light with depth.
    require(AttSW >= Real(0.0), "AttSW", "non-negative", AttSW);
    require(AttChl >= Real(0.0), "AttChl", "non-negative", AttChl);
    // Carbon excess uptake divides by slCexc.
    require(slCexc > Real(0.0), "slCexc", "positive", slCexc);
    // Grazing divides by K_Phy + phyt^2, which is K_Phy when phyt is zero.
    require(K_Phy > Real(0.0), "K_Phy", "positive", K_Phy);
    // A non-positive atmospheric pCO2 reverses the sign of the air-sea CO2 flux.
    require(pCO2air > Real(0.0), "pCO2air", "positive", pCO2air);
}

void
REMORABiology::ECBBiology::configure (const std::string& remora_prefix)
{
    m_parameters.init_params(remora_prefix);
    m_declaration.tracer_names.assign(std::begin(ecb_tracer_names), std::end(ecb_tracer_names));
    m_declaration.plot_aliases = {"ecb"};
}

void
REMORABiology::ECBBiology::advance (const BiologyAdvanceContext& context) const
{
    const auto parms = m_parameters;
    if (parms.BioIter <= 0) {
        return;
    }

    auto const& mf_cons_old = context.old_state;
    auto& mf_cons_new = context.new_state;
    const bool do_bulk_flux = context.solver_choice.bulk_fluxes;
    const int N = context.top_layer + 1;   // number of vertical layers
    const Real dt_lev = context.dt;
    const Real rho0 = context.solver_choice.rho0;

    if (context.hz == nullptr || context.z_w == nullptr || context.mask == nullptr ||
        context.shortwave == nullptr ||
        context.bottom_stress_u == nullptr || context.bottom_stress_v == nullptr ||
        (do_bulk_flux && (context.wind_u == nullptr || context.wind_v == nullptr)) ||
        (!do_bulk_flux && (context.stress_u == nullptr || context.stress_v == nullptr))) {
        amrex::Abort("ECB biology requires Hz, z_w, rmask, srflx, bottom stresses, and "
                     "either uwind/vwind (bulk_fluxes) or sustr/svstr");
    }
    if (N < 2) {
        amrex::Abort("ECB sinking requires at least two vertical layers");
    }
    if (dt_lev <= Real(0.0)) {
        amrex::Abort("ECB requires a positive timestep when BioIter > 0");
    }

    // Set time-stepping according to the number of iterations.
    const Real dtdays = dt_lev * Real(1.0 / 86400.0) / Real(parms.BioIter);

    // Air-sea exchange always uses the Wanninkhof (2014) relations (ROMS
    // RW14_OXYGEN_SC and RW14_CO2_SC).
    const SchmidtCoefficients o2_schmidt = GasExchange::O2Wanninkhof2014;
    const SchmidtCoefficients co2_schmidt = GasExchange::CO2Wanninkhof2014;

    // Set biological tracer component identifiers, looked up by name.
    auto const& decl = m_declaration;
    auto comp = [&] (const char* name) {
        const int c = decl.component_index(name, context.component_start);
        AMREX_ALWAYS_ASSERT(c >= 0);
        return c;
    };
    const int c_no3 = comp("NO3");
    const int c_nh4 = comp("NH4");
    const int c_chlo = comp("chlorophyll");
    const int c_phyt = comp("phytoplankton");
    const int c_zoop = comp("zooplankton");
    const int c_lden = comp("LdetritusN");
    const int c_sden = comp("SdetritusN");
    const int c_don = comp("DON");
    const int c_po4 = comp("PO4");
    const int c_ldec = comp("LdetritusC");
    const int c_sdec = comp("SdetritusC");
    const int c_tic = comp("TIC");
    const int c_talk = comp("alkalinity");
    const int c_doc = comp("DOC");
    const int c_oxyg = comp("oxygen");
    const int c_rdon = comp("Rdon");
    const int c_rdoc = comp("Rdoc");
    const int c_calc = comp("calcium");

    // Scratch arrays mirror the ROMS Bio, qc, bR, bL, WR, WL, and
    // inverse-thickness work arrays, but use per-tile FArrayBox storage.
    int ncell = 0;
    const int sc_no3 = ncell++;
    const int sc_nh4 = ncell++;
    const int sc_chlo = ncell++;
    const int sc_phyt = ncell++;
    const int sc_zoop = ncell++;
    const int sc_lden = ncell++;
    const int sc_sden = ncell++;
    const int sc_don = ncell++;
    const int sc_po4 = ncell++;
    const int sc_ldec = ncell++;
    const int sc_sdec = ncell++;
    const int sc_tic = ncell++;
    const int sc_talk = ncell++;
    const int sc_doc = ncell++;
    const int sc_oxyg = ncell++;
    const int sc_rdon = ncell++;
    const int sc_rdoc = ncell++;
    const int sc_calc = ncell++;
    const int sc_temp = ncell++;
    const int sc_salt = ncell++;
    const int sc_inv_hz = ncell++;
    const int sc_inv_hz2 = ncell++;
    const int sc_inv_hz3 = ncell++;
    const int sc_qc = ncell++;
    const int sc_bR = ncell++;
    const int sc_bL = ncell++;
    const int sc_WR = ncell++;
    const int sc_WL = ncell++;

    int nw = 0;
    const int sw_FC = nw++;

    for (MFIter mfi(mf_cons_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box bx = mfi.tilebox();
        Box bx2d = bx;
        bx2d.makeSlab(2, 0);
        Box wbx = surroundingNodes(bx, 2);

        FArrayBox fab_cell(bx, ncell, amrex::The_Async_Arena());
        fab_cell.template setVal<RunOn::Device>(Real(0.0));
        FArrayBox fab_w(wbx, nw, amrex::The_Async_Arena());
        fab_w.template setVal<RunOn::Device>(Real(0.0));

        auto no3 = fab_cell.array(sc_no3);
        auto nh4 = fab_cell.array(sc_nh4);
        auto chlo = fab_cell.array(sc_chlo);
        auto phyt = fab_cell.array(sc_phyt);
        auto zoop = fab_cell.array(sc_zoop);
        auto lden = fab_cell.array(sc_lden);
        auto sden = fab_cell.array(sc_sden);
        auto don = fab_cell.array(sc_don);
        auto po4 = fab_cell.array(sc_po4);
        auto ldec = fab_cell.array(sc_ldec);
        auto sdec = fab_cell.array(sc_sdec);
        auto tic = fab_cell.array(sc_tic);
        auto talk = fab_cell.array(sc_talk);
        auto doc = fab_cell.array(sc_doc);
        auto oxyg = fab_cell.array(sc_oxyg);
        auto rdon_s = fab_cell.array(sc_rdon);
        auto rdoc_s = fab_cell.array(sc_rdoc);
        auto calc = fab_cell.array(sc_calc);
        auto temp = fab_cell.array(sc_temp);
        auto salt = fab_cell.array(sc_salt);
        auto inv_hz = fab_cell.array(sc_inv_hz);
        auto inv_hz2 = fab_cell.array(sc_inv_hz2);
        auto inv_hz3 = fab_cell.array(sc_inv_hz3);
        auto qc = fab_cell.array(sc_qc);
        auto bR = fab_cell.array(sc_bR);
        auto bL = fab_cell.array(sc_bL);
        auto WR = fab_cell.array(sc_WR);
        auto WL = fab_cell.array(sc_WL);
        auto FC = fab_w.array(sw_FC);

        Array4<Real const> const& state_old = mf_cons_old.const_array(mfi);
        Array4<Real> const& state_new = mf_cons_new.array(mfi);
        Array4<Real const> const& Hz = context.hz->const_array(mfi);
        Array4<Real const> const& z_w = context.z_w->const_array(mfi);
        Array4<Real const> const& srflx = context.shortwave->const_array(mfi);
        Array4<Real const> const& mskr = context.mask->const_array(mfi);
        Array4<Real const> const& bustr = context.bottom_stress_u->const_array(mfi);
        Array4<Real const> const& bvstr = context.bottom_stress_v->const_array(mfi);
        // ROMS takes the gas-exchange transfer velocity from Uwind/Vwind
        // under BULK_FLUXES and from the surface stress otherwise; only the
        // selected pair is read.
        Array4<Real const> sustr;
        Array4<Real const> svstr;
        Array4<Real const> uwind;
        Array4<Real const> vwind;
        if (do_bulk_flux) {
            uwind = context.wind_u->const_array(mfi);
            vwind = context.wind_v->const_array(mfi);
        } else {
            sustr = context.stress_u->const_array(mfi);
            svstr = context.stress_v->const_array(mfi);
        }

        ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            constexpr Real zero = Real(0.0);
            // Land columns contribute nothing: the increment below is scaled by rmask.
            if (mskr(i,j,0) == zero) return;

            constexpr Real one = Real(1.0);
            constexpr Real two = Real(2.0);
            constexpr Real eps = Real(1.0e-20);
            constexpr Real minval = Real(1.0e-6);
            constexpr Real cff_weno = Real(1.0e-14);

            constexpr Real rOxNO3 = Real(8.625);       // 138/16
            constexpr Real rOxNH4 = Real(6.625);       // 106/16
            constexpr Real rOxNH4Denit = Real(115.0) / Real(16.0);
            constexpr Real denitrification_NH4_fraction = Real(4.0) / Real(16.0);
            constexpr Real K_Nitri = Real(1.0);        // O2 half-saturation, nitrification
            constexpr Real K_Denit = K_Nitri;          // O2 half-saturation, denitrification
            constexpr Real K_wNO3 = Real(0.5);         // NO3 half-saturation, denitrification
            constexpr Real K_bO2 = Real(26.5);         // O2 half-saturation, sediment denitrification
            constexpr Real TIC_min = Real(1.0);
            constexpr Real TIC_max = Real(3000.0);

            // Compute inverse thickness to avoid repeated divisions.
            for (int k = 0; k < N; ++k) {
                inv_hz(i,j,k) = one / Hz(i,j,k);
            }
            for (int k = 0; k < N - 1; ++k) {
                inv_hz2(i,j,k) = one / (Hz(i,j,k) + Hz(i,j,k+1));
            }
            for (int k = 1; k < N - 1; ++k) {
                inv_hz3(i,j,k) = one / (Hz(i,j,k-1) + Hz(i,j,k) + Hz(i,j,k+1));
            }

            // Extract biological variables from tracer arrays, place them into
            // scratch arrays, and restrict their values to be positive definite.
            // TIC is also bounded, as in ROMS.
            for (int k = 0; k < N; ++k) {
                no3(i,j,k) = amrex::max(zero, state_old(i,j,k,c_no3));
                nh4(i,j,k) = amrex::max(zero, state_old(i,j,k,c_nh4));
                chlo(i,j,k) = amrex::max(zero, state_old(i,j,k,c_chlo));
                phyt(i,j,k) = amrex::max(zero, state_old(i,j,k,c_phyt));
                zoop(i,j,k) = amrex::max(zero, state_old(i,j,k,c_zoop));
                lden(i,j,k) = amrex::max(zero, state_old(i,j,k,c_lden));
                sden(i,j,k) = amrex::max(zero, state_old(i,j,k,c_sden));
                don(i,j,k) = amrex::max(zero, state_old(i,j,k,c_don));
                po4(i,j,k) = amrex::max(zero, state_old(i,j,k,c_po4));
                ldec(i,j,k) = amrex::max(zero, state_old(i,j,k,c_ldec));
                sdec(i,j,k) = amrex::max(zero, state_old(i,j,k,c_sdec));
                tic(i,j,k) = amrex::max(amrex::min(amrex::max(zero, state_old(i,j,k,c_tic)),
                                                   TIC_max), TIC_min);
                talk(i,j,k) = amrex::max(zero, state_old(i,j,k,c_talk));
                doc(i,j,k) = amrex::max(zero, state_old(i,j,k,c_doc));
                oxyg(i,j,k) = amrex::max(zero, state_old(i,j,k,c_oxyg));
                rdon_s(i,j,k) = amrex::max(zero, state_old(i,j,k,c_rdon));
                rdoc_s(i,j,k) = amrex::max(zero, state_old(i,j,k,c_rdoc));
                calc(i,j,k) = amrex::max(zero, state_old(i,j,k,c_calc));
                temp(i,j,k) = amrex::min(state_old(i,j,k,Temp_comp), Real(35.0));
                salt(i,j,k) = amrex::max(state_old(i,j,k,Salt_comp), zero);
            }

            // Calculate surface Photosynthetically Available Radiation (PAR).
            // REMORA stores srflx in W/m2, so only PARfrac applies. ROMS ecb.h
            // multiplies by rho0*Cp because its srflx is kinematic (degC m/s).
            const Real PARsur = parms.PARfrac * srflx(i,j,0);

            // Subcycle the complete ordered ECB reaction/column calculation.
            for (int iter = 0; iter < parms.BioIter; ++iter) {

                // Primary production, exudation and nitrification. Light is
                // attenuated from the surface down, so PAR at a depth depends
                // on the whole column above it.
                if (PARsur > zero) {
                    Real PAR = PARsur;
                    for (int k = N - 1; k >= 0; --k) {
                        const Real attfac = amrex::max(doc(i,j,k) + rdoc_s(i,j,k) - Real(70.819), zero) *
                                            Real(3.786e-3);
                        const Real Att = (parms.AttSW + parms.AttChl * chlo(i,j,k) + attfac) *
                                         (z_w(i,j,k+1) - z_w(i,j,k));
                        const Real ExpAtt = std::exp(-Att);
                        const Real Itop = PAR;
                        // Mean PAR over the cell; Itop when nothing attenuates.
                        PAR = (Att > zero) ? Itop * (one - ExpAtt) / Att : Itop;

                        const Real cff = parms.PhyCN * Real(12.0);
                        const Real Chl2C = amrex::min(chlo(i,j,k) / (phyt(i,j,k) * cff + eps), parms.Chl2C_m);
                        const Real Vp = amrex::max(parms.Vp0 * std::exp(Real(0.08065) * temp(i,j,k)), Real(4.00));
                        const Real fac1 = PAR * parms.PhyIS;
                        const Real Epp = Vp / std::sqrt(Vp * Vp + fac1 * fac1);
                        const Real t_PPmax = Epp * fac1;

                        // Nutrient limitation.
                        Real cff1 = nh4(i,j,k) * parms.K_NH4;
                        const Real cff2 = no3(i,j,k) * parms.K_NO3;
                        const Real inhNH4 = one / (one + cff1);
                        const Real L_NH4 = cff1 / (one + cff1);
                        const Real L_NO3 = cff2 * inhNH4 / (one + cff2);
                        const Real LTOT = L_NO3 + L_NH4;
                        const Real cff3 = po4(i,j,k) * parms.K_PO4;
                        const Real L_PO4 = cff3 / (one + cff3);
                        const Real LMIN = amrex::min(LTOT, L_PO4);

                        // Nutrient uptake.
                        const Real mu = dtdays * t_PPmax * LMIN;
                        const Real cff4 = mu * phyt(i,j,k) * L_NO3 /
                                          amrex::max(minval, LTOT) / amrex::max(minval, no3(i,j,k));
                        const Real cff5 = mu * phyt(i,j,k) * L_NH4 /
                                          amrex::max(minval, LTOT) / amrex::max(minval, nh4(i,j,k));
                        const Real cff6 = parms.R_P2N * mu * phyt(i,j,k) / amrex::max(minval, po4(i,j,k));
                        no3(i,j,k) = no3(i,j,k) / (one + cff4);
                        nh4(i,j,k) = nh4(i,j,k) / (one + cff5);
                        po4(i,j,k) = po4(i,j,k) / (one + cff6);
                        const Real N_Flux_NewProd = no3(i,j,k) * cff4;
                        const Real N_Flux_RegProd = nh4(i,j,k) * cff5;

                        // Phytoplankton growth and DON exudation.
                        Real f_ntr = oxyg(i,j,k) / (oxyg(i,j,k) + K_Nitri);
                        const Real f_dnf = K_Denit / (oxyg(i,j,k) + K_Denit) * no3(i,j,k) /
                                           (no3(i,j,k) + K_wNO3);
                        phyt(i,j,k) = phyt(i,j,k) + (one - parms.EsDON - (f_ntr + f_dnf) * parms.ElDON) *
                                                    (N_Flux_NewProd + N_Flux_RegProd);
                        chlo(i,j,k) = chlo(i,j,k) + (dtdays * t_PPmax * t_PPmax * LMIN * LMIN *
                                                     parms.Chl2C_m * chlo(i,j,k)) /
                                                    (parms.PhyIS * amrex::max(Chl2C, eps) * PAR + eps);
                        don(i,j,k) = don(i,j,k) + parms.EsDON * (N_Flux_NewProd + N_Flux_RegProd);
                        nh4(i,j,k) = nh4(i,j,k) + (f_ntr + f_dnf) * parms.ElDON *
                                                  (N_Flux_NewProd + N_Flux_RegProd);
                        talk(i,j,k) = talk(i,j,k) + (f_ntr + f_dnf) * parms.ElDON *
                                                    (N_Flux_NewProd + N_Flux_RegProd);
                        oxyg(i,j,k) = oxyg(i,j,k) + rOxNO3 * N_Flux_NewProd + rOxNH4 * N_Flux_RegProd -
                                      rOxNH4 * f_ntr * parms.ElDON * (N_Flux_NewProd + N_Flux_RegProd);
                        oxyg(i,j,k) = amrex::max(oxyg(i,j,k), zero);

                        // Carbon uptake and carbon-excess DOC exudation.
                        cff1 = parms.PhyCN * (N_Flux_NewProd + N_Flux_RegProd);
                        cff1 = amrex::min(cff1, tic(i,j,k) - eps);
                        tic(i,j,k) = tic(i,j,k) - cff1;
                        Real C_exc_up = parms.gammaC * parms.PhyCN * dtdays * t_PPmax * phyt(i,j,k) *
                                        (one - LTOT);
                        C_exc_up = amrex::min(C_exc_up, amrex::max(tic(i,j,k), zero) / parms.slCexc);
                        tic(i,j,k) = tic(i,j,k) - C_exc_up * parms.slCexc;
                        doc(i,j,k) = doc(i,j,k) + C_exc_up * parms.slCexc;
                        oxyg(i,j,k) = oxyg(i,j,k) + C_exc_up;
                        doc(i,j,k) = doc(i,j,k) + parms.PhyCN * parms.EsDON * (N_Flux_NewProd + N_Flux_RegProd);
                        tic(i,j,k) = tic(i,j,k) + parms.PhyCN * parms.ElDON * (f_ntr + f_dnf) *
                                                  (N_Flux_NewProd + N_Flux_RegProd);
                        talk(i,j,k) = talk(i,j,k) + N_Flux_NewProd - N_Flux_RegProd;

                        // Light-inhibited nitrification.
                        const Real fac2 = amrex::max(oxyg(i,j,k), zero);
                        const Real fac3 = amrex::max(fac2 / (K_Nitri + fac2), zero);
                        const Real fac1_n = dtdays * parms.NitriR * fac3;
                        const Real cff1_n = (PAR - parms.I_thNH4) / (parms.D_p5NH4 + PAR - two * parms.I_thNH4);
                        const Real cff2_n = one - amrex::max(zero, cff1_n);
                        const Real cff3_n = fac1_n * cff2_n;
                        nh4(i,j,k) = nh4(i,j,k) / (one + cff3_n);
                        const Real N_Flux_Nitrifi = nh4(i,j,k) * cff3_n;
                        no3(i,j,k) = no3(i,j,k) + N_Flux_Nitrifi;
                        oxyg(i,j,k) = oxyg(i,j,k) - two * N_Flux_Nitrifi;
                        talk(i,j,k) = talk(i,j,k) - two * N_Flux_Nitrifi;

                        // Light attenuation at the bottom of the grid cell, the
                        // starting PAR value for the next (deeper) cell.
                        PAR = Itop * ExpAtt;
                    }
                } else {
                    // Night: nitrification without light inhibition.
                    for (int k = N - 1; k >= 0; --k) {
                        const Real f_ntr = oxyg(i,j,k) / (oxyg(i,j,k) + K_Nitri);
                        const Real cff3 = dtdays * parms.NitriR * f_ntr;
                        nh4(i,j,k) = nh4(i,j,k) / (one + cff3);
                        const Real N_Flux_Nitrifi = nh4(i,j,k) * cff3;
                        no3(i,j,k) = no3(i,j,k) + N_Flux_Nitrifi;
                        oxyg(i,j,k) = oxyg(i,j,k) - two * N_Flux_Nitrifi;
                        oxyg(i,j,k) = amrex::max(zero, oxyg(i,j,k));
                        talk(i,j,k) = talk(i,j,k) - two * N_Flux_Nitrifi;
                    }
                }

                // Phytoplankton grazing by zooplankton, sloppy feeding and
                // phytoplankton mortality.
                {
                    const Real cff2 = dtdays * parms.PhyMR;
                    for (int k = 0; k < N; ++k) {
                        const Real fac1 = dtdays * parms.ZooGR * std::exp(Real(0.0875) * temp(i,j,k));
                        const Real cff1 = fac1 * zoop(i,j,k) * phyt(i,j,k) /
                                          (parms.K_Phy + phyt(i,j,k) * phyt(i,j,k));
                        const Real cff3 = one / (one + cff1);
                        phyt(i,j,k) = cff3 * phyt(i,j,k);
                        chlo(i,j,k) = cff3 * chlo(i,j,k);
                        const Real N_Flux_Assim = cff1 * phyt(i,j,k) * parms.ZooAE_N;
                        const Real N_Flux_Egest = phyt(i,j,k) * cff1 * (one - parms.ZooAE_N);
                        zoop(i,j,k) = zoop(i,j,k) + N_Flux_Assim;
                        const Real G_Gmax = phyt(i,j,k) * phyt(i,j,k) /
                                            (parms.K_Phy + phyt(i,j,k) * phyt(i,j,k));
                        const Real qDON = Real(0.71) * G_Gmax;
                        lden(i,j,k) = lden(i,j,k) + N_Flux_Egest * (one - qDON);
                        const Real N_Flux_Sloppy_ldon = N_Flux_Egest * qDON * (one - parms.deltN);
                        const Real N_Flux_Sloppy_sldon = N_Flux_Egest * qDON * parms.deltN;
                        don(i,j,k) = don(i,j,k) + N_Flux_Sloppy_sldon;
                        nh4(i,j,k) = nh4(i,j,k) + N_Flux_Sloppy_ldon;
                        talk(i,j,k) = talk(i,j,k) + N_Flux_Sloppy_ldon;
                        const Real N_Flux_Pmortal = cff2 * amrex::max(phyt(i,j,k) - parms.PhyMin, zero);
                        phyt(i,j,k) = phyt(i,j,k) - N_Flux_Pmortal;
                        chlo(i,j,k) = chlo(i,j,k) - cff2 * amrex::max(chlo(i,j,k) - parms.ChlMin, zero);
                        sden(i,j,k) = sden(i,j,k) + N_Flux_Pmortal;
                        const Real qDOC = qDON;
                        ldec(i,j,k) = ldec(i,j,k) + parms.PhyCN * N_Flux_Egest * (one - qDOC);
                        tic(i,j,k) = tic(i,j,k) + parms.PhyCN * N_Flux_Egest * qDOC * (one - parms.deltC);
                        doc(i,j,k) = doc(i,j,k) + parms.PhyCN * N_Flux_Egest * qDOC * parms.deltC;
                        sdec(i,j,k) = sdec(i,j,k) + parms.PhyCN * N_Flux_Pmortal;
                        tic(i,j,k) = tic(i,j,k) + (parms.PhyCN - parms.ZooCN) * N_Flux_Assim;
                        oxyg(i,j,k) = oxyg(i,j,k) - rOxNH4 * N_Flux_Sloppy_ldon;
                        oxyg(i,j,k) = amrex::max(zero, oxyg(i,j,k));
                    }
                }

                // Zooplankton mortality, excretion and basal metabolism.
                {
                    const Real cff1 = dtdays * parms.ZooBM;
                    const Real fac2 = dtdays * parms.ZooMR;
                    const Real fac3 = dtdays * parms.ZooER;
                    for (int k = 0; k < N; ++k) {
                        const Real fac1 = fac3 * phyt(i,j,k) * phyt(i,j,k) /
                                          (parms.K_Phy + phyt(i,j,k) * phyt(i,j,k));
                        const Real cff2 = fac2 * zoop(i,j,k);
                        const Real cff3 = fac1 * parms.ZooAE_N;
                        zoop(i,j,k) = zoop(i,j,k) / (one + cff2 + cff3);
                        const Real N_Flux_Zmortal = cff2 * zoop(i,j,k);
                        const Real N_Flux_Zexcret = cff3 * zoop(i,j,k);
                        nh4(i,j,k) = nh4(i,j,k) + N_Flux_Zexcret;
                        po4(i,j,k) = po4(i,j,k) + parms.R_P2N * N_Flux_Zexcret;
                        sden(i,j,k) = sden(i,j,k) + N_Flux_Zmortal;
                        const Real N_Flux_Zmetabo = cff1 * amrex::max(zoop(i,j,k) - parms.ZooMin, zero);
                        zoop(i,j,k) = zoop(i,j,k) - N_Flux_Zmetabo;
                        nh4(i,j,k) = nh4(i,j,k) + N_Flux_Zmetabo;
                        po4(i,j,k) = po4(i,j,k) + parms.R_P2N * N_Flux_Zmetabo;
                        oxyg(i,j,k) = oxyg(i,j,k) - rOxNH4 * (N_Flux_Zmetabo + N_Flux_Zexcret);
                        oxyg(i,j,k) = amrex::max(zero, oxyg(i,j,k));
                        sdec(i,j,k) = sdec(i,j,k) + parms.ZooCN * N_Flux_Zmortal;
                        tic(i,j,k) = tic(i,j,k) + parms.ZooCN * (N_Flux_Zmetabo + N_Flux_Zexcret);
                        talk(i,j,k) = talk(i,j,k) + N_Flux_Zmetabo + N_Flux_Zexcret;
                    }
                }

                // Coagulation of phytoplankton and small detritus into large detritus.
                {
                    const Real fac1 = dtdays * parms.CoagR;
                    for (int k = 0; k < N; ++k) {
                        const Real cff1 = fac1 * (sden(i,j,k) + phyt(i,j,k));
                        const Real cff2 = one / (one + cff1);
                        phyt(i,j,k) = phyt(i,j,k) * cff2;
                        chlo(i,j,k) = chlo(i,j,k) * cff2;
                        sden(i,j,k) = sden(i,j,k) * cff2;
                        sdec(i,j,k) = sdec(i,j,k) * cff2;
                        const Real N_Flux_CoagP = phyt(i,j,k) * cff1;
                        const Real N_Flux_CoagD = sden(i,j,k) * cff1;
                        lden(i,j,k) = lden(i,j,k) + N_Flux_CoagP + N_Flux_CoagD;
                        ldec(i,j,k) = ldec(i,j,k) + parms.PhyCN * N_Flux_CoagP + sdec(i,j,k) * cff1;
                    }
                }

                // Nitrogen remineralization of detritus and semilabile DON,
                // aerobic and by denitrification.
                for (int k = 0; k < N; ++k) {
                    const Real f_ntr = oxyg(i,j,k) / (oxyg(i,j,k) + K_Nitri);
                    const Real f_dnf = K_Denit / (oxyg(i,j,k) + K_Denit) * no3(i,j,k) / (no3(i,j,k) + K_wNO3);
                    const Real cff = dtdays * std::exp(Real(0.0875) * temp(i,j,k));
                    const Real cff1 = cff * parms.SDeRRN;
                    const Real cff2 = cff1 * parms.deltN;
                    const Real cff3 = cff1 * (one - parms.deltN);
                    const Real cff4 = cff * parms.LDeRRN;
                    const Real cff5 = cff4 * parms.deltN;
                    const Real cff6 = cff4 * (one - parms.deltN);
                    const Real rdon = cff * parms.RDeRRN;
                    const Real fac = (f_ntr + f_dnf);
                    sden(i,j,k) = sden(i,j,k) / (one + cff2 + cff3 * fac);
                    lden(i,j,k) = lden(i,j,k) / (one + cff5 + cff6 * fac);
                    don(i,j,k) = don(i,j,k) / (one + rdon * fac);
                    const Real N_Flux_RemineS = sden(i,j,k) * cff3 * fac;
                    const Real N_Flux_RemineL = lden(i,j,k) * cff6 * fac;
                    const Real N_Flux_RemineR = don(i,j,k) * rdon * fac;
                    nh4(i,j,k) = nh4(i,j,k) + N_Flux_RemineS + N_Flux_RemineL + N_Flux_RemineR;
                    po4(i,j,k) = po4(i,j,k) + parms.R_P2N * (N_Flux_RemineS + N_Flux_RemineL + N_Flux_RemineR);
                    const Real cff7 = (don(i,j,k) * rdon + sden(i,j,k) * cff3 + lden(i,j,k) * cff6) * rOxNH4 /
                                      (oxyg(i,j,k) + K_Nitri);
                    oxyg(i,j,k) = oxyg(i,j,k) / (one + cff7);
                    talk(i,j,k) = talk(i,j,k) + N_Flux_RemineS + N_Flux_RemineL + N_Flux_RemineR;
                    const Real cff8 = (sden(i,j,k) * cff3 + lden(i,j,k) * cff6 + don(i,j,k) * rdon) * Real(5.3) /
                                      (oxyg(i,j,k) + K_Denit) * K_Denit / (no3(i,j,k) + K_wNO3);
                    no3(i,j,k) = no3(i,j,k) / (one + cff8);
                    talk(i,j,k) = talk(i,j,k) + no3(i,j,k) * cff8;
                    don(i,j,k) = don(i,j,k) + sden(i,j,k) * cff2 + lden(i,j,k) * cff5;
                }

                // Surface O2 gas exchange.
                const Real u10squ = u10_squared(do_bulk_flux, uwind, vwind, sustr, svstr, i, j, rho0);
                {
                    const int k = N - 1;
                    const Real cff2 = transfer_rate(o2_schmidt, dtdays);
                    const Real SchmidtN_Ox = schmidt_number(o2_schmidt, temp(i,j,k));
                    const Real cff3 = transfer_velocity(cff2, u10squ, SchmidtN_Ox);
                    const Real O2satu = o2_saturation(temp(i,j,k), salt(i,j,k));
                    const Real O2_Flux = cff3 * (O2satu - oxyg(i,j,k));
                    oxyg(i,j,k) = oxyg(i,j,k) + O2_Flux * inv_hz(i,j,k);
                }

                // Carbon remineralization of detritus and semilabile DOC.
                for (int k = 0; k < N; ++k) {
                    const Real f_ntr = oxyg(i,j,k) / (oxyg(i,j,k) + K_Nitri);
                    const Real f_dnf = K_Denit / (oxyg(i,j,k) + K_Denit) * no3(i,j,k) / (no3(i,j,k) + K_wNO3);
                    const Real cff = dtdays * std::exp(Real(0.0875) * temp(i,j,k));
                    const Real cff1 = cff * parms.SDeRRC;
                    const Real cff2 = cff1 * parms.deltC;
                    const Real cff3 = cff1 * (one - parms.deltC);
                    const Real cff4 = cff * parms.LDeRRC;
                    const Real cff5 = cff4 * parms.deltC;
                    const Real cff6 = cff4 * (one - parms.deltC);
                    const Real rdoc = cff * parms.RDeRRC;
                    const Real fac = (f_ntr + f_dnf);
                    sdec(i,j,k) = sdec(i,j,k) / (one + cff2 + cff3 * fac);
                    ldec(i,j,k) = ldec(i,j,k) / (one + cff5 + cff6 * fac);
                    doc(i,j,k) = doc(i,j,k) / (one + rdoc * fac);
                    const Real C_Flux_RemineS = sdec(i,j,k) * cff3 * fac;
                    const Real C_Flux_RemineL = ldec(i,j,k) * cff6 * fac;
                    const Real C_Flux_RemineR = doc(i,j,k) * rdoc * fac;
                    tic(i,j,k) = tic(i,j,k) + C_Flux_RemineS + C_Flux_RemineL + C_Flux_RemineR;
                    doc(i,j,k) = doc(i,j,k) + cff2 * sdec(i,j,k) + cff5 * ldec(i,j,k);
                }

                // Surface CO2 gas exchange.
                {
                    const int k = N - 1;
                    const Real pCO2 = pco2_water(temp(i,j,k), salt(i,j,k), tic(i,j,k), talk(i,j,k));
                    const Real cff2 = transfer_rate(co2_schmidt, dtdays);
                    const Real SchmidtN = schmidt_number(co2_schmidt, temp(i,j,k));
                    const Real cff3 = transfer_velocity(cff2, u10squ, SchmidtN);
                    const Real CO2_sol = co2_solubility(temp(i,j,k), salt(i,j,k));
                    const Real CO2_Flux = cff3 * CO2_sol * (parms.pCO2air - pCO2);
                    tic(i,j,k) = tic(i,j,k) + CO2_Flux * inv_hz(i,j,k);
                }

                // Vertical sinking terms: phytoplankton, chlorophyll, small and
                // large nitrogen detritus, and small and large carbon detritus,
                // by the ROMS PPM/WENO scheme with its monotone boundary closures.
                for (int isink = 0; isink < 6; ++isink) {
                    Array4<Real> bio = phyt;
                    Real wbio = parms.wPhy;
                    if (isink == 1) {
                        bio = chlo;
                        wbio = parms.wPhy;
                    } else if (isink == 2) {
                        bio = sden;
                        wbio = parms.wSDet;
                    } else if (isink == 3) {
                        bio = lden;
                        wbio = parms.wLDet;
                    } else if (isink == 4) {
                        bio = sdec;
                        wbio = parms.wSDet;
                    } else if (isink == 5) {
                        bio = ldec;
                        wbio = parms.wLDet;
                    }
                    const bool phyt_sink = (isink == 0);
                    const bool nitrogen_sink = (isink == 0 || isink == 2 || isink == 3);
                    const bool carbon_sink = (isink == 4 || isink == 5);

                    for (int k = 0; k < N; ++k) {
                        qc(i,j,k) = amrex::max(bio(i,j,k), zero);
                    }

                    // Reconstruct vertical profile of selected biological
                    // constituent "Bio" in terms of a set of parabolic segments
                    // within each grid box. Then, compute semi-Lagrangian flux
                    // due to sinking.
                    for (int iface = N - 1; iface >= 1; --iface) {
                        FC(i,j,iface) = (qc(i,j,iface) - qc(i,j,iface-1)) * inv_hz2(i,j,iface-1);
                    }
                    for (int k = 1; k < N - 1; ++k) {
                        Real dltR = Hz(i,j,k) * FC(i,j,k+1);
                        Real dltL = Hz(i,j,k) * FC(i,j,k);
                        Real cff = Hz(i,j,k-1) + two * Hz(i,j,k) + Hz(i,j,k+1);
                        const Real cffR = cff * FC(i,j,k+1);
                        const Real cffL = cff * FC(i,j,k);

                        // Apply PPM monotonicity constraint to prevent
                        // oscillations within the grid box.
                        if ((dltR * dltL) <= zero) {
                            dltR = zero;
                            dltL = zero;
                        } else if (std::abs(dltR) > std::abs(cffL)) {
                            dltR = cffL;
                        } else if (std::abs(dltL) > std::abs(cffR)) {
                            dltL = cffR;
                        }

                        // Compute right and left side values (bR,bL) of
                        // parabolic segments within grid box Hz(k); (WR,WL) are
                        // measures of quadratic variations.
                        cff = (dltR - dltL) * inv_hz3(i,j,k);
                        dltR = dltR - cff * Hz(i,j,k+1);
                        dltL = dltL + cff * Hz(i,j,k-1);
                        bR(i,j,k) = qc(i,j,k) + dltR;
                        bL(i,j,k) = qc(i,j,k) - dltL;
                        WR(i,j,k) = std::pow((two * dltR - dltL), 2);
                        WL(i,j,k) = std::pow((dltR - two * dltL), 2);
                    }
                    for (int k = 1; k < N - 2; ++k) {
                        const Real dltL = amrex::max(cff_weno, WL(i,j,k));
                        const Real dltR = amrex::max(cff_weno, WR(i,j,k+1));
                        bR(i,j,k) = (dltR * bR(i,j,k) + dltL * bL(i,j,k+1)) / (dltR + dltL);
                        bL(i,j,k+1) = bR(i,j,k);
                    }
                    FC(i,j,N) = zero;              // no-flux boundary condition
                    bR(i,j,N-1) = qc(i,j,N-1);     // default strictly monotonic
                    bL(i,j,N-1) = qc(i,j,N-1);     // conditions
                    bR(i,j,N-2) = qc(i,j,N-1);
                    bL(i,j,1) = qc(i,j,0);         // bottom grid boxes are
                    bR(i,j,0) = qc(i,j,0);         // re-assumed to be
                    bL(i,j,0) = qc(i,j,0);         // piecewise constant.

                    // Apply monotonicity constraint again, since the reconciled
                    // interfacial values may cause a non-monotonic behavior of
                    // the parabolic segments inside the grid box.
                    for (int k = 0; k < N; ++k) {
                        Real dltR = bR(i,j,k) - qc(i,j,k);
                        Real dltL = qc(i,j,k) - bL(i,j,k);
                        const Real cffR = two * dltR;
                        const Real cffL = two * dltL;
                        if ((dltR * dltL) < zero) {
                            dltR = zero;
                            dltL = zero;
                        } else if (std::abs(dltR) > std::abs(cffL)) {
                            dltR = cffL;
                        } else if (std::abs(dltL) > std::abs(cffR)) {
                            dltL = cffR;
                        }
                        bR(i,j,k) = qc(i,j,k) + dltR;
                        bL(i,j,k) = qc(i,j,k) - dltL;
                    }

                    // After this moment reconstruction is considered complete.
                    // The next stage is to compute vertical advective fluxes,
                    // FC. It is expected that sinking may occur relatively fast,
                    // the algorithm is designed to be free of CFL criterion,
                    // which is achieved by allowing integration bounds for
                    // semi-Lagrangian advective flux to use as many grid boxes
                    // in upstream direction as necessary.
                    const Real cff_w = dtdays * std::abs(wbio);
                    for (int k = 0; k < N; ++k) {
                        FC(i,j,k) = zero;
                        WL(i,j,k) = z_w(i,j,k) + cff_w;
                        WR(i,j,k) = Hz(i,j,k) * qc(i,j,k);
                    }
                    for (int k = 0; k < N; ++k) {
                        int ksource = k;
                        for (int ks = k; ks < N - 1; ++ks) {
                            if (WL(i,j,k) > z_w(i,j,ks+1)) {
                                ksource = ks + 1;
                                FC(i,j,k) = FC(i,j,k) + WR(i,j,ks);
                            }
                        }
                        const int ks = ksource;
                        const Real cu = amrex::min(one, (WL(i,j,k) - z_w(i,j,ks)) * inv_hz(i,j,ks));
                        FC(i,j,k) = FC(i,j,k) +
                                    Hz(i,j,ks) * cu *
                                    (bL(i,j,ks) + cu * (Real(0.5) * (bR(i,j,ks) - bL(i,j,ks)) -
                                     (Real(1.5) - cu) * (bR(i,j,ks) + bL(i,j,ks) - two * qc(i,j,ks))));
                    }
                    for (int k = 0; k < N; ++k) {
                        bio(i,j,k) = qc(i,j,k) + (FC(i,j,k+1) - FC(i,j,k)) * inv_hz(i,j,k);
                    }

                    // Stress-dependent resuspension and instantaneous sediment
                    // recycling and burial of the material reaching the bottom.
                    Real cff1 = zero;
                    Real cbe = zero;
                    if (nitrogen_sink || carbon_sink) {
                        const Real ustarb = std::sqrt(std::sqrt(
                            std::pow((Real(0.5) * (bustr(i,j,0) + bustr(i+1,j,0))), 2) +
                            std::pow((Real(0.5) * (bvstr(i,j,0) + bvstr(i,j+1,0))), 2)));
                        const Real resuspr = amrex::min(std::pow(ustarb, 2) / Real(9.61e-6), one);
                        cff1 = FC(i,j,0) * inv_hz(i,j,0) * (one - resuspr);
                        const Real cff5 = FC(i,j,0) * inv_hz(i,j,0) * resuspr;
                        if (nitrogen_sink) {
                            sden(i,j,0) = sden(i,j,0) + cff5;
                        }
                        if (carbon_sink) {
                            sdec(i,j,0) = sdec(i,j,0) + cff5;
                        }
                        if (phyt_sink) {
                            sdec(i,j,0) = sdec(i,j,0) + cff5 * parms.PhyCN;
                        }
                    }
                    if (nitrogen_sink) {
                        // Nitrogen burial efficiency, then sediment
                        // denitrification of the rest.
                        const Real CNbur = phyt_sink ? parms.PhyCN : Real(9.3);
                        const Real cff5 = Real(0.01);
                        const Real fac1 = CNbur * Real(12.0) / Real(1000.0) * cff1 * Hz(i,j,0) *
                                          (Real(365.0) / dtdays);
                        cbe = amrex::min(Real(0.50), Real(0.023) * std::pow(fac1, Real(0.5797)));
                        cff1 = (one - cbe) * cff1;
                        const Real O2satb = o2_saturation(temp(i,j,0), salt(i,j,0));
                        const Real LBO2 = (K_bO2 * (O2satb - oxyg(i,j,0))) / (O2satb * (oxyg(i,j,0) + K_bO2));
                        const Real fac = cff1 * (one + Real(3.0) * LBO2);
                        nh4(i,j,0) = nh4(i,j,0) + denitrification_NH4_fraction * fac;
                        talk(i,j,0) = talk(i,j,0) + denitrification_NH4_fraction * fac;
                        don(i,j,0) = don(i,j,0) + cff5 * fac;
                        oxyg(i,j,0) = oxyg(i,j,0) - cff1 * rOxNH4Denit * (one - LBO2);
                    }
                    if (phyt_sink || carbon_sink) {
                        // Carbon burial efficiency; phytoplankton reuses the
                        // flux before the nitrogen burial above.
                        if (phyt_sink) {
                            cff1 = cff1 / (one - cbe);
                        }
                        Real fac1 = Real(12.0) / Real(1000.0) * cff1 * Hz(i,j,0) * (Real(365.0) / dtdays);
                        if (phyt_sink) {
                            fac1 = fac1 * parms.PhyCN;
                        }
                        cbe = amrex::min(Real(0.50), Real(0.023) * std::pow(fac1, Real(0.5797)));
                        const Real fac2 = phyt_sink ? parms.PhyCN : one;
                        cff1 = (one - cbe) * cff1;
                        const Real cff5 = Real(0.01);
                        const Real fac = fac2 * cff1;
                        doc(i,j,0) = doc(i,j,0) + fac * cff5;
                        tic(i,j,0) = tic(i,j,0) + fac * (one - cff5);
                    }
                }
            }

            for (int k = 0; k < N; ++k) {
                tic(i,j,k) = amrex::min(tic(i,j,k), TIC_max);
                tic(i,j,k) = amrex::max(tic(i,j,k), TIC_min);
            }

            // Update global tracer variables: add the Hz-weighted increment,
            // relative to the same restricted old values the scratch started
            // from, to the incoming state.
            const Real rmask = mskr(i,j,0);
            for (int k = 0; k < N; ++k) {
                auto add_increment = [&] (Array4<Real> const& bio, int c) {
                    const Real old_value = amrex::max(zero, state_old(i,j,k,c));
                    state_new(i,j,k,c) += (bio(i,j,k) - old_value) * rmask * Hz(i,j,k);
                };
                add_increment(no3, c_no3);
                add_increment(nh4, c_nh4);
                add_increment(chlo, c_chlo);
                add_increment(phyt, c_phyt);
                add_increment(zoop, c_zoop);
                add_increment(lden, c_lden);
                add_increment(sden, c_sden);
                add_increment(don, c_don);
                add_increment(po4, c_po4);
                add_increment(ldec, c_ldec);
                add_increment(sdec, c_sdec);
                const Real old_tic = amrex::max(amrex::min(amrex::max(zero, state_old(i,j,k,c_tic)),
                                                           TIC_max), TIC_min);
                state_new(i,j,k,c_tic) += (tic(i,j,k) - old_tic) * rmask * Hz(i,j,k);
                add_increment(talk, c_talk);
                add_increment(doc, c_doc);
                add_increment(oxyg, c_oxyg);
                add_increment(rdon_s, c_rdon);
                add_increment(rdoc_s, c_rdoc);
                add_increment(calc, c_calc);
            }
        });
    }
}
