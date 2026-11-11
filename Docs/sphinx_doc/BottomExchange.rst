.. _sec:bottom-exchange:

Water and Tracer Exchange Through the Bottom
============================================

The default bottom is impermeable. No additional input or exchange arrays are
needed for this case. ``remora.bottom_exchange_type`` selects ``none`` (default),
``constant``, ``analytic``, ``netcdf``, or ``coupler``. Bottom exchange is independent
of bottom drag and the existing diffusive bottom heat and tracer fluxes.

Velocity and tracer conventions
-------------------------------

All input fields have the horizontal locations of rho points and one value per
water column. The velocity ``bottom_q`` is volume flow per unit **horizontal**
area, in m/s, positive into the ocean. On a stationary bed :math:`z=-h(x,y)`,

.. math::

   q_b = w_b + u_b\partial_x h + v_b\partial_y h.

Thus zero ``bottom_q`` means no flow through a sloping bed. A Cartesian vertical
velocity or a unit-normal velocity must first be converted to this convention.
For an inward unit-normal speed :math:`v_n`,
:math:`q_b=v_n\sqrt{1+|\nabla h|^2}`.

REMORA converts this to volume transport :math:`Q_b=q_b/(pm\,pn)` and uses it in
both continuity and vertical advection, including the free-surface update. For
each tracer the bottom advective transport is

.. math::

   F_{b,c}=\begin{cases}
      Q_b C_{\mathrm{in},c}, & Q_b>0,\\
      Q_b C_{\mathrm{ocean},c}, & Q_b<0,\\
      0, & Q_b=0.
   \end{cases}

Incoming concentrations are required for every transported tracer, including
temperature, salinity, passive tracers, and enabled biology tracers. Their units
and names match the ocean state. They are properties of incoming water and do
not clamp the bottom ocean cell. An outflow uses the ocean-side donor and needs
no external concentrations. Land cells contribute zero exchange.

Incoming water carries zero horizontal grid-aligned u and v momentum. Outflow
removes the ocean-side horizontal momentum. Conservative volume and momentum
updates therefore include dilution of the ocean velocity by incoming water.

Constant and analytic input
---------------------------

A freshwater example with one passive tracer is::

   remora.bottom_exchange_type = constant
   remora.bottom_q = 1.e-6
   remora.bottom_temp = 12
   remora.bottom_salt = 0
   remora.bottom_tracer = 0

Use ``bottom_<name>`` for each name in REMORA's tracer registry. Missing
concentrations are rejected wherever inflow occurs; zero is an explicit choice.

For analytic forcing, override ``ProblemBase::bottom_exchange``. Its MultiFab has
component 0 for ``q`` and components 1 through ``ncons`` for concentrations in
registry order. The hook fills the level-0 rho slab and is evaluated at each
level's slow-step midpoint, using elapsed model seconds. No hook is required for
the default impermeable boundary.

NetCDF input
------------

Build with PNetCDF support and set::

   remora.bottom_exchange_type = netcdf
   remora.bottom_exchange_file = bottom.nc
   remora.bottom_time_name = bottom_time

The files contain ``bottom_q`` and ``bottom_<tracer-name>`` on dimensions
``(bottom_time, eta_rho, xi_rho)``. There is no vertical dimension. The model rho
dimensions include the ROMS boundary ring, so they are
``(n_cell_y + 2, n_cell_x + 2)``. Supply data on the level-0 model grid and keep
variables and dimensions consistent across the file set.

All fields are read through ``NCTimeSeries``. Supply at least two time records in
total; records may be split across files. Equal values at two different times
specify steady forcing. Static arrays and isolated single-record files are not
supported. The constant provider is also available for steady forcing.

The existing calendar, reference-date, cyclic-time and interpolation rules apply.
Time units are days, with a reference date consistent with ``remora.time_ref``.
A tracer absent from the files may be supplied explicitly as
``remora.bottom_<name>``. Outflow uses ocean concentrations. ``bottom_q`` must be
in the files.

Each level samples the interpolated fields at its slow-step midpoint and holds
the resulting transports through that slow step, including all fast barotropic
steps. Resolve changes of sign and rapidly varying forcing with the slow
timestep. The product of interpolated velocity and concentration defines the
incoming tracer transport at the sample time.

Coupler interface
-----------------

Set ``remora.bottom_exchange_type = coupler``. After ``InitData``:

.. code-block:: cpp

   amrex::BoxArray ba;
   amrex::DistributionMapping dm;
   ocean.GetBottomExchangeLayout(ba, dm);
   const auto& names = ocean.GetTracerNames();
   amrex::MultiFab q(ba, dm, 1, 0);
   amrex::MultiFab c(ba, dm, names.size(), 0);
   // Fill q [m/s] and c [ocean concentration units].
   ocean.SetBottomExchangeData(q, c, begin_time, end_time);

REMORA copies the valid data and owns its one-cell horizontal halos. The producer may reuse its
buffers immediately. A snapshot is held over the stated interval in elapsed
model seconds; it must cover every slow step using it. Missing or expired data
cause an error. Supply a new snapshot after restarting a coupled model.

The coupler handles remapping to REMORA's level-0 rho layout. Coordinate and
mask queries are available alongside the existing coupling API. The same
snapshot supplies all refinement levels. ``PackBottomState(lev, c)`` exports the
current bottom-cell concentrations on that level's rho slab, with land zeroed.
``GetBottomExchangeTransport(lev)`` exposes read-only water transport followed by
incoming tracer transports for the current step; it returns null for an
impermeable boundary or an unprepared layout.

Snapshots describe constant rates and incoming concentrations. If a coupler has
interval totals, it must preserve both water and tracer transfers when preparing
these fields; independently averaging q and concentration generally does not
preserve their product. Subdivide intervals containing inflow/outflow reversals.

Refinement and restart
----------------------

The provider's level-0 water transport and incoming tracer transports are split
equally among descendant cells in index space. Their sums are preserved across
refinement, including with curvilinear metrics. The corresponding local m/s
values follow from each child's physical area. This is a piecewise constant
transport specification; fine spatial structure requires a finer forcing grid
or a future remapping adapter.

Wet masks must be consistent across refinement wherever exchange is nonzero;
REMORA rejects a refined land cell that would discard part of a source. Regridding
invalidates the affected exchange slabs, and each substep samples its own forcing
time. Constant transports are reused until the layout changes;
coupled transports are reused until the snapshot or layout changes.
Existing coarse/fine dynamics and conservation options still apply; see
:ref:`Mesh Refinement <MeshRefinement>`.

File, constant, and analytic forcing are reconstructed from the input settings
after restart. Coupled forcing is supplied again by the driver. Exchange
buffers are derived data and do not change the checkpoint format.

GLS and timestep diagnostics
----------------------------

The interior GLS closure, its stability functions, and its buoyancy-production
term are unchanged. Freshwater in the bottom cell creates the same buoyancy
tendency as an equivalent lateral input to that cell. The existing treatment
of unstable bottom boundary layers remains applicable.

The weak-seepage approximation retains the existing stress/roughness wall
conditions and the near-bed turbulence reconstructions. Turbulence entering the
resolved wall layer is assumed to be equilibrated with those wall conditions;
there is no independently prescribed incoming TKE or GLS length scale. The
existing particle boundary behavior remains reflecting; water exchange does
not inject or remove Lagrangian particles.

With ``remora.bottom_exchange_diagnostics = true``, each slow step reports its
maximum source Courant number and, for GLS, the mechanical measures

.. math::

   \Pi_{\mathrm{transport}} & = \max_{j=m,k,\psi}\frac{|q_b|d}{K_j},\\
   \Pi_{\mathrm{adjustment}} & = \frac{|q_b|k}{d\epsilon},\\
   \Pi_{\mathrm{wall}} & = \frac{|q_b|U_1}{u_*^2},
   \qquad \epsilon=(c_\mu^0)^3 k^{3/2}/\ell.

The diagnostic uses half the bottom-cell thickness for d, the first interior
GLS values for k and length, and the current mixing coefficients. The wall
measure is undefined at zero physical stress. TKE floors are numerical seeds
and do not establish physical turbulence strength. Values much smaller than
one support the weak mechanical seepage approximation; 0.1 is a useful
screening tolerance, not a calibrated closure limit.

For a neutral log layer these measures scale with ``abs(q_b)/u_star``; the
wall measure additionally contains ``log(d/z0)/kappa``. Small seepage velocity
does not imply weak buoyancy forcing: the buoyancy flux is
:math:`B_b=q_b g(\rho_{\mathrm{ocean}}-\rho_{\mathrm{in}})/\rho_0`.

The source Courant number ``abs(q_b)*dt/Hz_bottom`` must not exceed
``remora.bottom_exchange_max_courant`` (default 0.5, allowed range (0,1]). The
adaptive timestep estimate includes the current imposed exchange. The actual
midpoint forcing is checked again before advancing, including for a fixed
timestep. Reduce the slow timestep if that check fails. Resolve convection and
forcing transients according to the accuracy requirements of the application.
