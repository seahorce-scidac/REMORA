
 .. role:: cpp(code)
    :language: c++

.. _sec:VerticalMixing:

Vertical Mixing Parametrizations
================================

These equations are as given in the `ROMS documentation`_

.. _`ROMS documentation`: https://www.myroms.org/wiki/Vertical_Mixing_Parameterizations

The type of vertical mixing parametrization can be chosen in the :ref:`inputs file`<list-of-parameters-15>` with the ``remora.vertical_mixing_type`` option.

Custom Mixing
-------------

By default, the vertical mixing coefficients :math:`K_v` and :math:`K_t` (``vec_Akv`` and ``vec_Akt`` in the code)
are defined as a custom function in ``init_custom_vmix`` in the problem's ``prob.cpp`` file. This function can be defined to be a function of simulation variables, most commonly depth :math:`z_w`. It is re-evaluated at every time step.

Generic Length Scale (GLS)
--------------------------

The Generic Length Scale (GLS) model is a two-equation turbulence closure scheme which can be tuned to behave like many of the traditional schemes such as Mellor and Yamada 2.5. This class of schemes add additional prognostic equations for the turbulent kinetic energy :math:`k=\frac{q^2}{2}` and the product :math:`kl`, where :math:`l` is a length scale. The evolution equation for TKE is:

.. math::
    \frac{D}{Dt}\left(\frac{q^2}{2}\right)-\frac{\partial}{\partial z} \left[K_q \frac{\partial}{\partial z} \left(\frac{q^2}{2}\right)\right] = P_s + P_b + \epsilon

where :math:`P_s` is the shear production, :math:`P_b` is the buoyant production, and :math:`\epsilon` is the dissipation. In model coordinates, this evolution equation becomes:

.. math::
    \frac{\partial}{\partial t} \left(\frac{H_z q^2}{mn}\right) &+ \frac{\partial}{\partial \xi} \left(\frac{H_z u q^2}{n}\right) + \frac{\partial}{\partial \eta}\left(\frac{H_z v q^2}{m}\right) + \frac{\partial}{\partial s}\left(\frac{H_z \Omega q^2}{mn}\right) - \frac{\partial}{\partial s}\left(\frac{K_q}{mnH_z}\frac{\partial q^2}{\partial s}\right)\\
    &= \frac{2 H_z}{mn}\left(P_s + P_b + \epsilon\right).

The terms on the right-hand side are:

.. math::
    P_s &= K_v \left[\left(\frac{\partial u}{\partial z}\right)^2 + \left(\frac{\partial v}{\partial z}\right)^2\right]\\
    P_b &= -K_t N^2\\
    \epsilon &= \left(c_mu^0\right)^{3+p/n} k^{3/2+m/n} \psi^{-1/n}.


The parameter :math:`\psi` is used to extablished the turbulence length scale. Its evolution equation is:

.. math::
    \frac{D\psi}{Dt} = \frac{\partial}{\partial z} \left(K_{\psi} \frac{\partial\psi}{\partial z}\right) + \frac{\psi}{k} \left(c_1 P_s + c_3 P_b - c_2 \epsilon F_{\mathrm{wall}}\right).

The coefficients :math:`c_1` and :math:`c_2` are chosen to be consistent with observations of decaying homogeneous, isotropic turbulence. The other parameter :math:`c_3` has different values of stable (:math:`c_3^+`) and unstable (:math:`c_3^-`) stratification. Also,

.. math::
    \psi &= \left(c_{\mu}^0\right)^p k^m l^n\\
    l &= \left(c_{\mu}^0\right)^3 k^{3/2}\epsilon - 1

The indices :math:`p`, :math:`m`, and :math:`n` as well as coefficients are set in the :ref:`inputs file`<list-of-parameters-gls>`. The default values correspond to the :math:`k-\epsilon` turbulence model.

The equations for :math:`q` and :math:`\psi` are evolved much like the model tracer equations, including an implicit solve for vertical operations. We use a predictor-corrector scheme in which the predictor step only computes advection. The mixing coefficients are calculated from :math:`q` and :math:`l`:

.. math::
    K_v &= q l S_m + K_{v,\mathrm{background}}\\
    K_t &= q l S_h + K_{t,\mathrm{background}}\\
    K_k &= q l S_m / \sigma_k + K_{k,\mathrm{background}}\\
    K_{\psi} &= q l S_m / \sigma_{\psi} + K_{\psi,\mathrm{background}}

The constants :math:`\sigma_k = K_v / K_k` and :math:`\sigma_{\psi} = K_v / K_{\psi}` are also set in the inputs file. The stability coefficients, :math:`S_m` and :math:`S_h` are calculated using either the Galperin, Canuto-A or Canuto-B schemes.

