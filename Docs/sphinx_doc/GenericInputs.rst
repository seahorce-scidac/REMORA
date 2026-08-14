Generic Master Inputs
=====================

REMORA now ships a generated master template at ``Exec/Generic/inputs_generic``.
This file is intended as a comprehensive, ROMS-style reference inputs file that
lists all currently discovered ``ParmParse`` parameters.

This is reference material for case setup and parameter discovery, not a fixed
workflow contract.

Generation Flow
---------------

1. Build schema from REMORA + dependency sources:

   .. code-block:: bash

      bash Exec/Generic/build_schema_remora.sh

2. Render the generic inputs template:

   .. code-block:: bash

      python3 Exec/Generic/render_inputs_generic.py

Value Selection
---------------

For each parameter, values are selected with deterministic precedence:

1. schema default extracted from source,
2. ``Exec/Generic/remora_value_hints.json`` value hint,
3. placeholder ``<REQUIRED>`` or ``<UNKNOWN_DEFAULT>``.

Comment Selection
-----------------

For each parameter, comments are selected in this order:

1. source-adjacent description extracted with the schema,
2. table description from ``Docs/sphinx_doc/Inputs.rst``,
3. source location fallback.

CI Validation
-------------

GitHub Actions workflow ``generic-inputs-validate.yml`` regenerates
``Exec/Generic/inputs_generic`` and reports when the file is stale. The check
itself is ``.github/workflows/style/check_generic_inputs.sh``, which prints the
stale diff along with the command that fixes it. To clear the failure, run that
same script locally and commit the result:

.. code-block:: bash

   .github/workflows/style/check_generic_inputs.sh
   git add Exec/Generic/inputs_generic

Because the schema records source file line numbers, a branch that is behind
``development`` will report a stale file even with no local source changes.
Merge ``development`` first, then regenerate.
