Generic Master Inputs
=====================

REMORA now ships a generated master template at ``Exec/Generic/inputs_generic``.
This file is intended as a comprehensive, ROMS-style reference inputs file that
lists all currently discovered ``ParmParse`` parameters.

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
``Exec/Generic/inputs_generic`` and fails if the file is stale.
