Agentic Workflow
================

This page describes a practical AI-assisted workflow for REMORA inputs work.

Recommended Loop
----------------

1. Start from ``Exec/Generic/inputs_generic`` as the canonical parameter catalog.
2. Draft or update your case input file.
3. Validate parameter intent against ``Docs/sphinx_doc/Inputs.rst``.
4. If REMORA source changes, regenerate schema/template before review:

   .. code-block:: bash

      bash Exec/Generic/build_schema_remora.sh
      python3 Exec/Generic/render_inputs_generic.py

Review and Safety Guidance
--------------------------

- Treat generated inputs as a starting template, not as an automatically correct run configuration.
- Validate physics choices, boundary conditions, and forcing assumptions with domain experts.
- Review changes in ``Exec/Generic/inputs_generic`` with normal PR review process.
- Keep generated files in sync; CI will fail if the checked-in template is stale.

Notes
-----

This workflow is intentionally validate-only in CI. It does not auto-commit generated files.
