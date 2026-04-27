Agentic Workflow
================

This page describes a lightweight AI-assisted workflow for REMORA inputs work.

Recommended Loop
----------------

1. Start from ``Exec/Generic/inputs_generic`` as the parameter catalog.
2. Draft or update your case input file.
3. Validate parameter intent against ``Docs/sphinx_doc/Inputs.rst``.
4. If REMORA source changes, regenerate schema/template before review:

   .. code-block:: bash

      bash Exec/Generic/build_schema_remora.sh
      python3 Exec/Generic/render_inputs_generic.py

Usage Notes
-----------

- Treat generated inputs as a starting point and adjust for each case.
- Reviewing ``inputs_generic`` changes alongside source edits helps keep context clear.
- Keep generated files in sync; CI checks for stale output.

.. note::

   LLM-generated code and docs are draft edits; a quick verification of intent
   and behavior helps catch drift.

Notes
-----

This workflow is validate-only in CI. It does not auto-commit generated files.
