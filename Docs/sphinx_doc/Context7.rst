.. _Context7MCP:

Context7 for REMORA Workflows
=============================

`Context7 <https://context7.com>`_ can expose project documentation to
MCP-enabled coding assistants in real time. For REMORA work, this helps keep
parameter suggestions and API lookups grounded in current docs.

How It Fits
-----------

With MCP configured, an assistant can retrieve relevant sections from docs on
demand instead of relying only on model pretraining memory.

This is useful for:

- REMORA inputs editing and explanation,
- AMReX-level parameter and API questions,
- mixed C++/Python workflows that also touch pyAMReX.

Setup
-----

1. Configure Context7 in your client using:
   `Context7 client documentation <https://context7.com/docs/resources/all-clients>`_.
2. Add the documentation namespaces you want available to the assistant.
3. Restart/reload the client so MCP tools are active.

Recommended Namespaces
----------------------

For REMORA-focused usage, these two are the highest-value baseline:

- **AMReX**:
  `context7.com/amrex-codes/amrex <https://context7.com/amrex-codes/amrex>`__
- **pyAMReX**:
  `context7.com/amrex-codes/pyamrex <https://context7.com/amrex-codes/pyamrex>`__

If a REMORA-specific Context7 namespace is available in your environment, add
that as well.

Practical Usage Pattern
-----------------------

1. Start from ``Exec/Generic/inputs_generic`` for full parameter coverage.
2. Ask the assistant to cross-check a small set of keys against docs context.
3. Apply edits, then run a short REMORA test.
4. Iterate in small changes rather than rewriting the full file in one pass.

Example prompts:

- "Using Context7 docs, explain these parameters and suggest conservative
  values: ``remora.cfl``, ``amr.max_level``, ``geometry.prob_lo``."
- "Update this case from z-level style assumptions to REMORA sigma-coordinate
  expectations and list what changed."
