.. _AgenticWorkFlow:

Agentic Workflow
================

`AMReX-Agent <https://github.com/AMReX-Codes/amrex-agent>`_ can be used to
draft and iterate on REMORA input files with schema-aware context. This page
keeps the setup close to the existing REMORA generic-inputs flow in this
repository.

For MCP/Context7 setup in coding assistants, see :ref:`Context7MCP`.

Basic Setup
-----------

1. Clone ``amrex-agent`` and enter the repository.

   .. code-block:: bash

      git clone --recursive https://github.com/AMReX-Codes/amrex-agent.git
      cd amrex-agent

2. Create and activate the demo environment.

   .. code-block:: bash

      conda env create -f environment.yaml
      conda activate amrex-agent-dev

3. Export paths to local code checkouts.

   .. code-block:: bash

      export REMORA_REPO_PATH=/path/to/REMORA
      export AMREX_REPO_PATH=/path/to/amrex

4. Build REMORA schema and indices for the agent.

   .. code-block:: bash

      bash demo/setup_demo_database.sh --code remora --force-rebuild

   If you are testing pipeline mechanics without live embeddings/API access:

   .. code-block:: bash

      bash demo/setup_demo_database.sh --code remora --force-rebuild --mock

5. Run a REMORA smoke-test prompt.

   .. code-block:: bash

      python amrex_agent.py \
        --baseline-override "REMORA/Exec/Upwelling" \
        --indexing-strategy override_static \
        --prompt "Run the REMORA Upwelling case to demonstrate wind-driven upwelling over a periodic channel."

Using ``inputs_generic`` in the Loop
------------------------------------

``Exec/Generic/inputs_generic`` is a generated parameter catalog that works
well as a starting reference while iterating with an agent:

1. draft case-specific inputs from ``inputs_generic`` sections,
2. ask the agent to explain or edit targeted parameter groups,
3. run short REMORA tests and adjust.

When REMORA ``ParmParse`` usage changes, regenerate the catalog:

.. code-block:: bash

   bash Exec/Generic/build_schema_remora.sh
   python3 Exec/Generic/render_inputs_generic.py

Review Rhythm
-------------

- Keep generated-file updates and source changes in the same review context.
- Treat AI edits as draft patches and verify behavior on a short run before
  scaling up.
