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

1. source-adjacent description extracted with the schema: a trailing comment on
   the query line, a comment block directly above it (no blank or code line in
   between), or the same for the variable's declaration,
2. table description from ``Docs/sphinx_doc/Inputs.rst``,
3. source file fallback.

The source fallback names only the file, not a line number, so that unrelated
code motion does not change the generated file. A comment that is not
directly attached to a query is never borrowed for it.

CI Validation and Updates
-------------------------

Do not commit ``Exec/Generic/inputs_generic`` from a feature branch. After a PR
merges into ``development``, the ``generic-inputs-update.yml`` workflow
regenerates the file and opens (or updates) a single PR from the
``bot/inputs-generic`` branch. A maintainer reviews and merges it. Because
feature branches never touch the file, it does not conflict when merging
``development``.

On pull requests, ``generic-inputs-validate.yml``:

* fails if the PR changes ``inputs_generic``. Restore it with

  .. code-block:: bash

     git checkout origin/development -- Exec/Generic/inputs_generic

  To edit it deliberately, for example while changing the generator, label the
  PR ``inputs-generic-hand-edit``;
* runs ``.github/workflows/style/check_generic_inputs.sh --preview``, which
  warns, without failing, with the diff the merge will cause. Use it to check
  that new parameters and descriptions come out as intended. To change a
  description, edit the comment at the query or its declaration, or the row in
  ``Docs/sphinx_doc/Inputs.rst``.

The preview regenerates the file both at the PR's base and at its head, so it
shows only the PR's own effect, even while an earlier bot PR is unmerged. It
writes nothing to the working tree. To preview a local branch, pass the branch
it will merge into:

.. code-block:: bash

   .github/workflows/style/check_generic_inputs.sh --preview upstream/development

Without ``--preview``, the script regenerates the file in place and fails if it
changed; a manual run of the workflow does this on the chosen branch.

The file changes when a parameter is added or removed; when its default, type,
or attached description changes; when its row in ``Inputs.rst`` changes; or
when its query moves to another file. Moving code within a file normally does
not change it. The exception is a description taken from the variable's
declaration, which is looked for only within 140 lines above the query.
