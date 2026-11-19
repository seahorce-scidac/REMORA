<!--
Target `development`. PRs are squash-merged, so the title becomes the commit summary:
keep it short and imperative ("Read only each rank's part of a lat-lon forcing record").
One topic per PR; put unrelated fixes in their own PRs. See CONTRIBUTING.md.
Comments like this one are hidden in the rendered PR and can be left in place.
-->

## Summary

<!-- The problem first, then what this PR changes. For a ROMS port or a change that should
match ROMS, name the ROMS routine and list any intentional differences. -->

-

Closes #

## Testing

<!-- Build configuration (CMake or GNU make; MPI/GPU; PnetCDF or not; precision), the cases
run, and results with numbers. Say what was verified by hand and what was not run. -->

- Build:
- Cases:
- Results:

## Answer changes

<!-- Delete the line that does not apply. -->

- [ ] No answer changes: regression tests pass unchanged.
- [ ] Answers change: which tests and fields change, by how much, and why. Gold files
      updated:

## Checklist

- [ ] New or changed behavior has a test in `Tests/CTestList.cmake` that fails when the
      feature is broken.
- [ ] New options are in the `remora.` namespace, default to the ROMS behavior, and are
      documented in `Docs/sphinx_doc/`.
- [ ] New tests are listed in `Docs/sphinx_doc/testing.rst` / `RegressionTests.rst`.
- [ ] Builds with and without PnetCDF (if NetCDF code was touched) and with
      `PRECISION=FLOAT` (if particle or precision-sensitive code was touched).
- [ ] `Exec/Generic/inputs_generic` is not modified (a bot regenerates it after merge).
- [ ] Follow-up work or deferred findings have their own issues:
