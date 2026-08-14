#!/usr/bin/env bash

set -eu -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

cd "${REPO_ROOT}"

INPUTS=Exec/Generic/inputs_generic

bash Exec/Generic/build_schema_remora.sh
python3 Exec/Generic/render_inputs_generic.py

if git diff --quiet -- "${INPUTS}"
then
    exit 0
fi

TMPDIR_CHECK="$(mktemp -d)"
trap 'rm -rf "${TMPDIR_CHECK}"' EXIT

# Source line numbers shift whenever unrelated code moves, so compare with them
# stripped: that isolates real parameter changes from incidental churn.
strip_source_lines() {
    sed -E 's|^(# source: .*):[0-9]+$|\1|'
}

param_keys() {
    grep -E '^[A-Za-z_][A-Za-z0-9_.]* = ' "${1}" | sed -E 's/ = .*//' | sort || true
}

git show ":${INPUTS}" | strip_source_lines > "${TMPDIR_CHECK}/committed"
strip_source_lines < "${INPUTS}" > "${TMPDIR_CHECK}/regenerated"

if diff -q "${TMPDIR_CHECK}/committed" "${TMPDIR_CHECK}/regenerated" > /dev/null
then
    if [ -n "${GITHUB_ACTIONS:-}" ]
    then
        echo "::warning title=inputs_generic source line numbers are stale::Parameters are unchanged; refresh with .github/workflows/style/check_generic_inputs.sh when convenient"
    fi
    echo -e "\nExec/Generic/inputs_generic records stale source line numbers, but no"
    echo -e "parameter was added, removed, or changed. Not treating this as a failure."
    echo -e "Refresh it whenever convenient with"
    echo -e "  ${0}\n"
    exit 0
fi

param_keys "${TMPDIR_CHECK}/committed" > "${TMPDIR_CHECK}/keys_committed"
param_keys "${TMPDIR_CHECK}/regenerated" > "${TMPDIR_CHECK}/keys_regenerated"
added="$(comm -13 "${TMPDIR_CHECK}/keys_committed" "${TMPDIR_CHECK}/keys_regenerated")"
removed="$(comm -23 "${TMPDIR_CHECK}/keys_committed" "${TMPDIR_CHECK}/keys_regenerated")"

if [ -n "${GITHUB_ACTIONS:-}" ]
then
    echo "::error title=inputs_generic is out of date::Run .github/workflows/style/check_generic_inputs.sh and commit Exec/Generic/inputs_generic"
fi

echo -e "\nExec/Generic/inputs_generic is out of date. Changes suggested by"
echo -e "  ${0}\n"

if [ -n "${added}" ]
then
    echo -e "Parameters missing from the committed file:"
    echo "${added}" | sed 's/^/  /'
    echo ""
fi

if [ -n "${removed}" ]
then
    echo -e "Parameters no longer found in the sources:"
    echo "${removed}" | sed 's/^/  /'
    echo ""
fi

if [ -z "${added}" ] && [ -z "${removed}" ]
then
    echo -e "The parameter list is unchanged, but values, types, or descriptions differ.\n"
fi

echo -e "Regenerate and commit the result with"
echo -e "  ${0}"
echo -e "  git add ${INPUTS}\n"
echo -e "That script runs, equivalently,"
echo -e "  bash Exec/Generic/build_schema_remora.sh"
echo -e "  python3 Exec/Generic/render_inputs_generic.py\n"
echo -e "Diff below ignores source line numbers; run 'git diff -- ${INPUTS}' for the full diff.\n"
diff -u --label "a/${INPUTS}" --label "b/${INPUTS}" \
    "${TMPDIR_CHECK}/committed" "${TMPDIR_CHECK}/regenerated" || true
echo ""
exit 1
