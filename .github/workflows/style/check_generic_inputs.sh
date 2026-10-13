#!/usr/bin/env bash

set -eu -o pipefail

# Exec/Generic/inputs_generic is regenerated on development by the
# generic-inputs-update workflow, which opens a bot PR. Feature branches should
# not commit it.
#
# Without arguments, regenerate the file in place and fail if it changed.
# With --preview (used on pull requests), report how the PR will change the
# file, restore the committed copy, and succeed: the bot commits the change
# after merge.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

cd "${REPO_ROOT}"

PREVIEW=0
if [ "${1:-}" = "--preview" ]
then
    PREVIEW=1
fi

INPUTS=Exec/Generic/inputs_generic

TMPDIR_CHECK="$(mktemp -d)"
cleanup() {
    if [ "${PREVIEW}" = 1 ]
    then
        git show ":${INPUTS}" > "${INPUTS}"
    fi
    rm -rf "${TMPDIR_CHECK}"
}
trap cleanup EXIT

bash Exec/Generic/build_schema_remora.sh
python3 Exec/Generic/render_inputs_generic.py

if git diff --quiet -- "${INPUTS}"
then
    exit 0
fi

param_keys() {
    grep -E '^[A-Za-z_][A-Za-z0-9_.]* = ' "${1}" | sed -E 's/ = .*//' | sort || true
}

git show ":${INPUTS}" > "${TMPDIR_CHECK}/committed"
cp "${INPUTS}" "${TMPDIR_CHECK}/regenerated"

param_keys "${TMPDIR_CHECK}/committed" > "${TMPDIR_CHECK}/keys_committed"
param_keys "${TMPDIR_CHECK}/regenerated" > "${TMPDIR_CHECK}/keys_regenerated"
added="$(comm -13 "${TMPDIR_CHECK}/keys_committed" "${TMPDIR_CHECK}/keys_regenerated")"
removed="$(comm -23 "${TMPDIR_CHECK}/keys_committed" "${TMPDIR_CHECK}/keys_regenerated")"

if [ "${PREVIEW}" = 1 ]
then
    if [ -n "${GITHUB_ACTIONS:-}" ]
    then
        echo "::warning title=inputs_generic will change after merge::No action needed; a bot PR updates Exec/Generic/inputs_generic once this merges. Do not commit it here."
    fi
    echo -e "\nMerging this PR will change Exec/Generic/inputs_generic as shown below."
    echo -e "No action is needed: a bot PR regenerates and commits it after merge."
    echo -e "Do not commit it here. To change a description, edit the comment at the"
    echo -e "query or its declaration, or the row in Docs/sphinx_doc/Inputs.rst.\n"
else
    if [ -n "${GITHUB_ACTIONS:-}" ]
    then
        echo "::error title=inputs_generic is out of date::Run .github/workflows/style/check_generic_inputs.sh and commit Exec/Generic/inputs_generic"
    fi
    echo -e "\nExec/Generic/inputs_generic is out of date. Changes suggested by"
    echo -e "  ${0}\n"
fi

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

if [ "${PREVIEW}" = 1 ]
then
    diff -u --label "a/${INPUTS}" --label "b/${INPUTS}" \
        "${TMPDIR_CHECK}/committed" "${TMPDIR_CHECK}/regenerated" || true
    echo ""
    exit 0
fi

echo -e "Regenerate and commit the result with"
echo -e "  ${0}"
echo -e "  git add ${INPUTS}\n"
echo -e "That script runs, equivalently,"
echo -e "  bash Exec/Generic/build_schema_remora.sh"
echo -e "  python3 Exec/Generic/render_inputs_generic.py\n"
diff -u --label "a/${INPUTS}" --label "b/${INPUTS}" \
    "${TMPDIR_CHECK}/committed" "${TMPDIR_CHECK}/regenerated" || true
echo ""
exit 1
