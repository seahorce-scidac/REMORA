#!/usr/bin/env bash

set -eu -o pipefail

# Exec/Generic/inputs_generic is regenerated on development by the
# generic-inputs-update workflow, which opens a bot PR. Feature branches should
# not commit it.
#
# Usage:
#   check_generic_inputs.sh
#       Regenerate the file in place and fail if it changed.
#   check_generic_inputs.sh --preview BASE
#       Report how the changes since BASE (for a branch, the branch it merges
#       into, e.g. upstream/development) change the generated file, and
#       succeed.
#       Both sides are regenerated, each with its own copy of the generator,
#       so the report shows only this change's effect even when the committed
#       file is stale. The working tree is left untouched.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

cd "${REPO_ROOT}"

usage() {
    echo "usage: ${0} [--preview BASE]" >&2
    exit 2
}

PREVIEW=0
BASE=
case "$#:${1:-}" in
    0:) ;;
    2:--preview) PREVIEW=1; BASE="${2}" ;;
    *) usage ;;
esac

INPUTS=Exec/Generic/inputs_generic
AMREX=Submodules/AMReX

TMPDIR_CHECK="$(mktemp -d)"
BASE_TREE="${TMPDIR_CHECK}/base"
cleanup() {
    if [ -d "${BASE_TREE}" ]
    then
        git worktree remove --force "${BASE_TREE}" || true
    fi
    rm -rf "${TMPDIR_CHECK}"
}
trap cleanup EXIT

# Render the file for the source tree at $1 into $2, using that tree's own
# generator. Schemas go to the tree's gitignored Exec/Generic/.generated.
render_tree() {
    local root="${1}"
    bash "${root}/Exec/Generic/build_schema_remora.sh" > /dev/null
    python3 "${root}/Exec/Generic/render_inputs_generic.py" \
        --repo-root "${root}" --output "${2}" > /dev/null
}

param_keys() {
    grep -E '^[A-Za-z_][A-Za-z0-9_.]* = ' "${1}" | sed -E 's/ = .*//' | sort || true
}

if [ "${PREVIEW}" = 1 ]
then
    base_sha="$(git merge-base "${BASE}" HEAD)"
    git worktree add --quiet --detach "${BASE_TREE}" "${base_sha}"
    if [ "$(git rev-parse "${base_sha}:${AMREX}")" = "$(git rev-parse "HEAD:${AMREX}")" ]
    then
        rmdir "${BASE_TREE}/${AMREX}"
        ln -s "${REPO_ROOT}/${AMREX}" "${BASE_TREE}/${AMREX}"
    elif ! git -C "${BASE_TREE}" submodule update --init --depth 1 "${AMREX}"
    then
        # The preview is advisory; do not fail a PR because it is unavailable.
        if [ -n "${GITHUB_ACTIONS:-}" ]
        then
            echo "::warning title=inputs_generic preview unavailable::Could not fetch the base AMReX commit to compare against."
        fi
        echo "Could not fetch the base ${AMREX}; skipping the preview." >&2
        exit 0
    fi
    render_tree "${BASE_TREE}" "${TMPDIR_CHECK}/committed"
    render_tree "${REPO_ROOT}" "${TMPDIR_CHECK}/regenerated"
else
    bash Exec/Generic/build_schema_remora.sh
    python3 Exec/Generic/render_inputs_generic.py
    git show ":${INPUTS}" > "${TMPDIR_CHECK}/committed"
    cp "${INPUTS}" "${TMPDIR_CHECK}/regenerated"
fi

if cmp -s "${TMPDIR_CHECK}/committed" "${TMPDIR_CHECK}/regenerated"
then
    exit 0
fi

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
    echo -e "\nThese changes alter the generated Exec/Generic/inputs_generic as shown below."
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
    echo -e "Parameters added:"
    echo "${added}" | sed 's/^/  /'
    echo ""
fi

if [ -n "${removed}" ]
then
    echo -e "Parameters removed:"
    echo "${removed}" | sed 's/^/  /'
    echo ""
fi

if [ -z "${added}" ] && [ -z "${removed}" ]
then
    echo -e "The parameter list is unchanged, but values, types, or descriptions differ.\n"
fi

if [ "${PREVIEW}" = 0 ]
then
    echo -e "Regenerate and commit the result with"
    echo -e "  ${0}"
    echo -e "  git add ${INPUTS}\n"
    echo -e "That script runs, equivalently,"
    echo -e "  bash Exec/Generic/build_schema_remora.sh"
    echo -e "  python3 Exec/Generic/render_inputs_generic.py\n"
fi

diff -u --label "a/${INPUTS}" --label "b/${INPUTS}" \
    "${TMPDIR_CHECK}/committed" "${TMPDIR_CHECK}/regenerated" || true
echo ""

if [ "${PREVIEW}" = 1 ]
then
    exit 0
fi
exit 1
