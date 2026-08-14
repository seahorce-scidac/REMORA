#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

python3 "${SCRIPT_DIR}/build_schema_remora.py" \
  --repo-root "${REPO_ROOT}" \
  --output-dir "${SCRIPT_DIR}" \
  "$@"
