#!/usr/bin/env python3
"""Render a ROMS-style master REMORA inputs template from schema + local hints."""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parents[1]
    parser = argparse.ArgumentParser(description="Render Exec/Generic/inputs_generic from schema.")
    parser.add_argument(
        "--schema",
        type=Path,
        default=script_dir / ".generated" / "remora_schema_current.json",
        help="Schema JSON path",
    )
    parser.add_argument(
        "--value-hints",
        type=Path,
        default=script_dir / "remora_value_hints.json",
        help="Repo-local value hints JSON",
    )
    parser.add_argument(
        "--docs-inputs",
        type=Path,
        default=repo_root / "Docs" / "sphinx_doc" / "Inputs.rst",
        help="Docs Inputs.rst path",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=script_dir / "inputs_generic",
        help="Output inputs file path",
    )
    return parser.parse_args()


def load_schema(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        payload = json.load(fh)
    if "parameters" not in payload or not isinstance(payload["parameters"], dict):
        raise ValueError(f"Schema missing 'parameters' object: {path}")
    return payload


def load_value_hints(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as fh:
        payload = json.load(fh)

    if isinstance(payload, dict) and "hints" in payload and isinstance(payload["hints"], dict):
        raw_hints = payload["hints"]
    elif isinstance(payload, dict):
        raw_hints = payload
    else:
        return {}

    hints: dict[str, Any] = {}
    for key, val in raw_hints.items():
        if isinstance(val, dict) and "value" in val:
            hints[key] = val["value"]
        else:
            hints[key] = val
    return hints


def parse_docs_descriptions(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}

    row_re = re.compile(r"^\|\s+\*\*([^*]+)\*\*\s*\|\s*([^|]+?)\s*\|")
    descs: dict[str, str] = {}

    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        m = row_re.match(line)
        if not m:
            continue
        param = m.group(1).strip()
        desc = m.group(2).strip()
        if not desc:
            continue
        # Keep first appearance for stability.
        descs.setdefault(param, desc)

    return descs


def section_for_param(name: str) -> str:
    if "." in name:
        return name.split(".", 1)[0]
    return "<noprefix>"


def ordered_sections(keys: list[str]) -> list[str]:
    preferred = [
        "remora",
        "amr",
        "geometry",
        "fabarray",
        "vismf",
        "amrex",
        "eb2",
        "integration",
        "particles",
        "device",
        "tiny_profiler",
        "DistributionMapping",
        "plot",
        "fab",
        "<noprefix>",
    ]

    present = {section_for_param(k) for k in keys}
    out: list[str] = [p for p in preferred if p in present]
    leftovers = sorted(present - set(out))
    out.extend(leftovers)
    return out


def section_label(section: str) -> str:
    labels = {
        "remora": "REMORA PARAMETERS",
        "amr": "AMR PARAMETERS",
        "geometry": "GEOMETRY PARAMETERS",
        "fabarray": "FABARRAY PARAMETERS",
        "vismf": "VISMF PARAMETERS",
        "amrex": "AMREX BASE PARAMETERS",
        "eb2": "EB2 PARAMETERS",
        "integration": "INTEGRATION PARAMETERS",
        "particles": "PARTICLE PARAMETERS",
        "device": "DEVICE PARAMETERS",
        "tiny_profiler": "TINY PROFILER PARAMETERS",
        "DistributionMapping": "DISTRIBUTION MAPPING PARAMETERS",
        "plot": "PLOT TOOL PARAMETERS",
        "fab": "FAB PARAMETERS",
        "<noprefix>": "LEGACY / UNPREFIXED PARAMETERS",
    }
    return labels.get(section, f"{section.upper()} PARAMETERS")


def choose_value(name: str, param: dict[str, Any], hints: dict[str, Any]) -> Any:
    if param.get("default") is not None:
        return param["default"]
    if name in hints:
        return hints[name]
    return "<REQUIRED>" if bool(param.get("required")) else "<UNKNOWN_DEFAULT>"


def format_value(value: Any, is_array: bool) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    if isinstance(value, list):
        return " ".join(format_value(v, is_array=False) for v in value)

    s = str(value)
    if is_array:
        return s

    # Keep placeholders and simple identifiers unquoted.
    if s.startswith("<") and s.endswith(">"):
        return s
    if re.fullmatch(r"[A-Za-z0-9_./:+-]+", s):
        return s

    return f'"{s}"'


def choose_comment(name: str, param: dict[str, Any], docs_desc: dict[str, str]) -> str:
    if param.get("description"):
        return str(param["description"]).strip()
    if name in docs_desc:
        return docs_desc[name]
    src_file = param.get("source_file")
    src_line = param.get("source_line")
    if src_file and src_line:
        return f"source: {src_file}:{src_line}"
    if src_file:
        return f"source: {src_file}"
    return ""


def render_text(schema: dict[str, Any], hints: dict[str, Any], docs_desc: dict[str, str], schema_path: Path) -> str:
    params: dict[str, dict[str, Any]] = schema["parameters"]
    keys = sorted(params)

    by_section: dict[str, list[str]] = defaultdict(list)
    for key in keys:
        by_section[section_for_param(key)].append(key)

    lines: list[str] = []
    lines.append("#")
    lines.append("# REMORA Generic Master Inputs (auto-generated)")
    lines.append("#")
    lines.append("# Regenerate with:")
    lines.append("#   bash Exec/Generic/build_schema_remora.sh")
    lines.append("#   python3 Exec/Generic/render_inputs_generic.py")
    lines.append(f"# Schema: {schema_path}")
    lines.append("#")
    lines.append("# Value precedence:")
    lines.append("#   1) schema default")
    lines.append("#   2) Exec/Generic/remora_value_hints.json")
    lines.append("#   3) <REQUIRED> / <UNKNOWN_DEFAULT>")
    lines.append("#")

    for section in ordered_sections(keys):
        lines.append("#==============================================================================")
        lines.append(f"# {section_label(section)}")
        lines.append("#==============================================================================")
        for key in sorted(by_section[section]):
            param = params[key]
            comment = choose_comment(key, param, docs_desc)
            if comment:
                lines.append(f"# {comment}")

            ptype = param.get("type") or "unknown"
            method = param.get("method") or "unknown"
            lines.append(f"# type={ptype}; method={method}")

            value = choose_value(key, param, hints)
            rendered_value = format_value(value, bool(param.get("is_array")))
            lines.append(f"{key} = {rendered_value}")
            lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def main() -> int:
    args = parse_args()

    schema = load_schema(args.schema)
    hints = load_value_hints(args.value_hints)
    docs_desc = parse_docs_descriptions(args.docs_inputs)

    text = render_text(schema, hints, docs_desc, args.schema)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")

    print(f"Wrote: {args.output}")
    print(f"Parameters: {len(schema['parameters'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
