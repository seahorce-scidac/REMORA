#!/usr/bin/env python3
"""Build and compose a REMORA schema from ParmParse usage.

This script is repo-local and CI-safe:
- discovers dependencies from .gitmodules (currently Submodules/AMReX)
- scans ParmParse usage with broad AMReX method coverage
- composes dependency schemas + REMORA schema deterministically
- writes versioned schema JSON and updates remora_schema_current.json symlink
"""

from __future__ import annotations

import argparse
import contextlib
import json
import re
import subprocess
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

CPP_EXTS = (".cpp", ".cc", ".cxx", ".H", ".h", ".hpp")

UNNUMBERED_112_ID = "UNNUMBERED-112"
TIER2_STABILITY_PARAMETERS = ("amr.max_level", "amr.blocking_factor")

# Broad method set from AMReX ParmParse API.
METHODS = (
    "get",
    "query",
    "getarr",
    "queryarr",
    "getkth",
    "querykth",
    "getktharr",
    "queryktharr",
    "getline",
    "queryline",
    "queryAdd",
    "queryAddWithParser",
    "getWithParser",
    "queryWithParser",
    "getarrWithParser",
    "queryarrWithParser",
    "getAsDouble",
    "queryAsDouble",
    "getarrAsDouble",
    "queryarrAsDouble",
    "gettable",
    "querytable",
    "get_enum_case_insensitive",
    "query_enum_case_insensitive",
    "get_enum_sloppy",
    "query_enum_sloppy",
    "add",
    "addarr",
    "set",
)

REQUIRED_METHODS = {
    "get",
    "getarr",
    "getkth",
    "getktharr",
    "getline",
    "getWithParser",
    "getarrWithParser",
    "getAsDouble",
    "getarrAsDouble",
    "gettable",
    "get_enum_case_insensitive",
    "get_enum_sloppy",
    "add",
    "addarr",
}

ARRAY_METHODS = {
    "getarr",
    "queryarr",
    "getkth",
    "querykth",
    "getktharr",
    "queryktharr",
    "getline",
    "queryline",
    "queryAdd",
    "queryAddWithParser",
    "addarr",
    "getarrWithParser",
    "queryarrWithParser",
    "getarrAsDouble",
    "queryarrAsDouble",
    "gettable",
    "querytable",
}

TYPE_PATTERNS = (
    (r"\b(amrex::)?RealVect\b", "RealVect"),
    (r"\b(amrex::)?IntVect\b", "IntVect"),
    (r"\b(amrex::)?Real\b", "Real"),
    (r"\b(double|float)\b", "Real"),
    (r"\b(long\s+int|long|int)\b", "int"),
    (r"\b(bool)\b", "bool"),
    (r"\b(std::)?string\b", "string"),
    (r"\b(Vector|std::vector)\s*<", "Vector"),
)

ALLOWED_PARAM_NAME = re.compile(r"^[a-zA-Z0-9_.]+$")


@dataclass
class Declaration:
    var_name: str
    var_type: str
    file_path: Path
    line_num: int


@dataclass
class Dependency:
    name: str
    path: Path
    commit: str


@dataclass
class SolverConfig:
    schema_source_patterns: list[str]
    build_requirements: dict[str, list[str]]
    manual_schema_params: dict[str, dict[str, Any]]
    schema_version: int = 1


def sh(cmd: list[str], cwd: Path | None = None) -> str:
    out = subprocess.check_output(
        cmd,
        cwd=str(cwd) if cwd else None,
        text=True,
        stderr=subprocess.DEVNULL,
    )
    return out.strip()


def short_commit(repo_root: Path) -> str:
    try:
        return sh(["git", "-C", str(repo_root), "rev-parse", "--short", "HEAD"])
    except Exception:
        return "unknown"


def full_commit(repo_root: Path) -> str:
    try:
        return sh(["git", "-C", str(repo_root), "rev-parse", "HEAD"])
    except Exception:
        return "unknown"


def normalize_type(raw: str | None) -> str | None:
    if raw is None:
        return None
    cleaned = " ".join(raw.replace("\t", " ").split())
    cleaned = cleaned.replace("amrex::", "")
    mapping = {
        "double": "Real",
        "float": "Real",
        "Real": "Real",
        "int": "int",
        "long": "int",
        "long int": "int",
        "bool": "bool",
        "std::string": "string",
        "string": "string",
        "IntVect": "IntVect",
        "RealVect": "RealVect",
    }
    return mapping.get(cleaned, cleaned)


def parse_literal(expr: str) -> Any:
    expr = expr.strip()
    expr = expr.split("//", 1)[0].strip()
    expr = expr.replace("_rt", "").strip()

    if not expr:
        return None

    if expr.lower() in {"true", "false"}:
        return expr.lower() == "true"

    if len(expr) >= 2 and expr[0] == '"' and expr[-1] == '"':
        return expr[1:-1]

    if expr.startswith("{") and expr.endswith("}"):
        inner = expr[1:-1].strip()
        if not inner:
            return []
        vals = []
        for tok in split_top_level(inner):
            p = parse_literal(tok)
            vals.append(tok.strip() if p is None else p)
        return vals

    if re.fullmatch(r"[+-]?\d+", expr):
        with contextlib.suppress(ValueError):
            return int(expr)

    if re.fullmatch(r"[+-]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[eE][+-]?\d+)?", expr):
        with contextlib.suppress(ValueError):
            return float(expr)

    return None


def split_top_level(arg_text: str) -> list[str]:
    out: list[str] = []
    cur: list[str] = []
    depth_paren = 0
    depth_bracket = 0
    depth_brace = 0
    in_str = False
    esc = False

    for ch in arg_text:
        if in_str:
            cur.append(ch)
            if esc:
                esc = False
            elif ch == "\\":
                esc = True
            elif ch == '"':
                in_str = False
            continue

        if ch == '"':
            in_str = True
            cur.append(ch)
            continue

        if ch == "(":
            depth_paren += 1
        elif ch == ")":
            depth_paren = max(0, depth_paren - 1)
        elif ch == "[":
            depth_bracket += 1
        elif ch == "]":
            depth_bracket = max(0, depth_bracket - 1)
        elif ch == "{":
            depth_brace += 1
        elif ch == "}":
            depth_brace = max(0, depth_brace - 1)

        if ch == "," and depth_paren == 0 and depth_bracket == 0 and depth_brace == 0:
            out.append("".join(cur).strip())
            cur = []
            continue

        cur.append(ch)

    if cur:
        out.append("".join(cur).strip())

    return out


def is_valid_param_name(name: str) -> bool:
    if not name or not name.strip():
        return False
    name = name.strip()
    if name in {".", ".."}:
        return False
    if name.endswith("."):
        return False
    return bool(ALLOWED_PARAM_NAME.match(name))


def parse_ifdefs(lines: list[str]) -> dict[int, list[str]]:
    active: list[str] = []
    blocks: dict[int, list[str]] = {}

    for idx, line in enumerate(lines):
        stripped = line.strip()

        m = re.match(r"#ifdef\s+(\w+)", stripped)
        if m:
            macro = m.group(1)
            if macro.startswith("AMREX_"):
                macro = macro[6:]
            active.append(macro)

        m = re.match(r"#ifndef\s+(\w+)", stripped)
        if m:
            macro = "!" + m.group(1)
            if macro.startswith("!AMREX_"):
                macro = "!" + macro[7:]
            active.append(macro)

        m = re.match(r"#if\s+defined\((\w+)\)", stripped)
        if m:
            macro = m.group(1)
            if macro.startswith("AMREX_"):
                macro = macro[6:]
            active.append(macro)

        if stripped.startswith("#endif") and active:
            active.pop()

        if active:
            blocks[idx] = list(active)

    return blocks


def extract_description(lines: list[str], line_idx: int) -> str | None:
    if line_idx < 0 or line_idx >= len(lines):
        return None

    line = lines[line_idx]
    if "//" in line:
        tail = line.split("//", 1)[1].strip()
        if tail:
            return tail

    comments: list[str] = []
    for j in range(line_idx - 1, max(-1, line_idx - 6), -1):
        s = lines[j].strip()
        if not s:
            if comments:
                break
            continue
        if s.startswith("//"):
            comments.insert(0, s[2:].strip())
            continue
        if s.startswith("/*") or s.startswith("*") or s.endswith("*/"):
            cleaned = s.replace("/*", "").replace("*/", "").lstrip("*").strip()
            if cleaned:
                comments.insert(0, cleaned)
            continue
        if comments:
            break

    text = " ".join(c for c in comments if c)
    return text or None


def extract_var_description(var_name: str, lines: list[str], line_idx: int) -> str | None:
    if not var_name:
        return None
    decl_pat = re.compile(rf"\b{re.escape(var_name)}\b")
    for j in range(line_idx - 1, max(-1, line_idx - 140), -1):
        s = lines[j]
        if not decl_pat.search(s):
            continue
        if "//" in s:
            tail = s.split("//", 1)[1].strip()
            if tail:
                return tail
        if "/*" in s and "*/" in s:
            inner = s.split("/*", 1)[1].split("*/", 1)[0].strip()
            if inner:
                return inner
    return None


def detect_type_inline(var_name: str, lines: list[str], line_idx: int, full_text: str) -> str | None:
    if not var_name:
        return None

    for j in range(line_idx - 1, max(-1, line_idx - 140), -1):
        s = lines[j]
        for pattern, norm in TYPE_PATTERNS:
            if re.search(rf"{pattern}\s+[\*&]?\s*{re.escape(var_name)}\b", s):
                return normalize_type(norm)

    for pattern, norm in TYPE_PATTERNS:
        if re.search(rf"{pattern}\s+[\*&]?\s*{re.escape(var_name)}\b", full_text):
            return normalize_type(norm)

    return None


def extract_default(var_name: str, lines: list[str], line_idx: int) -> Any:
    if not var_name:
        return None

    assign_pat = re.compile(rf"\b{re.escape(var_name)}\s*=\s*([^;]+);")
    for j in range(line_idx - 1, max(-1, line_idx - 180), -1):
        m = assign_pat.search(lines[j])
        if not m:
            continue
        parsed = parse_literal(m.group(1))
        if parsed is not None:
            return parsed

    return None


class DeclarationExtractor:
    """Extract C++ variable declarations for better type inference."""

    @classmethod
    def extract_all(cls, files: list[Path]) -> dict[str, list[Declaration]]:
        out: dict[str, list[Declaration]] = defaultdict(list)
        for file_path in files:
            with contextlib.suppress(Exception):
                decls = cls._extract_from_file(file_path)
                for d in decls:
                    out[d.var_name].append(d)
        return out

    @classmethod
    def _extract_from_file(cls, file_path: Path) -> list[Declaration]:
        text = file_path.read_text(encoding="utf-8", errors="ignore")
        lines = text.splitlines()
        decls: list[Declaration] = []

        decl_re = re.compile(
            r"^\s*(?:static\s+|const\s+|constexpr\s+|mutable\s+|inline\s+|extern\s+)*"
            r"([A-Za-z_][\w:<>]*(?:\s*<[^;>{}()]+>)?(?:\s+[A-Za-z_][\w:<>]*)?)"
            r"\s+[\*&]?\s*([A-Za-z_]\w*)\s*(?:[=;,{[])"
        )

        for i, line in enumerate(lines, start=1):
            m = decl_re.match(line)
            if not m:
                continue
            raw_type = normalize_type(m.group(1).strip())
            var_name = m.group(2)
            if not raw_type:
                continue
            decls.append(Declaration(var_name=var_name, var_type=raw_type, file_path=file_path, line_num=i))

        return decls

    @staticmethod
    def match(var_name: str, file_path: Path, decl_map: dict[str, list[Declaration]]) -> str | None:
        if var_name not in decl_map:
            return None
        decls = decl_map[var_name]

        for d in decls:
            if d.file_path == file_path:
                return d.var_type

        headers = [d for d in decls if d.file_path.suffix in {".H", ".h", ".hpp"}]
        if headers:
            return headers[0].var_type

        return decls[0].var_type if decls else None


def resolve_string_expr(expr: str, string_vars: dict[str, str]) -> str | None:
    expr = expr.strip()
    if not expr:
        return ""

    # Raw literal.
    m = re.fullmatch(r'"([^"]*)"', expr)
    if m:
        return m.group(1)

    # Braced literal {"foo"}
    m = re.fullmatch(r'\{\s*"([^"]*)"\s*\}', expr)
    if m:
        return m.group(1)

    # Simple variable.
    if re.fullmatch(r"[A-Za-z_]\w*", expr):
        return string_vars.get(expr)

    # Concatenation of literals/vars.
    parts = [p.strip() for p in expr.split("+")]
    resolved: list[str] = []
    any_resolved = False
    for p in parts:
        m = re.fullmatch(r'"([^"]*)"', p)
        if m:
            resolved.append(m.group(1))
            any_resolved = True
            continue
        if re.fullmatch(r"[A-Za-z_]\w*", p):
            if p in string_vars:
                resolved.append(string_vars[p])
                any_resolved = True
            else:
                resolved.append(f"<{p}>")
            continue
        m = re.fullmatch(r"([A-Za-z_]\w*)\s*\[[^\]]+\]", p)
        if m:
            resolved.append(f"<{m.group(1)}[]>")
            continue
        # Unknown token.
        resolved.append(f"<{p.replace(' ', '')}>")

    if not resolved or not any_resolved:
        return None

    out = "".join(resolved)
    out = out.strip()
    out = out[:-1] if out.endswith(".") else out
    return out


def source_type_for_file(file_path: Path, root: Path) -> str:
    try:
        rel = str(file_path.relative_to(root))
    except ValueError:
        rel = str(file_path)
    if "_cpp_parameters" in rel or "generated" in rel.lower():
        return "generated"
    if file_path.suffix in {".H", ".h", ".hpp"}:
        return "header"
    return "manual"


def param_priority(name: str) -> str:
    tier1 = {"remora.cfl", "amr.n_cell", "remora.plot_int"}
    tier2 = set(TIER2_STABILITY_PARAMETERS)
    if name in tier1:
        return "tier1"
    if name in tier2:
        return "tier2"
    return "tier3"


class SchemaBuilder:
    def __init__(self, root: Path, component_name: str):
        self.root = root
        self.component_name = component_name
        self.schema: dict[str, dict[str, Any]] = {}
        self.declaration_map: dict[str, list[Declaration]] = {}

    def _iter_source_files(self, source_dirs: list[str]) -> list[Path]:
        files: list[Path] = []
        for rel in source_dirs:
            base = self.root / rel
            if not base.exists() or not base.is_dir():
                continue
            for ext in CPP_EXTS:
                files.extend(base.rglob(f"*{ext}"))
        return sorted(set(files))

    def scan_source_code(self, source_dirs: list[str], config: SolverConfig) -> dict[str, dict[str, Any]]:
        # Phase 0: _cpp_parameters parsing, if present.
        for rel in source_dirs:
            candidate = self.root / rel / "_cpp_parameters"
            if candidate.exists() and candidate.is_file():
                self.schema.update(self._parse_cpp_parameters(candidate))

        pele_like = self.root / "Source" / "Params" / "_cpp_parameters"
        if pele_like.exists() and pele_like.is_file():
            self.schema.update(self._parse_cpp_parameters(pele_like))

        # Phase 1: declaration map.
        files = self._iter_source_files(source_dirs)
        self.declaration_map = DeclarationExtractor.extract_all(files)

        # Phase 2: ParmParse scan.
        for f in files:
            self._scan_file(f)

        # Apply config build requirements.
        for param_name, flags in config.build_requirements.items():
            if param_name in self.schema and not self.schema[param_name].get("build_flags"):
                self.schema[param_name]["build_flags"] = list(flags)

        # Apply manual schema params.
        for param_name, info in config.manual_schema_params.items():
            if param_name in self.schema:
                continue
            self.schema[param_name] = {
                "type": info.get("type", "string"),
                "required": bool(info.get("required", False)),
                "is_array": bool(info.get("is_array", False)),
                "build_flags": list(info.get("build_flags", [])),
                "default": info.get("default"),
                "source_file": info.get("source_file", "manual_config"),
                "source_line": info.get("source_line"),
                "source_type": "manual",
                "priority": param_priority(param_name),
                "description": info.get("description"),
                "method": info.get("method", "manual"),
            }

        self._emit_tier2_warning_uq_recommendation()
        self._log_non_blocking_tier34_issues()
        return self.schema

    def _parse_cpp_parameters(self, file_path: Path) -> dict[str, dict[str, Any]]:
        out: dict[str, dict[str, Any]] = {}
        text = file_path.read_text(encoding="utf-8", errors="ignore")
        m = re.search(r"@namespace:\s+(\w+)", text)
        if not m:
            return out
        namespace = m.group(1)

        for line in text.splitlines():
            s = line.strip()
            if not s or s.startswith("#") or s.startswith("@"):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue

            param = parts[0]
            ptype = parts[1]
            default_s = " ".join(parts[2:]) if len(parts) > 2 else None

            if param.startswith("(") and ")" in param:
                param = param.strip("()").split(",")[0].strip()

            full = f"{namespace}.{param}"
            mapped_type = {
                "int": "int",
                "Real": "Real",
                "bool": "bool",
                "string": "string",
                "dim_array": "RealVect",
            }.get(ptype, ptype)

            default_v = default_s
            if default_s:
                if ptype == "Real":
                    with contextlib.suppress(ValueError):
                        default_v = float(default_s)
                elif ptype == "int":
                    with contextlib.suppress(ValueError):
                        default_v = int(default_s)
                elif ptype == "bool":
                    default_v = default_s.lower() in {"true", "1"}

            out[full] = {
                "type": mapped_type,
                "required": False,
                "is_array": ptype in {"dim_array"},
                "build_flags": [],
                "default": default_v,
                "description": None,
                "source_file": str(file_path.relative_to(self.root)),
                "source_line": None,
                "source_type": "generated",
                "priority": param_priority(full),
                "method": "_cpp_parameters",
            }

        return out

    def _scan_file(self, file_path: Path) -> None:
        try:
            text = file_path.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            return

        lines = text.splitlines()
        ifdef_blocks = parse_ifdefs(lines)

        # Map string variables for namespace resolution.
        string_vars: dict[str, str] = {
            # Common REMORA member prefix variables.
            "pp_prefix": "remora",
        }

        string_assign = re.compile(
            r"(?:const\s+)?(?:std::string|string)\s+([A-Za-z_]\w*)\s*(?:=|\{)\s*([^;]+?)(?:;|\})"
        )
        for line in lines:
            m = string_assign.search(line)
            if not m:
                continue
            var_name = m.group(1)
            expr = m.group(2).strip()
            resolved = resolve_string_expr(expr, string_vars)
            if resolved is not None:
                string_vars[var_name] = resolved

        # Build ParmParse var -> namespace history map.
        pp_history: dict[str, list[tuple[int, str]]] = defaultdict(list)
        pp_decl = re.compile(r"(?:amrex::)?ParmParse\s+([A-Za-z_]\w*)\s*(?:\(\s*([^)]*?)\s*\))?\s*;")

        for idx, line in enumerate(lines):
            for m in pp_decl.finditer(line):
                pp_var = m.group(1)
                arg = m.group(2)
                if arg is None or not arg.strip():
                    pp_history[pp_var].append((idx, ""))
                    continue

                arg = arg.strip()
                resolved = resolve_string_expr(arg, string_vars)

                # If not resolvable but looks like bare var, try string var map.
                if resolved is None and re.fullmatch(r"[A-Za-z_]\w*", arg):
                    resolved = string_vars.get(arg)

                if resolved:
                    pp_history[pp_var].append((idx, resolved))

            # Track inline assignments that might appear after declarations.
            m = string_assign.search(line)
            if m:
                var_name = m.group(1)
                expr = m.group(2).strip()
                resolved = resolve_string_expr(expr, string_vars)
                if resolved is not None:
                    string_vars[var_name] = resolved

        def namespace_for(pp_var: str, line_num: int) -> str:
            history = pp_history.get(pp_var)
            if not history:
                return ""
            for decl_line, prefix in reversed(history):
                if decl_line <= line_num:
                    return prefix
            return ""

        method_group = "|".join(re.escape(m) for m in METHODS)
        call_start = re.compile(rf"([A-Za-z_]\w*)\.({method_group})\s*\(")

        i = 0
        while i < len(lines):
            line = lines[i]
            start = call_start.search(line)
            if not start:
                i += 1
                continue

            start_line = i
            call_buf = line[start.start():]
            while ");" not in call_buf and i + 1 < len(lines):
                i += 1
                call_buf += "\n" + lines[i]

            m = re.search(rf"([A-Za-z_]\w*)\.({method_group})\s*\((.*)\)\s*;", call_buf, re.DOTALL)
            if not m:
                i += 1
                continue

            pp_var = m.group(1)
            method = m.group(2)
            args = m.group(3).strip()
            parts = split_top_level(args)
            if len(parts) < 2:
                i += 1
                continue

            p1 = parts[0].strip()
            if not re.fullmatch(r'"[^"]+"', p1):
                i += 1
                continue

            name1 = p1.strip('"')
            part_idx = 1
            alias_name: str | None = None
            if part_idx < len(parts) and re.fullmatch(r'"[^"]+"', parts[part_idx].strip()):
                alias_name = parts[part_idx].strip().strip('"')
                part_idx += 1

            if part_idx >= len(parts):
                i += 1
                continue

            var_expr = parts[part_idx]
            token_match = re.search(r"([A-Za-z_]\w*)", var_expr)
            var_name = token_match.group(1) if token_match else ""

            param_name = alias_name if alias_name else name1

            prefix = namespace_for(pp_var, start_line)
            if "." in param_name:
                full_name = param_name
            elif prefix:
                full_name = f"{prefix}.{param_name}"
            else:
                full_name = param_name

            # Keep dynamic prefixes deterministic but valid.
            full_name = full_name.replace("<", "").replace(">", "")
            full_name = full_name.replace("[", "_").replace("]", "")

            if not is_valid_param_name(full_name):
                i += 1
                continue

            var_type = DeclarationExtractor.match(var_name, file_path, self.declaration_map)
            if not var_type:
                var_type = detect_type_inline(var_name, lines, start_line, text)

            line_num = start_line
            source_type = source_type_for_file(file_path, self.root)

            entry = {
                "type": normalize_type(var_type),
                "required": method in REQUIRED_METHODS,
                "is_array": method in ARRAY_METHODS,
                "build_flags": ifdef_blocks.get(line_num, []),
                "default": extract_default(var_name, lines, start_line),
                "source_file": str(file_path.relative_to(self.root)),
                "source_line": line_num + 1,
                "source_type": source_type,
                "priority": param_priority(full_name),
                "description": extract_description(lines, start_line)
                or extract_var_description(var_name, lines, start_line),
                "method": method,
            }

            self._merge_entry(full_name, entry)
            i += 1

        # parseUserKey helper handling.
        helper_re = re.compile(r'parseUserKey\(\s*(\w+)\s*,\s*"([^"]+)"\s*,\s*\w+\s*,\s*(\w+)')
        for m in helper_re.finditer(text):
            pp_var = m.group(1)
            param = m.group(2)
            line_num = text[:m.start()].count("\n")
            prefix = namespace_for(pp_var, line_num)
            full_name = f"{prefix}.{param}" if prefix else param
            if not is_valid_param_name(full_name):
                continue
            entry = {
                "type": "string",
                "required": False,
                "is_array": True,
                "build_flags": ifdef_blocks.get(line_num, []),
                "default": None,
                "source_file": str(file_path.relative_to(self.root)),
                "source_line": line_num + 1,
                "source_type": source_type_for_file(file_path, self.root),
                "priority": param_priority(full_name),
                "description": extract_description(lines, line_num),
                "method": "parseUserKey",
            }
            self._merge_entry(full_name, entry)

    def _merge_entry(self, name: str, entry: dict[str, Any]) -> None:
        if name not in self.schema:
            self.schema[name] = entry
            return

        cur = self.schema[name]

        # Prefer generated over manual fallback entries.
        if cur.get("source_type") == "manual" and entry.get("source_type") == "generated":
            self.schema[name] = entry
            return

        cur["required"] = bool(cur.get("required")) or bool(entry.get("required"))
        cur["is_array"] = bool(cur.get("is_array")) or bool(entry.get("is_array"))

        if not cur.get("type") and entry.get("type"):
            cur["type"] = entry["type"]
        if cur.get("default") is None and entry.get("default") is not None:
            cur["default"] = entry["default"]
        if not cur.get("description") and entry.get("description"):
            cur["description"] = entry["description"]

        if not cur.get("build_flags") and entry.get("build_flags"):
            cur["build_flags"] = entry["build_flags"]

    def _emit_tier2_warning_uq_recommendation(self) -> None:
        missing = [p for p in TIER2_STABILITY_PARAMETERS if p not in self.schema]
        if missing:
            print(
                f"warning: {UNNUMBERED_112_ID}: missing tier2 stability parameters: {missing}. "
                "recommendation: include these in UQ sweep planning."
            )

    def _log_non_blocking_tier34_issues(self) -> None:
        tier3 = 0
        tier4 = 0
        for data in self.schema.values():
            p = str(data.get("priority", "tier4")).strip().lower()
            if p == "tier3":
                tier3 += 1
            elif p == "tier4":
                tier4 += 1
            elif p not in {"tier1", "tier2"}:
                tier4 += 1

        if tier3:
            print(f"info: tier3 non-blocking parameters: {tier3}")
        if tier4:
            print(f"info: tier4 non-blocking parameters: {tier4}")


def discover_dependencies(repo_root: Path) -> list[Dependency]:
    gitmodules = repo_root / ".gitmodules"
    if not gitmodules.exists():
        return []

    deps: list[Dependency] = []
    current_name: str | None = None
    current_path: str | None = None

    for line in gitmodules.read_text(encoding="utf-8", errors="ignore").splitlines():
        s = line.strip()
        m = re.match(r'\[submodule\s+"([^"]+)"\]', s)
        if m:
            if current_name and current_path:
                dep_path = (repo_root / current_path).resolve()
                deps.append(
                    Dependency(
                        name=Path(current_name).name,
                        path=dep_path,
                        commit=full_commit(dep_path),
                    )
                )
            current_name = m.group(1)
            current_path = None
            continue

        m = re.match(r"path\s*=\s*(.+)", s)
        if m:
            current_path = m.group(1).strip()

    if current_name and current_path:
        dep_path = (repo_root / current_path).resolve()
        deps.append(
            Dependency(
                name=Path(current_name).name,
                path=dep_path,
                commit=full_commit(dep_path),
            )
        )

    return deps


def resolve_solver_config(repo_path: Path, include_amrex_for_main: bool = True) -> SolverConfig:
    name = repo_path.name.lower()

    # AMReX submodule schema sources.
    if name == "amrex":
        return SolverConfig(
            schema_source_patterns=["Src/Base", "Src/Amr", "Src/AmrCore", "Src/Boundary"],
            build_requirements={},
            manual_schema_params={},
            schema_version=1,
        )

    # REMORA defaults.
    source_patterns = ["Source"]
    if include_amrex_for_main:
        source_patterns.extend(
            [
                "Submodules/AMReX/Src/Base",
                "Submodules/AMReX/Src/Amr",
                "Submodules/AMReX/Src/AmrCore",
                "Submodules/AMReX/Src/Boundary",
            ]
        )

    build_requirements = {
        "amr.max_level": [],
        "amr.blocking_factor": [],
    }

    manual_schema_params = {
        # Keep compatibility for commonly expected legacy keys when parser misses them.
        "max_step": {
            "type": "int",
            "required": False,
            "default": None,
            "description": "Legacy top-level AMReX step limit (manual fallback).",
        },
        "stop_time": {
            "type": "Real",
            "required": False,
            "default": None,
            "description": "Legacy top-level AMReX stop time (manual fallback).",
        },
    }

    return SolverConfig(
        schema_source_patterns=source_patterns,
        build_requirements=build_requirements,
        manual_schema_params=manual_schema_params,
        schema_version=1,
    )


def compose_schemas(component_schemas: list[tuple[str, dict[str, Any]]]) -> OrderedDict[str, dict[str, Any]]:
    merged: dict[str, dict[str, Any]] = {}
    for _, schema in component_schemas:
        for k, v in schema.items():
            merged[k] = v

    out = OrderedDict()
    for k in sorted(merged):
        out[k] = merged[k]
    return out


def save_schema_payload(payload: dict[str, Any], output_dir: Path, prefix: str) -> tuple[Path, Path]:
    generated = output_dir / ".generated"
    generated.mkdir(parents=True, exist_ok=True)

    repo_commit_short = payload.get("metadata", {}).get("repo_commit", "unknown")
    versioned = generated / f"{prefix}_{repo_commit_short}.json"
    current = generated / "remora_schema_current.json"

    versioned.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    if current.exists() or current.is_symlink():
        current.unlink()
    current.symlink_to(versioned.name)

    return versioned, current


def build_component_schema(component_root: Path, component_name: str, include_amrex_for_main: bool = False) -> dict[str, Any]:
    cfg = resolve_solver_config(component_root, include_amrex_for_main=include_amrex_for_main)
    builder = SchemaBuilder(component_root, component_name)
    schema = builder.scan_source_code(cfg.schema_source_patterns, cfg)
    return schema


def build_single(repo_root: Path, output_dir: Path, include_amrex: bool) -> tuple[Path, Path]:
    cfg = resolve_solver_config(repo_root, include_amrex_for_main=include_amrex)
    builder = SchemaBuilder(repo_root, "REMORA")
    schema = builder.scan_source_code(cfg.schema_source_patterns, cfg)

    payload = {
        "metadata": {
            "solver": "REMORA",
            "schema_version": cfg.schema_version,
            "repo_commit": short_commit(repo_root),
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "composition_method": "single",
            "total_parameters": len(schema),
            "source_roots": cfg.schema_source_patterns,
        },
        "parameters": OrderedDict((k, schema[k]) for k in sorted(schema)),
    }

    return save_schema_payload(payload, output_dir, "remora_schema")


def build_with_auto_compose(repo_root: Path, output_dir: Path) -> tuple[Path, Path]:
    deps = discover_dependencies(repo_root)

    component_schemas: list[tuple[str, dict[str, Any]]] = []
    component_commits: dict[str, str] = {}

    for dep in deps:
        if not dep.path.exists():
            print(f"warning: skipping dependency {dep.name} (missing path {dep.path})")
            continue
        schema = build_component_schema(dep.path, dep.name)
        component_schemas.append((dep.name, schema))
        component_commits[dep.name.lower()] = dep.commit
        print(f"built dependency schema: {dep.name} ({len(schema)} parameters)")

    remora_schema = build_component_schema(repo_root, "REMORA", include_amrex_for_main=False)
    component_schemas.append(("REMORA", remora_schema))
    component_commits["remora"] = full_commit(repo_root)

    composed = compose_schemas(component_schemas)
    short_hash = short_commit(repo_root)

    payload = {
        "metadata": {
            "solver": "REMORA",
            "schema_version": 1,
            "repo_commit": short_hash,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "composition_method": "auto-compose (.gitmodules)",
            "composed_from": [name for name, _ in component_schemas],
            "repo_commits": component_commits,
            "total_parameters": len(composed),
        },
        "parameters": composed,
    }

    return save_schema_payload(payload, output_dir, "remora_complete_schema")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build REMORA schema from ParmParse usage.")
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="REMORA repository root",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory for generated schema artifacts",
    )
    parser.add_argument(
        "--auto-compose",
        dest="auto_compose",
        action="store_true",
        default=True,
        help="Discover dependencies from .gitmodules and compose full schema (default)",
    )
    parser.add_argument(
        "--no-compose",
        dest="auto_compose",
        action="store_false",
        help="Disable dependency composition and build a single schema",
    )
    parser.add_argument(
        "--no-amrex",
        action="store_true",
        help="When --no-compose is set, skip AMReX submodule scan",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    output_dir = args.output_dir.resolve()

    if not args.auto_compose:
        versioned, current = build_single(repo_root, output_dir, include_amrex=not args.no_amrex)
    else:
        versioned, current = build_with_auto_compose(repo_root, output_dir)

    print(f"wrote schema: {versioned}")
    print(f"updated symlink: {current}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
