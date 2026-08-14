#!/usr/bin/env python3

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import build_schema_remora as schema_mod


def write_file(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


class SchemaBuilderTests(unittest.TestCase):
    def build_schema(self, cpp_source: str) -> dict[str, dict]:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            write_file(root / "Source" / "sample.cpp", cpp_source)
            cfg = schema_mod.SolverConfig(
                schema_source_patterns=["Source"],
                build_requirements={},
                manual_schema_params={},
                schema_version=1,
            )
            builder = schema_mod.SchemaBuilder(root, "REMORA")
            return builder.scan_source_code(["Source"], cfg)

    def test_namespace_history_prefers_latest_declaration(self) -> None:
        src = """
        #include <AMReX_ParmParse.H>
        void f () {
            int max_step = 10;
            amrex::Real cfl = 0.9;
            std::string pp_prefix {"remora"};
            amrex::ParmParse pp;
            pp.queryAdd("max_step", max_step);
            amrex::ParmParse pp(pp_prefix);
            pp.queryAdd("cfl", cfl);
        }
        """
        schema = self.build_schema(src)
        self.assertIn("max_step", schema)
        self.assertIn("remora.cfl", schema)
        self.assertNotIn("cfl", schema)

    def test_method_semantics_required_and_array(self) -> None:
        src = """
        #include <AMReX_ParmParse.H>
        void f () {
            int i = 0;
            std::vector<int> arr;
            amrex::ParmParse pp("amr");
            pp.get("a", i);
            pp.query("b", i);
            pp.getarr("c", arr);
            pp.queryarr("d", arr);
            pp.getWithParser("e", i);
            pp.querytable("f", arr);
        }
        """
        schema = self.build_schema(src)

        self.assertTrue(schema["amr.a"]["required"])
        self.assertFalse(schema["amr.b"]["required"])
        self.assertTrue(schema["amr.c"]["required"])
        self.assertTrue(schema["amr.e"]["required"])

        self.assertTrue(schema["amr.c"]["is_array"])
        self.assertTrue(schema["amr.d"]["is_array"])
        self.assertTrue(schema["amr.f"]["is_array"])

    def test_type_default_and_comment_extraction(self) -> None:
        src = """
        #include <AMReX_ParmParse.H>
        void f () {
            amrex::ParmParse pp("remora");
            amrex::Real cfl = 0.45; // target cfl value
            pp.queryAdd("cfl", cfl);
        }
        """
        schema = self.build_schema(src)
        entry = schema["remora.cfl"]

        self.assertEqual(entry["type"], "Real")
        self.assertEqual(entry["default"], 0.45)
        self.assertIn("target cfl value", entry["description"])

    def test_compose_prefers_later_component_override(self) -> None:
        composed = schema_mod.compose_schemas(
            [
                (
                    "AMReX",
                    {
                        "amr.max_level": {
                            "type": "int",
                            "source_file": "amrex.cpp",
                        }
                    },
                ),
                (
                    "REMORA",
                    {
                        "amr.max_level": {
                            "type": "int",
                            "source_file": "remora.cpp",
                        },
                        "remora.cfl": {
                            "type": "Real",
                            "source_file": "remora.cpp",
                        },
                    },
                ),
            ]
        )

        self.assertEqual(composed["amr.max_level"]["source_file"], "remora.cpp")
        self.assertIn("remora.cfl", composed)

    def test_discover_dependencies_from_gitmodules(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            write_file(
                root / ".gitmodules",
                """
                [submodule \"Submodules/AMReX\"]
                path = Submodules/AMReX
                url = https://github.com/AMReX-Codes/amrex.git
                """,
            )
            (root / "Submodules" / "AMReX").mkdir(parents=True, exist_ok=True)

            deps = schema_mod.discover_dependencies(root)
            self.assertEqual(len(deps), 1)
            self.assertEqual(deps[0].name, "AMReX")
            self.assertEqual(deps[0].path, (root / "Submodules" / "AMReX").resolve())


if __name__ == "__main__":
    unittest.main()
