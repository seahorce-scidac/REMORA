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
import render_inputs_generic as render_mod


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

    def test_comment_must_directly_precede_query(self) -> None:
        src = """
        #include <AMReX_ParmParse.H>
        void f () {
            amrex::Real a = 0.0;
            amrex::Real b = 0.0;
            amrex::Real c = 0.0;
            amrex::ParmParse pp("remora");
            // describes a
            pp.queryAdd("a", a);
            pp.queryAdd("b", b);
            // separated by a blank line

            pp.queryAdd("c", c);
        }
        """
        schema = self.build_schema(src)

        self.assertEqual(schema["remora.a"]["description"], "describes a")
        self.assertIsNone(schema["remora.b"]["description"])
        self.assertIsNone(schema["remora.c"]["description"])

    def test_declaration_comment_describes_query(self) -> None:
        src = """
        #include <AMReX_ParmParse.H>
        void f () {
            amrex::ParmParse pp("remora");
            // Flat-bottom depth. Set prob_lo's z to -depth to match.
            amrex::Real depth = 10.0;
            pp.query("depth", depth);
        }
        """
        schema = self.build_schema(src)

        self.assertEqual(
            schema["remora.depth"]["description"],
            "Flat-bottom depth. Set prob_lo's z to -depth to match.",
        )

    def test_only_nearest_declaration_describes_query(self) -> None:
        src = """
        #include <AMReX_ParmParse.H>
        void g () {
            // unrelated: depth of the sponge layer
            amrex::Real depth = 3.0;
        }
        void f () {
            amrex::Real depth = 10.0;
            amrex::ParmParse pp("remora");
            pp.query("depth", depth);
        }
        """
        schema = self.build_schema(src)
        self.assertIsNone(schema["remora.depth"]["description"])

    def test_trailing_block_comment_on_code_is_not_borrowed(self) -> None:
        src = """
        #include <AMReX_ParmParse.H>
        void f () {
            amrex::ParmParse pp("remora");
            int a = 0; /* a thing */
            int b = 0;
            pp.query("b", b);
        }
        """
        schema = self.build_schema(src)
        self.assertIsNone(schema["remora.b"]["description"])

    def test_long_and_block_comments_are_read_whole(self) -> None:
        src = """
        #include <AMReX_ParmParse.H>
        void f () {
            int a = 0;
            int b = 0;
            amrex::ParmParse pp("remora");
            // one
            // two
            // three
            // four
            // five
            // six
            // seven
            pp.query("a", a);
            /* first line
               second line */
            pp.query("b", b);
        }
        """
        schema = self.build_schema(src)
        self.assertEqual(
            schema["remora.a"]["description"],
            "one two three four five six seven",
        )
        self.assertEqual(
            schema["remora.b"]["description"], "first line second line"
        )

    def test_preprocessor_guard_keeps_description(self) -> None:
        src = """
        #include <AMReX_ParmParse.H>
        void f () {
            std::string x;
            amrex::ParmParse pp("remora");
            // NetCDF file name for x
        #ifdef REMORA_USE_NETCDF
            pp.query("nc_x", x);
        #endif
        }
        """
        schema = self.build_schema(src)
        self.assertEqual(
            schema["remora.nc_x"]["description"], "NetCDF file name for x"
        )

    def test_closing_preprocessor_line_separates_comment(self) -> None:
        src = """
        #include <AMReX_ParmParse.H>
        void f () {
            int a = 0;
            int b = 0;
            amrex::ParmParse pp("remora");
        #ifdef REMORA_USE_NETCDF
            pp.query("a", a);
            // TODO remove once X lands
        #endif
            pp.query("b", b);
        }
        """
        schema = self.build_schema(src)
        self.assertIsNone(schema["remora.b"]["description"])

    def test_unrecognized_type_declaration_is_nearest(self) -> None:
        src = """
        #include <AMReX_ParmParse.H>
        void g () {
            // unrelated string of the same name
            std::string ic_type;
        }
        void f () {
            IC_Type ic_type = IC_Type::analytic;
            amrex::ParmParse pp("remora");
            pp.query("ic_type", ic_type);
        }
        """
        schema = self.build_schema(src)
        self.assertIsNone(schema["remora.ic_type"]["description"])

    def test_declaration_of_other_types_describes_query(self) -> None:
        src = """
        #include <AMReX_ParmParse.H>
        void f () {
            // number of passive tracers
            amrex::Long ntr = 0;
            amrex::ParmParse pp("remora");
            pp.query("ntr", ntr);
        }
        """
        schema = self.build_schema(src)
        self.assertEqual(
            schema["remora.ntr"]["description"], "number of passive tracers"
        )

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


class DocsDescriptionTests(unittest.TestCase):
    def test_grid_table_cells_wrap_across_lines(self) -> None:
        table = """
+-----------------------+----------------------+-------------+
| Parameter             | Definition           | Values      |
+=======================+======================+=============+
| **remora.do_rivers**  | Whether rivers are   | true / false|
|                       |                      |             |
|                       | temperature sources  |             |
+-----------------------+----------------------+-------------+
| **remora.cfl**        | Courant number       | Real        |
+-----------------------+----------------------+-------------+
| **remora.empty**      |                      | Real        |
+-----------------------+----------------------+-------------+
"""
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "Inputs.rst"
            write_file(path, table)
            descs = render_mod.parse_docs_descriptions(path)

        self.assertEqual(
            descs["remora.do_rivers"], "Whether rivers are temperature sources"
        )
        self.assertEqual(descs["remora.cfl"], "Courant number")
        self.assertNotIn("remora.empty", descs)


if __name__ == "__main__":
    unittest.main()
