#!/usr/bin/env python3
"""Tabulate the runs made by perlmutter_phase1.sh.

    python3 Tests/perf/summarize.py <workdir>/runs

Each run directory is <arch>_<case>_<config>_<variant>/log with a TinyProfiler table; the
Evolve inclusive time is the metric. The fcompare logs give the bitwise verdict.
"""
import glob
import os
import re
import sys


def evolve_time(log):
    """Evolve time from the inclusive TinyProfiler table; failing that, the run's total."""
    t = None
    total = None
    incl = False
    try:
        with open(log) as f:
            for line in f:
                if "Incl. Min" in line:          # TinyProfiler prints the exclusive table first
                    incl = True
                elif incl and line.startswith("REMORA::Evolve()"):
                    t = float(line.split()[-2])
                elif "TinyProfiler total time" in line:
                    total = float(line.split()[-1])
    except (OSError, ValueError):
        pass
    if t is None and total is not None:
        print(f"note: {log}: no inclusive Evolve line, using the TinyProfiler total")
        t = total
    return t


def main(runs):
    rows = {}
    for d in sorted(glob.glob(os.path.join(runs, "*_*_*_*"))):
        m = re.match(r"(cpu|gpu)_(\w+?)_(\w+?)_(base|base_j16|new)$", os.path.basename(d))
        if not m:
            continue
        arch, case, cfg, var = m.groups()
        rows.setdefault((arch, case, cfg), {})[var] = evolve_time(os.path.join(d, "log"))
    print(f"{'arch':4} {'case':10} {'config':7} {'base':>9} {'base_j16':>9} {'new':>9} "
          f"{'new/base':>9}  bitwise")
    for (arch, case, cfg), v in rows.items():
        fcl = os.path.join(runs, f"{arch}_{case}_{cfg}_fcompare.log")
        verdict = "-"
        if os.path.exists(fcl):
            txt = open(fcl).read()
            verdict = "IDENTICAL" if "AGREE" in txt else "DIFFER"
        def fmt(x):
            return f"{x:9.2f}" if x is not None else f"{'--':>9}"
        ratio = (f"{v['new']/v['base']:9.2f}"
                 if v.get("new") and v.get("base") else f"{'--':>9}")
        print(f"{arch:4} {case:10} {cfg:7} {fmt(v.get('base'))} {fmt(v.get('base_j16'))} "
              f"{fmt(v.get('new'))} {ratio}  {verdict}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
