#!/usr/bin/env python3
"""Tabulate the runs made by perlmutter_phase1.sh.

    python3 Tests/perf/summarize.py <workdir>/runs

Each run directory is <arch>_<case>_<config>_<variant>/log with a TinyProfiler table; the
Evolve inclusive time is the metric. The variants are the revision labels given to the
driver (the first is the baseline) plus <baseline>_j16 on CPU. The fcompare logs
(<arch>_<case>_<config>_fcompare_<variant>.log) give the bitwise verdict against the
baseline.
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
                elif incl and t is None and line.startswith("REMORA::Evolve()"):
                    t = float(line.split()[-2])  # first row only: GPU builds add a memory table after
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
    variants = []
    for d in sorted(glob.glob(os.path.join(runs, "*_*_*_*"))):
        m = re.match(r"(cpu|gpu)_([A-Za-z0-9]+)_([A-Za-z0-9]+)_(\w+)$", os.path.basename(d))
        if not m or not os.path.exists(os.path.join(d, "log")):
            continue
        arch, case, cfg, var = m.groups()
        rows.setdefault((arch, case, cfg), {})[var] = evolve_time(os.path.join(d, "log"))
        if var not in variants:
            variants.append(var)
    if not rows:
        print("no runs found under", runs)
        return
    # the baseline is the variant that has a fcompare log for every other one
    base = next((v for v in variants
                 if glob.glob(os.path.join(runs, f"*_fcompare_*.log")) and
                 not glob.glob(os.path.join(runs, f"*_fcompare_{v}.log"))), variants[0])
    others = [v for v in variants if v != base]

    def fmt(x):
        return f"{x:9.2f}" if x is not None else f"{'--':>9}"

    def verdict(fcl):
        # IDENTICAL, the largest relative difference fcompare reported, or DIFFER when the
        # log has no table (a run that did not finish).
        if not os.path.exists(fcl):
            return "-"
        text = open(fcl).read()
        if "AGREE" in text:
            return "IDENTICAL"
        worst = None
        for ln in text.splitlines():
            parts = ln.split()
            if len(parts) == 3 and parts[0] != "level":
                try:
                    rel = float(parts[2])
                except ValueError:
                    continue
                worst = rel if worst is None else max(worst, rel)
        return f"{worst:.1e}" if worst is not None else "DIFFER"

    head = f"{'arch':4} {'case':10} {'config':7} {base:>9}"
    for v in others:
        head += f" {v:>9} {'/'+base:>7} {'bitwise':>9}"
    print(head)
    for (arch, case, cfg), t in rows.items():
        line = f"{arch:4} {case:10} {cfg:7} {fmt(t.get(base))}"
        for v in others:
            ratio = (f"{t[v]/t[base]:7.2f}" if t.get(v) and t.get(base) else f"{'--':>7}")
            fcl = os.path.join(runs, f"{arch}_{case}_{cfg}_fcompare_{v}.log")
            line += f" {fmt(t.get(v))} {ratio} {verdict(fcl):>9}"
        print(line)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
