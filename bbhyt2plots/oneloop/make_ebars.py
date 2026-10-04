#!/usr/bin/env python3
"""Bin-centre statistical error bars, to overlay on the histograms
(same format as the bbH GMVFNS plots, make_ebars.py there).

The eta band (columns 4/6 of the .dat) and the MC error (column 3) are two
different uncertainties: a shaded band for the eta variation, a vertical bar
on the central point for the statistics.

Reads the merged central .top files in topfiles/<run>/ and writes
    ebars/<obs>__<TAG>-<tag>.ebar    x_centre  value      stat_err      (fb, density)
    ebars/<obs>__<TAG>-<tag>.rebar   x_centre  value/REF  stat_err/REF  (ratio panel)
with REF the exact curve (the ratio panel is normalised to it).

    python3 make_ebars.py [topfiles/<TAG>-lhe-rebin] [TAG]
"""
import os, sys
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "ebars")
os.makedirs(OUT, exist_ok=True)
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "topfiles", "HT4-lhe-rebin")
TAG = sys.argv[2] if len(sys.argv) > 2 else "HT4"
TAGS = {"exact": "exact-central.top", "mass": "mass-central.top"}
REF = "exact"


def parse(p):
    b, c = OrderedDict(), None
    for ln in open(p):
        s = ln.strip()
        if s.startswith("#") and "index" in s:
            c = s[1:].split("index")[0].strip(); b.setdefault(c, []); continue
        if not s or s[0] not in "-.0123456789":
            continue
        t = s.replace("D", "E").split()
        if c is None or len(t) < 4:
            continue
        b[c].append([float(x) for x in t])
    return b


blocks = {k: parse(os.path.join(SRC, f)) for k, f in TAGS.items()}
R = blocks[REF]
n = 0
for tag, B in blocks.items():
    for name, rows in B.items():
        if name == "xsec":
            continue
        rows_main, rows_ratio = [], []
        for j, (lo, hi, v, e) in enumerate(rows):
            xc = 0.5 * (lo + hi)
            rows_main.append((xc, v * 1000.0, e * 1000.0))
            m = R[name][j][2]
            if m:
                rows_ratio.append((xc, v / m, e / m))
        for ext, rr in ((".ebar", rows_main), (".rebar", rows_ratio)):
            with open(os.path.join(OUT, f"{name}__{TAG}-{tag}{ext}"), "w") as f:
                for x, y, ee in rr:
                    f.write("%.8g %.8g %.8g\n" % (x, y, ee))
            n += 1
print(f"  wrote {n} files in {OUT}  (from {SRC})")
