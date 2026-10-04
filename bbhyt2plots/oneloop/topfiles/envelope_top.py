#!/usr/bin/env python3
"""Bin-by-bin envelope of several POWHEG .top files (same binning).

    python3 envelope_top.py out_min.top out_max.top in1.top in2.top ...

For every bin the variation with the smallest (largest) value is taken,
together with its own statistical error.
"""
import sys

def read(f):
    blocks, cur = [], None
    for ln in open(f):
        s = ln.strip()
        if s.startswith("#"):
            cur = [ln.rstrip("\n"), []]; blocks.append(cur); continue
        if s and cur is not None:
            cur[1].append([float(x.replace("D", "E")) for x in s.split()])
    return blocks

fmin, fmax, ins = sys.argv[1], sys.argv[2], [read(f) for f in sys.argv[3:]]
with open(fmin, "w") as gmin, open(fmax, "w") as gmax:
    for ib, (head, rows) in enumerate(ins[0]):
        for g in (gmin, gmax):
            g.write(head + "\n")
        for j, r in enumerate(rows):
            vals = [b[ib][1][j] for b in ins]
            lo = min(vals, key=lambda x: x[2]); hi = max(vals, key=lambda x: x[2])
            gmin.write(" %.8E %.8E %.8E %.8E\n" % tuple(lo))
            gmax.write(" %.8E %.8E %.8E %.8E\n" % tuple(hi))
        for g in (gmin, gmax):
            g.write("\n\n")
