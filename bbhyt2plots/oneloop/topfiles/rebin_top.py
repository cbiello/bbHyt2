#!/usr/bin/env python3
"""Merge adjacent bins of a POWHEG .top file (density histograms).

    python3 rebin_top.py in.top out.top  name=factor [name=factor ...]

value_new = sum(v_i w_i)/W ,  err_new = sqrt(sum (e_i w_i)^2)/W ,  W = sum w_i
Histograms not listed are copied unchanged; a trailing incomplete group
is merged as it is.
"""
import sys, math

fin, fout = sys.argv[1], sys.argv[2]
fac = dict((k, int(v)) for k, v in (a.split("=") for a in sys.argv[3:]))

blocks, cur = [], None
for ln in open(fin):
    s = ln.strip()
    if s.startswith("#"):
        cur = [ln.rstrip("\n"), s[1:].split("index")[0].strip(), []]
        blocks.append(cur); continue
    if s and cur is not None:
        cur[2].append([float(x.replace("D", "E")) for x in s.split()])

with open(fout, "w") as f:
    for head, name, rows in blocks:
        f.write(head + "\n")
        k = fac.get(name, 1)
        for i in range(0, len(rows), k):
            g = rows[i:i + k]
            W = sum(r[1] - r[0] for r in g)
            v = sum(r[2] * (r[1] - r[0]) for r in g) / W
            e = math.sqrt(sum((r[3] * (r[1] - r[0])) ** 2 for r in g)) / W
            f.write(" %.8E %.8E %.8E %.8E\n" % (g[0][0], g[-1][1], v, e))
        f.write("\n\n")
