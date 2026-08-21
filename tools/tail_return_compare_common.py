#!/usr/bin/env python3
"""Shared readers for H10 tail-return acceptance tools."""

from __future__ import print_function

import glob
import math
import os


def read_table(path):
    with open(path, "r") as stream:
        lines = [line.split() for line in stream if line.strip() and
                 not line.startswith("#")]
    if not lines:
        raise ValueError("empty table: " + path)
    header = lines[0]
    rows = []
    for fields in lines[1:]:
        if len(fields) != len(header):
            continue
        try:
            rows.append(dict(zip(header, [float(x) for x in fields])))
        except ValueError:
            continue
    return header, rows


def diagnostics(root):
    path = os.path.join(root, "vpfp_step_diagnostics.dat")
    if not os.path.isfile(path):
        raise ValueError("missing diagnostics: " + path)
    return read_table(path)[1]


def latest_snapshot(root):
    candidates = [p for p in glob.glob(os.path.join(root, "snapshot_*"))
                  if os.path.isdir(p)]
    if not candidates:
        raise ValueError("no snapshot directory under " + root)
    def snapshot_time(path):
        name = os.path.basename(path)
        try:
            return float(name.split("_t", 1)[1].split("fs", 1)[0])
        except (IndexError, ValueError):
            return -math.inf
    return max(candidates, key=snapshot_time)


def rank_series(snapshot, stem, column):
    points = []
    for path in glob.glob(os.path.join(snapshot, stem + "_rank*.dat")):
        header, rows = read_table(path)
        if column not in header:
            raise ValueError("missing column %s in %s" % (column, path))
        for row in rows:
            points.append((row[header[0]], row[column]))
    if not points:
        raise ValueError("no %s rank files under %s" % (stem, snapshot))
    points.sort()
    return [p[1] for p in points]


def rank_sum_series(snapshot, stem, column):
    totals = {}
    for path in glob.glob(os.path.join(snapshot, stem + "_rank*.dat")):
        header, rows = read_table(path)
        if column not in header:
            raise ValueError("missing column %s in %s" % (column, path))
        for row in rows:
            key = row[header[0]]
            totals[key] = totals.get(key, 0.0) + row[column]
    if not totals:
        raise ValueError("no %s rank files under %s" % (stem, snapshot))
    return [totals[key] for key in sorted(totals)]


def relative_l2(a, b):
    if len(a) != len(b) or not a:
        return math.inf
    numerator = math.sqrt(sum((x-y)*(x-y) for x, y in zip(a, b)))
    denominator = max(math.sqrt(sum(x*x for x in a)), 1.0e-300)
    return numerator / denominator


def max_value(rows, names, default=0.0):
    values = [abs(row.get(name, default)) for row in rows for name in names]
    return max(values) if values else default


def last_value(rows, name, default=math.nan):
    return rows[-1].get(name, default) if rows else default


def failed(root):
    path = os.path.join(root, "vpfp_failure.dat")
    return os.path.isfile(path) and os.path.getsize(path) > 0
