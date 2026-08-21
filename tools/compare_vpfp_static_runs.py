#!/usr/bin/env python3
"""Compare production checkpoint bytes for the flux-audit reproducibility gate.

The first comparison is static-A versus static-B.  Those runs must be
bitwise identical before any difference involving flux-audit is interpreted.
The optional audit checkpoint comparison is deliberately informational: an
audit failure may still make the audit run fail its audit gate, but it must
not change the physical state bytes.
"""

from __future__ import print_function

import argparse
import hashlib
import os
import sys


def rank_hashes(path):
    if not os.path.isdir(path):
        raise RuntimeError("checkpoint directory does not exist: %s" % path)
    result = {}
    for name in sorted(os.listdir(path)):
        if not name.startswith("rank_") or not name.endswith(".bin"):
            continue
        full = os.path.join(path, name)
        if not os.path.isfile(full):
            continue
        digest = hashlib.sha256()
        with open(full, "rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        result[name] = digest.hexdigest()
    if not result:
        raise RuntimeError("no rank_*.bin files found in %s" % path)
    return result


def compare(label, left, right):
    if left == right:
        print("%s=1" % label)
        return True
    print("%s=0" % label)
    left_names = set(left)
    right_names = set(right)
    for name in sorted(left_names | right_names):
        if left.get(name) != right.get(name):
            print("mismatch=%s" % name)
    return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--static-a", required=True)
    parser.add_argument("--static-b", required=True)
    parser.add_argument("--audit")
    args = parser.parse_args()
    try:
        static_a = rank_hashes(args.static_a)
        static_b = rank_hashes(args.static_b)
        static_ok = compare("static_reproducible", static_a, static_b)
        audit_ok = True
        if args.audit:
            audit = rank_hashes(args.audit)
            audit_ok = compare("audit_physical_state_matches_static", static_a, audit)
        print("passes=%d" % (1 if static_ok and audit_ok else 0))
        return 0 if static_ok and audit_ok else 2
    except RuntimeError as error:
        print("error=%s" % error, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
