#!/usr/bin/env python3
"""Regression check: the Python bindings and the CLI must agree.

Because CGAL decimation and the flow/re-inflation are deterministic, the cages
returned by ``nested_cages.nested_cages(...)`` must match the ``output_*.obj``
files written by the ``nested_cages`` command-line tool for the same inputs
(up to the precision of the OBJ files the CLI writes).

Usage:
    python tests/regression.py --cli build/nested_cages --input gargo.off
"""
import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np


def read_off(path):
    with open(path) as f:
        toks = f.read().split()
    assert toks[0].upper().startswith("OFF")
    # header may be "OFF" alone or "OFF nv nf ne" fused; handle the common case
    idx = 1
    nv, nf = int(toks[idx]), int(toks[idx + 1])
    idx += 3  # skip nv nf ne
    V = np.array(toks[idx: idx + 3 * nv], dtype=np.float64).reshape(nv, 3)
    idx += 3 * nv
    F = []
    for _ in range(nf):
        c = int(toks[idx])
        F.append([int(x) for x in toks[idx + 1: idx + 1 + c]])
        idx += 1 + c
    return V, np.array(F, dtype=np.int32)


def read_obj(path):
    V, F = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("v "):
                V.append([float(x) for x in line.split()[1:4]])
            elif line.startswith("f "):
                F.append([int(t.split("/")[0]) - 1 for t in line.split()[1:4]])
    return np.array(V, dtype=np.float64), np.array(F, dtype=np.int32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cli", default="build/nested_cages")
    ap.add_argument("--input", default="gargo.off")
    ap.add_argument("--faces", default="1000,500")
    ap.add_argument("--quad", type=int, default=2)
    ap.add_argument("--regular", action="store_true", default=True)
    ap.add_argument("--energy-expansion", default="None")
    ap.add_argument("--energy-final", default="Volume")
    args = ap.parse_args()

    import nested_cages

    faces = [int(x) for x in args.faces.split(",")]
    V0, F0 = read_off(args.input)

    print(f"[py]  running nested_cages on {args.input}, faces={faces} ...")
    cages = nested_cages.nested_cages(
        V0, F0, faces,
        quad_order=args.quad,
        energy_expansion=args.energy_expansion,
        energy_final=args.energy_final,
        regular=[args.regular] * len(faces),
    )
    assert len(cages) == len(faces), f"expected {len(faces)} cages, got {len(cages)}"

    with tempfile.TemporaryDirectory() as td:
        prefix = os.path.join(td, "cli")
        face_args = [f"{n}r" if args.regular else str(n) for n in faces]
        cmd = [args.cli, args.input, str(args.quad), *face_args,
               args.energy_expansion, args.energy_final, prefix]
        print("[cli]", " ".join(cmd))
        # NOTE: the nested_cages CLI returns 1 on success and 0 on failure
        # (a long-standing quirk of the original program), so we do not treat
        # the exit code as pass/fail; we verify the output files instead.
        subprocess.run(cmd)

        ok = True
        for i, (Vp, Fp) in enumerate(cages):
            out = f"{prefix}_{i + 1}.obj"
            if not os.path.exists(out):
                print(f"  cage {i + 1}: CLI did not produce {out}  [MISMATCH]")
                ok = False
                continue
            Vc, Fc = read_obj(out)
            same_shape = Vp.shape == Vc.shape and Fp.shape == Fc.shape
            faces_eq = same_shape and np.array_equal(Fp, Fc)
            verts_eq = same_shape and np.allclose(Vp, Vc, rtol=1e-5, atol=1e-6)
            maxerr = np.abs(Vp - Vc).max() if same_shape else float("nan")
            status = "OK" if (faces_eq and verts_eq) else "MISMATCH"
            print(f"  cage {i + 1}: py {Vp.shape[0]}v/{Fp.shape[0]}f  "
                  f"cli {Vc.shape[0]}v/{Fc.shape[0]}f  max|dV|={maxerr:.2e}  [{status}]")
            ok = ok and faces_eq and verts_eq

    if not ok:
        print("REGRESSION FAILED: Python and CLI cages differ")
        sys.exit(1)
    print("REGRESSION PASSED: Python and CLI cages agree")


if __name__ == "__main__":
    main()
