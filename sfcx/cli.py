# -*- coding: utf-8 -*-
"""Command line interface."""

from __future__ import annotations

import argparse
import json
import os
import sys

from .coords import Coordinate, METRICS, describe
from .data import Cohort
from . import pipeline, demo, report


def _load(args):
    if args.demo:
        c = demo.make_cohort(n_subjects=args.subjects, n_nodes=args.nodes, seed=args.seed)
        return c
    if not args.npz:
        raise SystemExit("give --npz FILE or --demo")
    return Cohort.from_npz(args.npz, state=args.state)


def cmd_describe(args):
    print(describe())


def cmd_demo(args):
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    coh = demo.make_cohort(n_subjects=args.subjects, n_nodes=args.nodes, seed=args.seed)
    demo.write_npz(os.path.join(out, "demo_cohort.npz"), coh)
    print("[demo] cohort %r" % (coh,))

    res = pipeline.run(coh, Coordinate.from_preset(args.preset), outdir=os.path.join(out, "route"),
                       n_perm=args.perm, seed=args.seed)
    print("[demo] route axis r=%.3f p=%.4f null=%s"
          % (res.get("axis", {}).get("r", float("nan")),
             res.get("axis", {}).get("p", float("nan")),
             res.get("null_model")))

    sw = pipeline.sweep(coh, outdir=os.path.join(out, "sweep"), n_perm=max(200, args.perm // 2),
                        seed=args.seed)
    print("[demo] sweep r range %.2f .. %.2f over %d definitions x %d estimators"
          % (float(sw["r"].min()), float(sw["r"].max()), len(sw["rows"]), len(sw["cols"])))

    notes = ("<b>Synthetic cohort.</b> Structure is modular, function is an individually weighted "
             "noisy copy of structure, and the noise load rises towards association cortex. "
             "The spread of axis values below is an outcome of those assumptions, not a scripted result.")
    html_path = report.build(res, os.path.join(out, "report.html"),
                             title="SFCX demo - one dataset, many coordinates",
                             sweep=sw, notes=notes)
    print("[demo] report -> %s" % html_path)
    return 0


def cmd_run(args):
    coh = _load(args)
    coord = Coordinate.from_preset(args.preset) if args.preset else Coordinate(
        metric=args.metric, pipeline=args.pipeline, estimator=args.estimator, scale=args.scale)
    res = pipeline.run(coh, coord, outdir=os.path.abspath(args.out), n_perm=args.perm,
                       null=args.null, seed=args.seed)
    report.build(res, os.path.join(os.path.abspath(args.out), "report.html"),
                 title="SFCX route report")
    print(json.dumps({k: res[k] for k in ("coordinate", "axis", "layers", "reliability")
                      if k in res}, indent=2)[:1800])
    return 0


def cmd_sweep(args):
    coh = _load(args)
    sw = pipeline.sweep(coh, outdir=os.path.abspath(args.out), n_perm=args.perm,
                        null=args.null, seed=args.seed)
    print("[sweep] r matrix")
    print("        " + "  ".join("%-9s" % c for c in sw["cols"]))
    for i, m in enumerate(sw["rows"]):
        print("%-9s" % m + "  ".join("%-9.3f" % sw["r"][i, j] for j in range(len(sw["cols"]))))
    print("[sweep] -> %s" % os.path.abspath(args.out))
    return 0


def build_parser():
    p = argparse.ArgumentParser(prog="sfcx", description="Structure-function coupling explorer")
    sub = p.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("describe", help="print the coordinate system")
    d.set_defaults(func=cmd_describe)

    d = sub.add_parser("demo", help="generate a synthetic cohort and run the full route")
    d.add_argument("--out", default="sfcx_demo")
    d.add_argument("--subjects", type=int, default=320)
    d.add_argument("--nodes", type=int, default=180)
    d.add_argument("--perm", type=int, default=500)
    d.add_argument("--seed", type=int, default=0)
    d.add_argument("--preset", default="default")
    d.set_defaults(func=cmd_demo)

    r = sub.add_parser("run", help="run one route on a cohort")
    r.add_argument("--npz")
    r.add_argument("--demo", action="store_true")
    r.add_argument("--out", default="sfcx_out")
    r.add_argument("--preset")
    r.add_argument("--metric", default="communicability")
    r.add_argument("--pipeline", default="pearson")
    r.add_argument("--estimator", default="regional")
    r.add_argument("--scale", default="regional")
    r.add_argument("--state", default="rest")
    r.add_argument("--null", default="auto", choices=["auto", "spin", "msr", "permutation"])
    r.add_argument("--perm", type=int, default=500)
    r.add_argument("--subjects", type=int, default=320)
    r.add_argument("--nodes", type=int, default=180)
    r.add_argument("--seed", type=int, default=0)
    r.set_defaults(func=cmd_run)

    s = sub.add_parser("sweep", help="axis strength across definitions x estimators")
    s.add_argument("--npz")
    s.add_argument("--demo", action="store_true")
    s.add_argument("--out", default="sfcx_sweep")
    s.add_argument("--state", default="rest")
    s.add_argument("--null", default="auto", choices=["auto", "spin", "msr", "permutation"])
    s.add_argument("--perm", type=int, default=300)
    s.add_argument("--subjects", type=int, default=320)
    s.add_argument("--nodes", type=int, default=180)
    s.add_argument("--seed", type=int, default=0)
    s.set_defaults(func=cmd_sweep)

    w = sub.add_parser("web", help="open the browser front-end")
    w.add_argument("--host", default="127.0.0.1")
    w.add_argument("--port", type=int, default=8765)
    w.add_argument("--workdir", default=None, help="where uploads and results are kept")
    w.add_argument("--open", action="store_true", help="open the page in the default browser")
    w.set_defaults(func=cmd_web)
    return p


def cmd_web(args):
    from .webapp import serve
    serve(host=args.host, port=args.port, root=args.workdir, open_browser=args.open)
    return 0


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
