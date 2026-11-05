# -*- coding: utf-8 -*-
"""Zero-dependency web front-end.

    python -m sfcx.webapp --port 8765 --open

Serves a single-page interface and a small JSON API.  Uploads are sent as raw
bodies (no multipart parsing), analyses run in a worker thread and are polled
through ``/api/job``.
"""

from __future__ import annotations

import argparse
import json
import mimetypes
import os
import re
import shutil
import threading
import time
import traceback
import urllib.parse
import uuid
import webbrowser
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import numpy as np

from . import engine, io, report
from .coords import catalog
from .custom import CustomCodeError, dry_run, register_custom, unregister_custom, list_custom

HERE = os.path.dirname(os.path.abspath(__file__))
INDEX = os.path.join(HERE, "web", "index.html")

SLOTS = {
    "sc": dict(label="结构连接 SC", hint="(P, P) 或 (N, P, P) 或 (N, E)",
               required=False, accept=".npy,.npz,.mat,.csv,.tsv,.txt"),
    "fc": dict(label="功能连接 FC", hint="(P, P) 或 (N, P, P) 或 (N, E)",
               required=False, accept=".npy,.npz,.mat,.csv,.tsv,.txt"),
    "bold": dict(label="BOLD 时间序列", hint="(N, P, T) 或 (P, T)",
                 required=False, accept=".npy,.npz,.mat,.csv,.tsv,.txt"),
    "networks": dict(label="网络标签", hint="(P,) 整数", required=False,
                     accept=".npy,.npz,.mat,.csv,.tsv,.txt"),
    "hierarchy": dict(label="皮层层级", hint="(P,) 数值", required=False,
                      accept=".npy,.npz,.mat,.csv,.tsv,.txt"),
    "coords": dict(label="球面坐标", hint="(P, 3)，启用 spin 零模型", required=False,
                   accept=".npy,.npz,.mat,.csv,.tsv,.txt"),
    "behaviour": dict(label="行为变量", hint="(N,) 数值", required=False,
                      accept=".npy,.npz,.mat,.csv,.tsv,.txt"),
    "covariates": dict(label="协变量", hint="(N, k) 数值", required=False,
                       accept=".npy,.npz,.mat,.csv,.tsv,.txt"),
}

ALIASES = {
    "sc": ["sc", "SC", "struct", "structural", "sc_matrix"],
    "fc": ["fc", "FC", "func", "functional", "fc_matrix"],
    "bold": ["bold", "BOLD", "timeseries", "ts", "rss"],
    "networks": ["networks", "network", "yeo", "communities", "labels"],
    "hierarchy": ["hierarchy", "hier", "level", "laminar"],
    "coords": ["coords", "xyz", "sphere", "centres", "centers"],
    "behaviour": ["behaviour", "behavior", "behav", "y", "target"],
    "covariates": ["covariates", "cov", "covars", "nuisance"],
}


class Store:
    """One active dataset per server instance."""

    def __init__(self, root):
        self.root = os.path.abspath(root)
        self.upload_dir = os.path.join(self.root, "uploads")
        self.out_dir = os.path.join(self.root, "results")
        os.makedirs(self.upload_dir, exist_ok=True)
        os.makedirs(self.out_dir, exist_ok=True)
        self.lock = threading.Lock()
        self.files = {}
        self.batches = {}
        self.cohort = None
        self.cohort_report = None
        self.bundle = None
        self.job = dict(status="idle", progress=0.0, message="", error=None, result=None)
        self.custom = {}

    # ---------------------------------------------------------- uploads --
    def clear(self):
        with self.lock:
            shutil.rmtree(self.upload_dir, ignore_errors=True)
            os.makedirs(self.upload_dir, exist_ok=True)
            self.files = {}
            self.batches = {}
            self.cohort = None
            self.cohort_report = None
            self.bundle = None
            self.job = dict(status="idle", progress=0.0, message="", error=None, result=None)
        for name in list(self.custom):
            unregister_custom(name)
        self.custom = {}

    def put_file(self, slot, name, payload):
        if slot not in SLOTS:
            raise ValueError("unknown slot %r" % slot)
        safe = os.path.basename(name or "upload.bin")
        ext = os.path.splitext(safe)[1].lower()
        target = os.path.join(self.upload_dir, slot + ext)
        for old in os.listdir(self.upload_dir):
            if os.path.splitext(old)[0] == slot:
                try:
                    os.remove(os.path.join(self.upload_dir, old))
                except OSError:
                    pass
        with open(target, "wb") as f:
            f.write(payload)
        info = self.inspect(slot, target)
        with self.lock:
            self.files[slot] = dict(path=target, name=safe, size=len(payload), info=info)
            self.cohort = None
            self.bundle = None
        return info

    def put_auto(self, name, payload):
        """Store one upload in every slot its contents can serve.

        A single ``.npz``/``.mat`` holding ``sc`` and ``fc`` fills both slots;
        a bare matrix is routed by file-name hint and then by shape.
        """
        safe = os.path.basename(name or "upload.bin")
        ext = os.path.splitext(safe)[1].lower() or ".bin"
        probe = os.path.join(self.upload_dir, "_auto" + ext)
        with open(probe, "wb") as f:
            f.write(payload)
        try:
            obj = io.load_any(probe)
        finally:
            try:
                os.remove(probe)
            except OSError:
                pass
        with self.lock:
            existing = set(self.files)
        slots = _detect_slots(obj, safe, existing)
        if not slots:
            raise ValueError("无法判断文件属于哪个槽位；请手动选择槽位 / Cannot determine the slot for this file; pick a slot manually")
        info = {}
        for slot in slots:
            info[slot] = self.put_file(slot, safe, payload)
        return dict(slots=slots, info=info, dataset=self.summary())

    def inspect(self, slot, path):
        try:
            obj = _load_slot(path, slot)
        except Exception as exc:
            return dict(error=str(exc))
        try:
            arr = np.asarray(obj)
            if slot in ("sc", "fc"):
                vec, p = io.to_edge_vectors(arr, slot)
                return dict(shape=list(arr.shape), n_subjects=int(vec.shape[0]),
                            n_nodes=int(p), n_edges=int(vec.shape[1]))
            if slot == "bold":
                shape = list(arr.shape)
                if len(shape) == 2:
                    n_nodes, n_time = (shape if shape[1] >= shape[0] else shape[::-1])
                    return dict(shape=shape, n_subjects=1, n_nodes=int(n_nodes),
                                n_timepoints=int(n_time))
                return dict(shape=shape, n_subjects=int(shape[0]),
                            n_nodes=int(shape[1]), n_timepoints=int(shape[2]))
            if slot == "coords":
                return dict(shape=list(arr.shape), n_nodes=int(arr.shape[0]))
            return dict(shape=list(arr.shape), n=int(arr.size))
        except Exception as exc:
            return dict(error=str(exc))

    def summary(self):
        with self.lock:
            files = dict(self.files)
        nodes = {k: v.get("info", {}).get("n_nodes") for k, v in files.items()}
        subjects = {k: v.get("info", {}).get("n_subjects") for k, v in files.items()}
        problems = []
        p_vals = {v for v in nodes.values() if v}
        if len(p_vals) > 1:
            detail = ", ".join("%s=%s" % (k, nodes[k]) for k in sorted(nodes) if nodes[k])
            problems.append({"zh": "各文件脑区数不一致：%s" % detail,
                             "en": "Node counts differ across files: %s" % detail})
        s_vals = {v for v in subjects.values() if v}
        if len(s_vals) > 1:
            detail = ", ".join("%s=%s" % (k, subjects[k]) for k in sorted(subjects) if subjects[k])
            problems.append({"zh": "各文件被试数不一致：%s" % detail,
                             "en": "Subject counts differ across files: %s" % detail})
        if s_vals and max(s_vals) < 2:
            problems.append({"zh": "只检测到 1 个被试：变异度、可靠性等统计需要多被试"
                                   "（建议 ≥ 8）；单被试只能算组平均耦合强度",
                             "en": "Only 1 subject detected: variability and reliability "
                                   "statistics need multiple subjects (≥ 8 recommended); "
                                   "a single subject supports coupling values only"})
        if self.cohort is not None and getattr(self.cohort, "n", 0) < 8 and \
                getattr(self.cohort, "n", 0) >= 2:
            problems.append({"zh": "被试数 %d < 8：变异度与个体信号估计噪声较大，结果仅作演示参考"
                                    % self.cohort.n,
                             "en": "n = %d < 8: dispersion and individual-signal estimates are "
                                   "noisy; treat results as a demo" % self.cohort.n})
        return dict(files={k: dict(v) for k, v in files.items()},
                    n_nodes=(p_vals.pop() if len(p_vals) == 1 else None),
                    n_subjects=(s_vals.pop() if len(s_vals) == 1 else None),
                    problems=problems, slots=SLOTS,
                    custom=list_custom())

    # ------------------------------------------------------- batch upload --
    def start_batch(self, slot, name, payload):
        """Stage one subject file for a stacked multi-subject import.

        The browser sends each per-subject file here in turn; when all are
        staged the frontend calls :meth:`finish_batch`, which stacks them
        into one N×P×P (or N×P×T) array and assigns it to the slot.
        """
        if slot not in SLOTS:
            raise ValueError("unknown slot %r" % slot)
        if not payload:
            raise ValueError("空文件 / empty file")
        if len(payload) > 512 * 1024 * 1024:
            raise ValueError("文件超过 512 MB / file exceeds 512 MB")
        safe = os.path.basename(name or "subject.bin")
        ext = os.path.splitext(safe)[1].lower() or ".bin"
        tmp = os.path.join(self.upload_dir, "_batch_%s_%d%s" % (slot, len(payload), ext))
        with open(tmp, "wb") as f:
            f.write(payload)
        try:
            obj = io.load_any(tmp)
        finally:
            try:
                os.remove(tmp)
            except OSError:
                pass
        if isinstance(obj, dict):
            # a per-subject .npz/.mat may hold several arrays; take the first 2-D one
            arrs = [np.asarray(v) for v in obj.values()
                    if isinstance(v, np.ndarray) and v.ndim >= 2]
            if not arrs:
                raise ValueError("%s: 未找到二维数组 / no 2-D array found" % safe)
            obj = arrs[0]
        arr = np.asarray(obj, dtype=np.float64)
        if slot == "bold" and arr.ndim == 2 and arr.shape[0] > arr.shape[1]:
            arr = arr.T  # time x nodes -> nodes x time
        if arr.ndim < 2:
            raise ValueError("%s: 期望二维矩阵 / expected a 2-D matrix" % safe)
        with self.lock:
            batch = self.batches.setdefault(slot, [])
            if batch and batch[0].shape != arr.shape:
                raise ValueError(
                    "%s 形状 %s 与已导入文件形状 %s 不一致 / %s has shape %s, but the "
                    "staged files have shape %s"
                    % (safe, list(arr.shape), list(batch[0].shape),
                       safe, list(arr.shape), list(batch[0].shape)))
            batch.append(arr)
            n = len(batch)
        return dict(slot=slot, staged=n, shape=list(arr.shape),
                    name=os.path.basename(safe))

    def finish_batch(self, slot):
        with self.lock:
            batch = self.batches.get(slot, [])
            if len(batch) < 2:
                self.batches.pop(slot, None)
                raise ValueError(
                    "批量导入至少需要 2 个同形状文件 / batch import needs at least "
                    "2 files of the same shape")
            stack = np.stack(batch, axis=0)
            self.batches.pop(slot, None)
        stamp = uuid.uuid4().hex[:8]
        name = "batch_%d_subjects_%s.npy" % (stack.shape[0], stamp)
        target = os.path.join(self.upload_dir, slot + "_batch.npy")
        np.save(target, stack)
        with open(target, "rb") as f:
            payload = f.read()
        os.remove(target)
        self.put_file(slot, name, payload)
        return dict(slot=slot, subjects=int(stack.shape[0]),
                    dataset=self.summary())

    def cancel_batch(self, slot):
        with self.lock:
            self.batches.pop(slot, None)
        return dict(slot=slot, cancelled=True)

    # ----------------------------------------------------------- build --
    def build(self, fc_kind="pearson", time_axis="auto", state="rest", n_networks=None):
        with self.lock:
            files = dict(self.files)
        if "bold" not in files and "fc" not in files:
            raise ValueError("请先导入 FC 矩阵或 BOLD 时间序列 / "
                             "Import an FC matrix or a BOLD series first")
        if "sc" not in files and "bold" not in files:
            raise ValueError("请至少导入结构连接矩阵 SC / "
                             "Import at least a structural connectivity (SC) matrix")

        def get(slot):
            if slot not in files:
                return None
            return _load_slot(files[slot]["path"], slot)

        sc = get("sc")
        fc = get("fc")
        bold = get("bold")
        networks = get("networks")
        hierarchy = get("hierarchy")
        coords = get("coords")
        behaviour = get("behaviour")
        covariates = get("covariates")
        if networks is not None:
            networks = np.asarray(networks).ravel().astype(int)
        if hierarchy is not None:
            hierarchy = np.asarray(hierarchy, dtype=float).ravel()
        if behaviour is not None:
            behaviour = np.asarray(behaviour, dtype=float).ravel()
        if coords is not None:
            coords = np.asarray(coords, dtype=float).reshape(-1, 3)

        coh, rep = io.build_cohort(
            sc=sc, fc=fc, bold=bold, networks=networks, hierarchy=hierarchy,
            coords=coords, behaviour=behaviour, covariates=covariates,
            state=state, fc_kind=fc_kind, time_axis=time_axis,
            n_networks=n_networks)
        if getattr(coh, "n", 0) < 2:
            rep.setdefault("warnings", []).append({
                "zh": "当前数据只有 %d 个被试：变异度统计不可用，请使用「耦合值计算」模式；"
                      "若要跑变异分析，请上传多被试数据（SC/FC 为 n_subjects × P × P，"
                      "BOLD 为 n_subjects × P × T）" % getattr(coh, "n", 0),
                "en": "This dataset has %d subject(s): variability statistics are not "
                      "available — use the coupling-values mode; for variability analysis "
                      "upload multi-subject data (SC/FC as n_subjects × P × P, BOLD as "
                      "n_subjects × P × T)" % getattr(coh, "n", 0)})
        with self.lock:
            self.cohort = coh
            self.cohort_report = rep
        return coh, rep

    # ------------------------------------------------------------- job --
    def start(self, settings):
        with self.lock:
            if self.job.get("status") == "running":
                raise ValueError("已有分析在运行 / An analysis is already running")
            self.job = dict(status="running", progress=0.0, message="准备数据",
                            error=None, result=None)
        t = threading.Thread(target=self._run_job, args=(settings,), daemon=True)
        t.start()
        return dict(status="running")

    def _run_job(self, settings):
        try:
            for c in settings.get("custom") or []:
                name, source = c.get("name"), c.get("source")
                if not name or not source:
                    continue
                register_custom(name, source, c.get("description", ""))
                self.custom[name] = source

            coh, rep = self.build(
                fc_kind=settings.get("fc_kind", "pearson"),
                time_axis=settings.get("time_axis", "auto"),
                state=settings.get("state", "rest"),
                n_networks=settings.get("n_networks"))
            self._progress(0.02, "data: %d subjects × %d regions · 已载入 %d 被试 × %d 脑区"
                           % (coh.n, coh.n_nodes, coh.n, coh.n_nodes))

            metrics = settings.get("metrics") or ["communicability"]
            mode = settings.get("mode", "variability")

            if mode == "magnitude":
                bundle = engine.magnitude(coh, metrics, progress=self._progress)
                files = engine.export_coupling(
                    bundle, self.out_dir, fmt=("csv", "npz"))
                import zipfile
                zip_path = os.path.join(self.out_dir, "coupling_csv.zip")
                with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
                    for f in files:
                        if f.endswith(".csv"):
                            zf.write(f, os.path.basename(f))
                with self.lock:
                    self.bundle = bundle
                payload = _magnitude_for_web(bundle, self.out_dir)
                payload["downloads"] = [
                    dict(kind="coupling_npz", label="耦合值 .npz",
                         url="/api/download?kind=coupling_npz"),
                    dict(kind="coupling_csv", label="逐被试 CSV 打包",
                         url="/api/download?kind=coupling_csv"),
                ]
                payload["warnings"] = list(rep.get("warnings", []))
                with self.lock:
                    self.job.update(status="done", progress=1.0, message="完成",
                                    result=payload, error=None)
                return

            if coh.n < 2:
                raise ValueError(
                    "变异分析需要 ≥2 名被试（当前 %d）。单个被试请切换到「耦合值计算」模式；"
                    "多被试数据请上传 n_subjects × P × P 的 SC/FC 或 n_subjects × P × T 的 BOLD。"
                    " / Variability analysis needs ≥2 subjects (got %d). For a single subject "
                    "switch to the coupling-values mode; multi-subject data: SC/FC as "
                    "n_subjects × P × P or BOLD as n_subjects × P × T." % (coh.n, coh.n))
            estimators = tuple(settings.get("estimators") or ["regional", "edge"])
            n_perm = int(settings.get("n_perm", 500))
            null = settings.get("null", "auto")
            seed = int(settings.get("seed", 0))
            label = settings.get("behaviour_label") or "behaviour"

            bundle = engine.analyze(
                coh, metrics, estimators=estimators, n_perm=n_perm, null=null,
                seed=seed, outdir=self.out_dir, make_figures=True,
                behaviour_label=label, progress=self._progress,
                do_sweep=bool(settings.get("sweep", True)))

            html = report.build_bundle(
                bundle, os.path.join(self.out_dir, "report.html"),
                title="SFCX 分析报告")

            # export the cohort in both interchange formats
            io.export_npz(os.path.join(self.out_dir, "cohort.npz"), coh)
            io.export_mat(os.path.join(self.out_dir, "cohort.mat"), coh)

            with self.lock:
                self.bundle = bundle
            payload = _bundle_for_web(bundle, self.out_dir)
            payload["report_url"] = "/api/download?kind=html"
            with self.lock:
                self.job.update(status="done", progress=1.0, message="完成",
                                result=payload, error=None)
        except CustomCodeError as exc:
            self._fail("自定义函数错误：%s" % exc)
        except Exception as exc:  # noqa: BLE001 - surfaced to the user
            self._fail("%s: %s" % (type(exc).__name__, exc))
            traceback.print_exc()

    def _progress(self, frac, msg):
        with self.lock:
            self.job["progress"] = float(frac)
            self.job["message"] = str(msg)

    def _fail(self, msg):
        with self.lock:
            self.job.update(status="error", message="失败", error=msg, result=None)


def _detect_slots(obj, name, existing):
    base = os.path.basename(name or "").lower()
    if isinstance(obj, dict):
        hits = []
        for slot in SLOTS:
            if any(key in obj for key in ALIASES.get(slot, [slot])):
                hits.append(slot)
        if hits:
            return hits
    arr = np.asarray(obj)
    hint = lambda *words: any(w in base for w in words)  # noqa: E731
    if arr.ndim == 2:
        if arr.shape[0] == arr.shape[1]:
            if hint("sc", "struct"):
                return ["sc"]
            if hint("fc", "func", "corr", "conn"):
                return ["fc"]
            return ["fc"] if "sc" in existing else ["sc"]
        return ["bold"]
    if arr.ndim == 3:
        if hint("bold", "timeseries", "ts_"):
            return ["bold"]
        if arr.shape[1] == arr.shape[2]:
            return ["fc"] if "sc" in existing else ["sc"]
        return ["bold"]
    if arr.ndim == 1:
        if hint("net", "yeo", "label", "comm"):
            return ["networks"]
        if hint("hier", "level", "grad"):
            return ["hierarchy"]
        if hint("behav", "target"):
            return ["behaviour"]
        return ["hierarchy"] if "networks" in existing else ["networks"]
    return []


def _load_slot(path, slot):
    obj = io.load_any(path)
    if isinstance(obj, dict):
        for key in ALIASES.get(slot, [slot]):
            if key in obj:
                return obj[key]
        if len(obj) == 1:
            return next(iter(obj.values()))
        arrays = [v for v in obj.values() if isinstance(v, np.ndarray)]
        if len(arrays) == 1:
            return arrays[0]
        raise ValueError("文件中有多个数组，找不到 %r；可用键：%s"
                         % (slot, ", ".join(sorted(obj))))
    return obj


def _magnitude_for_web(bundle, out_dir):
    """JSON-safe summary of a magnitude bundle (full fields stay on disk)."""
    out = {"mode": "magnitude",
           "meta": bundle.get("meta", {}),
           "runtime_s": bundle.get("runtime_s"),
           "metrics": {}}
    for name, m in bundle.get("metrics", {}).items():
        entry = {"kind": m.get("kind"), "subject_ids": m.get("subject_ids")}
        reg = np.asarray(m.get("regional"), dtype=np.float64)
        entry["regional_shape"] = list(reg.shape)
        # keep at most 60 subjects x 400 parcels in the browser payload
        entry["regional_preview"] = reg[:60, :400].tolist()
        mm = m.get("mean_matrix")
        if mm is not None:
            mm = np.asarray(mm, dtype=np.float64)
            entry["mean_matrix_preview"] = mm[:400, :400].tolist()
            entry["mean_matrix_url"] = "/api/file/coupling_values_mean_%s.csv" % name
            path = os.path.join(out_dir, "coupling_values_mean_%s.csv" % name)
            with open(path, "w", encoding="utf-8", newline="") as fh:
                for row in mm:
                    fh.write(",".join("%.6g" % v for v in row) + "\n")
        else:
            entry["mean_matrix_preview"] = None
        entry["note"] = ("regional output only; no P×P matrix"
                         if mm is None else "")
        out["metrics"][name] = entry
    return out


def _bundle_for_web(bundle, out_dir):
    """Replace figure file paths with URLs the browser can fetch."""
    out = json.loads(json.dumps(bundle, ensure_ascii=False))
    for block in out.get("metrics", {}).values():
        figs = block.get("figures") or {}
        block["figures"] = {}
        for k, p in figs.items():
            if not p:
                continue
            rel = os.path.relpath(p, out_dir).replace("\\", "/")
            block["figures"][k] = dict(name=k, url="/api/file/" + rel)
    sw = out.get("sweep") or {}
    if sw.get("figure"):
        rel = os.path.relpath(sw["figure"], out_dir).replace("\\", "/")
        sw["figure"] = dict(name="coordinate_matrix", url="/api/file/" + rel)
    out["downloads"] = [
        dict(kind="html", label="分析报告 HTML", url="/api/download?kind=html"),
        dict(kind="mat", label="MATLAB 数据 .mat", url="/api/download?kind=mat"),
        dict(kind="npz", label="Python 数据 .npz", url="/api/download?kind=npz"),
    ]
    return out


# ------------------------------------------------------------- handler ----

class Handler(BaseHTTPRequestHandler):
    hub = None
    public = False
    server_version = "SFCX/0.3"

    def log_message(self, fmt, *args):  # quieter console
        if os.environ.get("SFCX_WEB_VERBOSE"):
            super().log_message(fmt, *args)

    # ---------------------------------------------------------- session --
    def _resolve_session(self):
        """Attach the per-visitor Store and mark the response if a fresh
        session id was issued (the cookie is written by ``_json``)."""
        sid = None
        raw = self.headers.get("Cookie") or ""
        for part in raw.split(";"):
            k, _, v = part.strip().partition("=")
            if k == "sfcx_sid" and v and self.hub.SID_RE.fullmatch(v):
                sid = v
        self._sid_new = False
        if sid is None:
            sid = uuid.uuid4().hex[:16]
            self._sid_new = True
        self.store = self.hub.get(sid)
        self._sid = sid
        if os.environ.get("SFCX_SESSION_LOG"):
            try:
                with open(os.environ["SFCX_SESSION_LOG"], "a", encoding="utf-8") as fh:
                    fh.write("%s %s %s files=%s\n" % (
                        time.strftime("%H:%M:%S"), self.command, self.path.split("?")[0],
                        sorted(self.store.files.keys())))
            except OSError:
                pass

    # ------------------------------------------------------------ utils --
    def _json(self, payload, status=HTTPStatus.OK):
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        if getattr(self, "_sid_new", False):
            self.send_header("Set-Cookie",
                             "sfcx_sid=%s; Path=/; Max-Age=604800; SameSite=Lax"
                             % self._sid)
        self.end_headers()
        self.wfile.write(body)

    def _error(self, message, status=HTTPStatus.BAD_REQUEST):
        self._json(dict(error=str(message)), status=status)

    def _body(self):
        length = int(self.headers.get("Content-Length") or 0)
        if length <= 0:
            return b""
        return self.rfile.read(length)

    def _json_body(self):
        raw = self._body()
        if not raw:
            return {}
        return json.loads(raw.decode("utf-8"))

    # ------------------------------------------------------------- routes --
    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        query = urllib.parse.parse_qs(parsed.query)
        try:
            self._resolve_session()
            if path in ("/", "/index.html"):
                return self._serve_index()
            if path == "/api/catalog":
                return self._json(_catalog())
            if path == "/api/dataset":
                return self._json(self.store.summary())
            if path == "/api/job":
                return self._json(self.store.job)
            if path == "/api/download":
                return self._download(query.get("kind", ["html"])[0])
            if path.startswith("/api/file/"):
                return self._serve_file(path[len("/api/file/"):])
            return self._error("not found", HTTPStatus.NOT_FOUND)
        except Exception as exc:  # noqa: BLE001
            self._error("%s: %s" % (type(exc).__name__, exc), HTTPStatus.INTERNAL_SERVER_ERROR)

    def do_POST(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        query = urllib.parse.parse_qs(parsed.query)
        try:
            self._resolve_session()
            if path == "/api/upload":
                slot = query.get("slot", [""])[0]
                name = query.get("name", ["upload.bin"])[0]
                payload = self._body()
                if not payload:
                    return self._error("空文件 / empty file")
                if len(payload) > 512 * 1024 * 1024:
                    return self._error("文件超过 512 MB / file exceeds 512 MB")
                if slot == "auto":
                    return self._json(self.store.put_auto(name, payload))
                info = self.store.put_file(slot, name, payload)
                return self._json(dict(slot=slot, info=info, dataset=self.store.summary()))
            if path == "/api/batch/upload":
                slot = query.get("slot", [""])[0]
                name = query.get("name", ["subject.bin"])[0]
                return self._json(self.store.start_batch(slot, name, self._body()))
            if path == "/api/batch/finish":
                slot = query.get("slot", [""])[0]
                return self._json(self.store.finish_batch(slot))
            if path == "/api/batch/cancel":
                slot = query.get("slot", [""])[0]
                return self._json(self.store.cancel_batch(slot))
            if path == "/api/demo":
                return self._json(self._demo())
            if path == "/api/custom/validate":
                if self.public:
                    return self._error("公开部署模式下已禁用自定义代码", HTTPStatus.FORBIDDEN)
                body = self._json_body()
                try:
                    res = dry_run(body.get("source", ""), n_subjects=8, n_nodes=24)
                    return self._json(dict(ok=True, **res))
                except CustomCodeError as exc:
                    return self._json(dict(ok=False, error=str(exc)))
            if path == "/api/run":
                body = self._json_body()
                if self.public:
                    body.pop("custom", None)
                return self._json(self.store.start(body))
            if path == "/api/reset":
                self.store.clear()
                return self._json(dict(ok=True))
            return self._error("not found", HTTPStatus.NOT_FOUND)
        except CustomCodeError as exc:
            self._error(str(exc))
        except ValueError as exc:
            self._error(str(exc), HTTPStatus.BAD_REQUEST)
        except Exception as exc:  # noqa: BLE001
            self._error("%s: %s" % (type(exc).__name__, exc), HTTPStatus.INTERNAL_SERVER_ERROR)

    # ------------------------------------------------------------ serving --
    def _serve_index(self):
        if not os.path.exists(INDEX):
            return self._error("index.html missing", HTTPStatus.NOT_FOUND)
        with open(INDEX, "rb") as f:
            body = f.read()
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        if getattr(self, "_sid_new", False):
            self.send_header("Set-Cookie",
                             "sfcx_sid=%s; Path=/; Max-Age=604800; SameSite=Lax"
                             % self._sid)
        self.end_headers()
        self.wfile.write(body)

    def _serve_file(self, rel):
        rel = urllib.parse.unquote(rel).replace("\\", "/")
        target = os.path.abspath(os.path.join(self.store.out_dir, rel))
        if not target.startswith(os.path.abspath(self.store.out_dir)):
            return self._error("forbidden", HTTPStatus.FORBIDDEN)
        if not os.path.exists(target):
            return self._error("missing file", HTTPStatus.NOT_FOUND)
        ctype = mimetypes.guess_type(target)[0] or "application/octet-stream"
        with open(target, "rb") as f:
            body = f.read()
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _download(self, kind):
        names = {"html": ("report.html", "sfcx_report.html", "text/html"),
                 "mat": ("cohort.mat", "sfcx_cohort.mat", "application/octet-stream"),
                 "npz": ("cohort.npz", "sfcx_cohort.npz", "application/octet-stream"),
                 "coupling_npz": ("coupling_values.npz", "sfcx_coupling_values.npz",
                                  "application/octet-stream"),
                 "coupling_csv": ("coupling_csv.zip", "sfcx_coupling_csv.zip",
                                  "application/zip")}
        if kind not in names:
            return self._error("unknown download kind")
        fname, out_name, ctype = names[kind]
        target = os.path.join(self.store.out_dir, fname)
        if not os.path.exists(target):
            return self._error("请先运行一次分析", HTTPStatus.NOT_FOUND)
        with open(target, "rb") as f:
            body = f.read()
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Content-Disposition", 'attachment; filename="%s"' % out_name)
        self.end_headers()
        self.wfile.write(body)

    def _demo(self):
        from .demo import make_cohort
        body = self._json_body()
        n_subjects = int(body.get("n_subjects", 240))
        n_nodes = int(body.get("n_nodes", 120))
        seed = int(body.get("seed", 0))
        coh = make_cohort(n_subjects=n_subjects, n_nodes=n_nodes, seed=seed)
        self.store.clear()
        npz = os.path.join(self.store.upload_dir, "demo.npz")
        io.export_npz(npz, coh)
        for slot in ("sc", "fc"):
            self.store.put_file(slot, "demo.npz", open(npz, "rb").read())
        for slot, key in (("networks", "networks"), ("hierarchy", "hierarchy"),
                          ("coords", "coords")):
            arr = getattr(coh, key)
            path = os.path.join(self.store.upload_dir, slot + ".npy")
            np.save(path, arr)
            self.store.put_file(slot, slot + ".npy", open(path, "rb").read())
        return dict(dataset=self.store.summary())


def _catalog():
    cat = catalog()
    cat["api_version"] = 2
    cat["slots"] = SLOTS
    cat["aliases"] = ALIASES
    cat["limits"] = dict(max_upload_mb=512, default_perm=500, max_perm=5000,
                         default_subjects=240, default_nodes=120)
    cat["custom_template"] = (
        "def couple(sc, fc, iu, iv, n_nodes):\n"
        "    \"\"\"sc, fc: (n_subjects, n_edges); iu, iv: edge endpoints.\"\"\"\n"
        "    zs = zscore_edges(sc, iu, n_nodes)\n"
        "    zf = zscore_edges(fc, iu, n_nodes)\n"
        "    return zs * zf\n")
    cat["notes"] = [
        "SC 与 FC 需为同一脑区划分、同一上三角顺序；本工具会自动展开矩阵。",
        "BOLD 会按所选 pipeline（pearson / precision）转成 FC 边向量。",
        "没有提供网络标签或层级时，会给出代理值并在结果中标注。",
        "p 值必须与零模型一起报告；无坐标时自动退化为 msr。",
    ]
    return cat


class SessionHub:
    """One :class:`Store` per browser session (cookie ``sfcx_sid``).

    Lets several people use the same deployment without touching each
    other's uploads or results.  Session ids are validated as 16 hex
    characters before they ever touch the filesystem.
    """

    SID_RE = re.compile(r"^[0-9a-f]{16}$")
    MAX_SESSIONS = 200

    def __init__(self, base):
        self.base = os.path.abspath(base)
        os.makedirs(self.base, exist_ok=True)
        self.lock = threading.Lock()
        self.stores = {}

    def get(self, sid):
        with self.lock:
            st = self.stores.get(sid)
            if st is None:
                if len(self.stores) >= self.MAX_SESSIONS:
                    raise ValueError("服务器会话数已达上限，请稍后再试 / server session limit reached, try again later")
                st = Store(os.path.join(self.base, "s_" + sid))
                self.stores[sid] = st
            return st


def serve(host="127.0.0.1", port=8765, root=None, open_browser=False, public=False):
    root = root or os.path.join(os.getcwd(), "sfcx_web_sessions")
    os.makedirs(root, exist_ok=True)
    Handler.hub = SessionHub(root)
    Handler.public = bool(public)
    httpd = ThreadingHTTPServer((host, int(port)), Handler)
    url = "http://%s:%d/" % (host, int(port))
    print("SFCX web UI  ->  %s" % url)
    print("session root ->  %s" % root)
    if public:
        print("public mode  ->  per-visitor sessions, custom code disabled")
    if open_browser:
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nbye")
    finally:
        httpd.server_close()


def main(argv=None):
    p = argparse.ArgumentParser(prog="sfcx-web", description="SFCX web interface")
    p.add_argument("--host", default="127.0.0.1")
    p.add_argument("--port", type=int, default=8765)
    p.add_argument("--root", default=None, help="session directory")
    p.add_argument("--open", action="store_true", help="open the browser")
    p.add_argument("--public", action="store_true",
                   help="per-visitor sessions; disable custom-code execution")
    args = p.parse_args(argv)

    # cloud/sandbox deployments inject PORT and expect 0.0.0.0
    env_port = os.environ.get("PORT")
    host = args.host
    if env_port:
        host = "0.0.0.0"
    port = int(env_port) if env_port else args.port
    serve(host, port, args.root, args.open, args.public)


if __name__ == "__main__":
    main()
