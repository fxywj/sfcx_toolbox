# -*- coding: utf-8 -*-
"""Self-contained HTML report (figures embedded, no external assets)."""

from __future__ import annotations

import base64
import html
import os


def _img(path, caption=""):
    if not path or not os.path.exists(path):
        return ""
    with open(path, "rb") as f:
        b64 = base64.b64encode(f.read()).decode("ascii")
    cap = '<figcaption>%s</figcaption>' % html.escape(caption) if caption else ""
    return ('<figure><img src="data:image/png;base64,%s" alt="%s"/>%s</figure>'
            % (b64, html.escape(caption), cap))


def _row(k, v, cls=""):
    return '<tr class="%s"><th>%s</th><td>%s</td></tr>' % (cls, html.escape(str(k)), v)


def _fmt(x, nd=3):
    try:
        return ("%." + str(nd) + "f") % float(x)
    except Exception:
        return str(x)


CSS = """
:root{--ink:#1f3a5f;--teal:#2a9d8f;--orange:#e76f51;--grey:#8d99ae;--bg:#ffffff;--panel:#f6f8fa;--line:#d8dee4}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.65 -apple-system,"Segoe UI",Roboto,"Helvetica Neue",Arial,"Microsoft YaHei",sans-serif}
.wrap{max-width:980px;margin:0 auto;padding:40px 28px 80px}
h1{font-size:25px;line-height:1.3;margin:0 0 6px}
h2{font-size:16px;margin:38px 0 12px;padding-bottom:6px;border-bottom:1px solid var(--line)}
.sub{color:var(--grey);font-size:13px;margin-bottom:26px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:16px 18px;margin:14px 0}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin:16px 0}
.kpi{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:14px 16px}
.kpi .n{font-size:23px;font-weight:600;letter-spacing:-.3px;display:block}
.kpi .l{font-size:12px;color:var(--grey)}
table{width:100%;border-collapse:collapse;font-size:13px}
th,td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--line);vertical-align:top}
tr>th{width:34%;font-weight:600}
figure{margin:16px 0}
img{width:100%;border:1px solid var(--line);border-radius:8px;background:#fff}
figcaption{color:var(--grey);font-size:12px;margin-top:6px}
code{background:#eef1f4;padding:1px 5px;border-radius:4px;font-size:12.5px}
.note{border-left:3px solid var(--orange);background:#fff7f4;padding:10px 14px;border-radius:0 8px 8px 0;font-size:13px}
footer{margin-top:44px;color:var(--grey);font-size:12px;border-top:1px solid var(--line);padding-top:14px}
"""


def build(res, out_html, title="SFCX route report", sweep=None, notes=None):
    c = res.get("coordinate", {})
    data = res.get("data", {})
    kpi = []
    ax = res.get("axis")
    if ax:
        kpi.append('<div class="kpi"><span class="n">%s</span><span class="l">axis r (hierarchy)</span></div>'
                   % _fmt(ax.get("r")))
        kpi.append('<div class="kpi"><span class="n">%s</span><span class="l">axis p (%s)</span></div>'
                   % (_fmt(ax.get("p")), html.escape(str(ax.get("null")))))
    lay = res.get("layers", {})
    if lay:
        kpi.append('<div class="kpi"><span class="n">%s</span><span class="l">within - between</span></div>'
                   % _fmt(lay.get("delta")))
    rel = res.get("reliability", {})
    if rel:
        kpi.append('<div class="kpi"><span class="n">%s</span><span class="l">split-half r</span></div>'
                   % _fmt(rel.get("split_half_r")))
    epi = res.get("epicentres")
    if epi:
        kpi.append('<div class="kpi"><span class="n">%d</span><span class="l">focal epicentres</span></div>'
                   % epi.get("n", 0))

    parts = []
    parts.append('<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"/>')
    parts.append('<meta name="viewport" content="width=device-width, initial-scale=1"/>')
    parts.append('<title>%s</title><style>%s</style></head><body><div class="wrap">'
                 % (html.escape(title), CSS))
    parts.append("<h1>%s</h1>" % html.escape(title))
    parts.append('<div class="sub">%s &nbsp;|&nbsp; n = %s subjects &nbsp;|&nbsp; %s regions &nbsp;|&nbsp; '
                 '%s edges &nbsp;|&nbsp; state: %s</div>'
                 % (html.escape(str(c.get("label", ""))), data.get("n_subjects"),
                    data.get("n_nodes"), data.get("n_edges"), html.escape(str(data.get("state")))))

    if kpi:
        parts.append('<div class="grid">%s</div>' % "".join(kpi))

    parts.append("<h2>Coordinate</h2>")
    parts.append('<table>')
    parts.append(_row("Definition (layer 1)", "<code>%s</code> &mdash; %s"
                      % (html.escape(str(c.get("metric"))), html.escape(str(c.get("metric_description"))))))
    parts.append(_row("Functional pipeline (layer 2)", "<code>%s</code>" % html.escape(str(c.get("pipeline")))))
    parts.append(_row("Dispersion estimator (layer 3)", "<code>%s</code>" % html.escape(str(c.get("estimator")))))
    parts.append(_row("Observation scale (layer 4)", "<code>%s</code>" % html.escape(str(c.get("scale")))))
    parts.append(_row("Null model", "%s, %s draws" % (html.escape(str(res.get("null_model"))),
                                                      (ax or {}).get("n_perm", "-"))))
    parts.append('</table>')

    if notes:
        parts.append('<div class="note">%s</div>' % notes)

    figs = res.get("figures", {})
    if figs:
        parts.append("<h2>Figures</h2>")
        caps = {
            "route_card": "Active coordinate card.",
            "spectrum": "Descending dispersion spectrum and its distribution.",
            "epicentres": "Network-ordered dispersion strip; triangles mark focal epicentres "
                          "(ordering, not cortical anatomy).",
            "layers": "Within- versus between-network mean edge dispersion.",
            "reliability": "Split-half agreement of the edge-level dispersion map.",
            "individual": "Individual deviation-spectrum slope against the behavioural target.",
        }
        for k in ("route_card", "spectrum", "epicentres", "layers", "reliability", "individual"):
            if k in figs:
                parts.append(_img(figs[k], caps.get(k, k)))

    if sweep:
        parts.append("<h2>Coordinate sweep</h2>")
        parts.append("<p>Same data, same subjects, same null. Only the estimator changes. "
                     "Cells give Spearman r between the dispersion map and the cortical hierarchy; "
                     "* marks p &lt; 0.05.</p>")
        rows = sweep.get("rows", [])
        cols = sweep.get("cols", [])
        r = sweep.get("r")
        p = sweep.get("p")
        parts.append('<table><tr><th style="width:26%%">definition \\ estimator</th>%s</tr>'
                     % "".join("<th>%s</th>" % html.escape(str(x)) for x in cols))
        for i, m in enumerate(rows):
            cells = []
            for j, _ in enumerate(cols):
                try:
                    star = "*" if p is not None and p[i, j] < 0.05 else ""
                    cells.append("<td>%s%s</td>" % (_fmt(r[i, j], 2), star))
                except Exception:
                    cells.append("<td>-</td>")
            parts.append("<tr><th>%s</th>%s</tr>" % (html.escape(str(m)), "".join(cells)))
        parts.append("</table>")
        if sweep.get("figure"):
            parts.append(_img(sweep["figure"], "Axis strength across definitions and estimators."))

    if res.get("individual"):
        parts.append("<h2>Observation scale</h2>")
        parts.append("<p>The same deviation-spectrum slope, read at four scales. "
                     "Only the edge scale is expected to carry individual signal.</p>")
        parts.append("<table><tr><th>scale</th><th>r</th><th>p</th><th>n</th></tr>")
        for k, st in res["individual"].items():
            if k == "best_scale":
                continue
            parts.append("<tr><th>%s</th><td>%s</td><td>%s</td><td>%s</td></tr>"
                         % (html.escape(k), _fmt(st.get("r")), _fmt(st.get("p"), 4), st.get("n")))
        parts.append("</table>")

    parts.append('<footer>Generated by <code>sfcx</code>. The null model matters: '
                 'spin and Moran spectral randomisation preserve spatial autocorrelation, '
                 'plain permutation does not and gives anti-conservative p values. '
                 'Report the null with the p.</footer>')
    parts.append("</div></body></html>")

    os.makedirs(os.path.dirname(os.path.abspath(out_html)), exist_ok=True)
    with open(out_html, "w", encoding="utf-8") as f:
        f.write("".join(parts))
    return out_html


# ------------------------------------------------------ bundle report ------

def _fmt_or_dash(x, nd=3):
    if x is None:
        return "&mdash;"
    return _fmt(x, nd)


def _pct(x):
    return "&mdash;" if x is None else ("%.3g" % float(x))


def build_bundle(bundle, out_html, title="SFCX analysis report", notes=None):
    """Render the full multi-metric bundle produced by :func:`sfcx.analyze`."""
    meta = bundle.get("meta", {})
    sweep = bundle.get("sweep") or {}
    parts = ['<!DOCTYPE html><html lang="zh"><head><meta charset="utf-8"/>',
             '<meta name="viewport" content="width=device-width, initial-scale=1"/>',
             '<title>%s</title><style>%s</style></head><body><div class="wrap">'
             % (html.escape(title), CSS)]
    parts.append("<h1>%s</h1>" % html.escape(title))
    parts.append('<div class="sub">n = %s subjects &nbsp;|&nbsp; %s regions &nbsp;|&nbsp; '
                 '%s edges &nbsp;|&nbsp; %s networks &nbsp;|&nbsp; state: %s &nbsp;|&nbsp; '
                 'null: %s (%s draws)</div>'
                 % (meta.get("n_subjects"), meta.get("n_nodes"), meta.get("n_edges"),
                    meta.get("n_networks"), html.escape(str(meta.get("state"))),
                    html.escape(str(meta.get("null_model"))), meta.get("n_perm")))

    for w in bundle.get("warnings", []):
        if isinstance(w, dict):
            w = w.get("en") or w.get("zh") or str(w)
        parts.append('<div class="note">%s</div>' % html.escape(str(w)))

    if notes:
        parts.append('<div class="note">%s</div>' % notes)

    if sweep:
        rows, cols = sweep.get("rows", []), sweep.get("cols", [])
        parts.append("<h2>Coordinate sweep</h2>")
        parts.append("<p>Spearman r between the dispersion map and the cortical hierarchy, "
                     "for every definition &times; estimator cell. * marks p &lt; 0.05.</p>")
        parts.append('<table><tr><th style="width:26%%">definition \\ estimator</th>%s</tr>'
                     % "".join("<th>%s</th>" % html.escape(str(x)) for x in cols))
        for i, m in enumerate(rows):
            cells = []
            for j in range(len(cols)):
                r = sweep["r"][i][j]
                p = sweep["p"][i][j]
                star = "*" if p is not None and p < 0.05 else ""
                cells.append("<td>%s%s</td>" % (_fmt_or_dash(r, 2), star))
            parts.append("<tr><th>%s</th>%s</tr>" % (html.escape(str(m)), "".join(cells)))
        parts.append("</table>")
        if sweep.get("figure"):
            parts.append(_img(sweep["figure"], "Axis strength across definitions and estimators."))

    for metric, block in bundle.get("metrics", {}).items():
        parts.append("<h2>%s <code>%s</code></h2>" % (html.escape(block.get("label", metric)),
                                                       html.escape(metric)))
        if block.get("description"):
            parts.append('<div class="sub" style="margin:-6px 0 12px">%s</div>'
                         % html.escape(block["description"]))
        kpi = []
        for est, item in block.get("estimators", {}).items():
            ax = item.get("axis", {})
            kpi.append('<div class="kpi"><span class="n">%s</span>'
                       '<span class="l">axis r &middot; %s (%s)</span></div>'
                       % (_fmt_or_dash(ax.get("r")),
                          html.escape(str(est)), _pct(ax.get("p"))))
        lay = block.get("layers", {})
        kpi.append('<div class="kpi"><span class="n">%s</span>'
                   '<span class="l">within &minus; between (p=%s)</span></div>'
                   % (_fmt_or_dash(lay.get("delta")), _pct(lay.get("p"))))
        rel = block.get("reliability", {})
        kpi.append('<div class="kpi"><span class="n">%s</span>'
                   '<span class="l">split-half r</span></div>'
                   % _fmt_or_dash(rel.get("split_half_r")))
        parts.append('<div class="grid">%s</div>' % "".join(kpi))

        figs = block.get("figures", {})
        caps = {"route_card": "Active coordinate card.",
                "spectrum": "Descending dispersion spectrum and its distribution.",
                "epicentres": "Network-ordered dispersion strip; triangles mark focal epicentres.",
                "layers": "Within- versus between-network mean edge dispersion.",
                "reliability": "Split-half agreement of the edge-level dispersion map.",
                "individual": "Individual deviation-spectrum slope against the behavioural target."}
        for k in ("route_card", "spectrum", "epicentres", "layers", "reliability", "individual"):
            if k in figs:
                parts.append(_img(figs[k], caps.get(k, k)))

        indiv = block.get("individual") or {}
        if indiv:
            parts.append("<h3>Deviation-spectrum slopes</h3>")
            parts.append("<table><tr><th>scale</th><th>mean</th><th>r</th><th>p</th><th>n</th></tr>")
            for k, st in indiv.items():
                if k == "best_scale":
                    continue
                sm = st.get("summary", {})
                parts.append("<tr><th>%s</th><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>"
                             % (html.escape(str(k)), _fmt_or_dash(sm.get("mean")),
                                _fmt_or_dash(st.get("r")), _pct(st.get("p")), st.get("n", "&mdash;")))
            parts.append("</table>")

    parts.append('<footer>Generated by <code>sfcx</code>. p values are read against the '
                 'declared null model; a p without its null is not a result.</footer>')
    parts.append("</div></body></html>")
    os.makedirs(os.path.dirname(os.path.abspath(out_html)), exist_ok=True)
    with open(out_html, "w", encoding="utf-8") as f:
        f.write("".join(parts))
    return out_html
