# -*- coding: utf-8 -*-
"""Custom coupling definitions.

A user function is compiled from source with an AST whitelist, then executed
against a restricted namespace.  The contract is deliberately narrow:

    def couple(sc, fc, iu, iv, n_nodes):
        # sc, fc : (n_subjects, n_edges) float64
        # iu, iv : (n_edges,) endpoint indices, upper triangle, C order
        # returns (n_subjects, n_edges) | (n_edges,) | (n_nodes,) | scalar
        return sc * fc

Helpers exposed to the namespace: ``np``, ``zscore_edges``, ``mat_to_vec``,
``vec_to_mat``, ``upper_tri_index``.  No imports, no file access, no attribute
names starting with an underscore.
"""

from __future__ import annotations

import ast
import textwrap

import numpy as np

from .data import _zscore_by_group, mat_to_vec, upper_tri_index, vec_to_mat

CUSTOM_METRICS = {}


class CustomCodeError(ValueError):
    """Raised when a user function fails the static check."""


_ALLOWED_NODES = {
    "Module", "FunctionDef", "AsyncFunctionDef", "arguments", "arg", "Return",
    "Assign", "AugAssign", "AnnAssign", "Expr", "BinOp", "UnaryOp", "BoolOp",
    "Compare", "If", "IfExp", "For", "While", "Break", "Continue", "Pass",
    "Call", "Attribute", "Name", "Constant", "Subscript", "Slice", "Tuple",
    "List", "Dict", "Set", "Starred", "keyword", "Load", "Store", "Del",
    "Lambda", "ListComp", "SetComp", "DictComp", "GeneratorExp", "comprehension",
    "JoinedStr", "FormattedValue", "Import", "alias", "Index",
    # operators
    "Add", "Sub", "Mult", "MatMult", "Div", "FloorDiv", "Mod", "Pow",
    "LShift", "RShift", "BitOr", "BitXor", "BitAnd", "FloorDiv",
    "UAdd", "USub", "Invert", "Not", "And", "Or",
    "Eq", "NotEq", "Lt", "LtE", "Gt", "GtE", "Is", "IsNot", "In", "NotIn",
}

_FORBIDDEN_NAMES = {
    "eval", "exec", "open", "compile", "__import__", "getattr", "setattr",
    "delattr", "globals", "locals", "vars", "input", "breakpoint", "exit",
    "quit", "help", "dir", "memoryview", "object", "super", "classmethod",
    "staticmethod", "property", "type", "id", "hash", "repr", "format",
    "__builtins__", "__loader__", "__spec__", "__file__",
}

_ALLOWED_IMPORTS = {"numpy"}


def _check_node(node):
    kind = type(node).__name__
    if kind not in _ALLOWED_NODES:
        raise CustomCodeError("construct not allowed: %s" % kind)

    if kind == "Import":
        for alias in node.names:
            if alias.name not in _ALLOWED_IMPORTS:
                raise CustomCodeError("only 'import numpy' is allowed, got %r" % alias.name)
            if alias.asname not in (None, "np"):
                raise CustomCodeError("numpy must be imported as np")
    if kind in ("ImportFrom", "ClassDef", "Global", "Nonlocal", "With", "Delete",
                "Raise", "Assert", "Yield", "YieldFrom", "Await", "AsyncWith",
                "AsyncFor", "Try", "TryStar", "Match", "NamedExpr"):
        raise CustomCodeError("construct not allowed: %s" % kind)

    if kind == "Name":
        if node.id in _FORBIDDEN_NAMES:
            raise CustomCodeError("name not allowed: %s" % node.id)
        if node.id.startswith("__"):
            raise CustomCodeError("dunder name not allowed: %s" % node.id)

    if kind == "Attribute":
        if node.attr.startswith("_"):
            raise CustomCodeError("private attribute access not allowed: .%s" % node.attr)

    if kind == "Call":
        f = node.func
        if isinstance(f, ast.Name) and f.id in _FORBIDDEN_NAMES:
            raise CustomCodeError("call not allowed: %s" % f.id)
        if isinstance(f, ast.Attribute) and f.attr in ("system", "popen", "spawn",
                                                       "remove", "unlink", "rmtree",
                                                       "write", "save", "load",
                                                       "fromfile", "tofile"):
            raise CustomCodeError("call not allowed: .%s" % f.attr)

    for child in ast.iter_child_nodes(node):
        _check_node(child)


def validate_source(source, name="custom"):
    """Parse and statically check a user function; returns the compiled code."""
    source = textwrap.dedent(source or "").strip()
    if not source:
        raise CustomCodeError("empty source")
    try:
        tree = ast.parse(source, mode="exec")
    except SyntaxError as exc:
        raise CustomCodeError("syntax error: %s" % exc)
    _check_node(tree)
    fns = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]
    if not fns:
        raise CustomCodeError("source must define a function")
    for fn in fns:
        if isinstance(fn, ast.AsyncFunctionDef):
            raise CustomCodeError("async functions are not supported")
    chosen = None
    for fn in fns:
        if fn.name in ("couple", "coupling", "custom_coupling", name):
            chosen = fn
            break
    if chosen is None:
        if len(fns) == 1:
            chosen = fns[0]
        else:
            raise CustomCodeError(
                "define a function named 'couple' (found: %s)"
                % ", ".join(fn.name for fn in fns))
    return source, chosen.name


def _namespace():
    def _no_print(*_a, **_k):
        return None

    return {
        "__builtins__": {
            "abs": abs, "min": min, "max": max, "sum": sum, "len": len,
            "range": range, "enumerate": enumerate, "zip": zip,
            "float": float, "int": int, "bool": bool, "round": round,
            "sorted": sorted, "list": list, "dict": dict, "tuple": tuple,
            "print": _no_print, "True": True, "False": False, "None": None,
        },
        "np": np,
        "zscore_edges": lambda x, iu, n: _zscore_by_group(
            np.asarray(x, dtype=np.float64), np.asarray(iu), int(n)),
        "mat_to_vec": mat_to_vec,
        "vec_to_mat": vec_to_mat,
        "upper_tri_index": upper_tri_index,
    }


def compile_source(source, name="custom"):
    """Compile a checked source string into a callable ``f(sc, fc, iu, iv, n_nodes)``."""
    src, fname = validate_source(source, name=name)
    ns = _namespace()
    exec(compile(src, "<custom:%s>" % name, "exec"), ns)  # noqa: S102 - AST checked
    fn = ns[fname]
    if not callable(fn):
        raise CustomCodeError("%s is not callable" % fname)
    return fn


def _coerce(out, cohort):
    """Map a user return value onto the canonical (n_subjects, n_edges) field."""
    arr = np.asarray(out, dtype=np.float64)
    n, e, p = cohort.n, cohort.n_edges, cohort.n_nodes
    if arr.ndim == 0:
        return np.full((n, e), float(arr))
    if arr.shape == (n, e):
        return arr
    if arr.shape == (e,):
        return np.tile(arr[None, :], (n, 1))
    if arr.shape == (n, p):
        return ("regional", arr)
    if arr.shape == (p,):
        return ("regional", np.tile(arr[None, :], (n, 1)))
    if arr.shape == (n, p, p):
        return arr[:, cohort.iu, cohort.iv]
    if arr.shape == (p, p):
        return np.tile(arr[cohort.iu, cohort.iv][None, :], (n, 1))
    raise CustomCodeError(
        "return shape %s is not usable; expected (n_subjects, n_edges) = (%d, %d), "
        "(n_edges,), (n_nodes,) or (n_nodes, n_nodes)" % (arr.shape, n, e))


def evaluate(fn, cohort):
    """Run a compiled custom function against a cohort.

    Returns ``(kind, array)`` where kind is ``"edge"`` or ``"regional"``.
    """
    out = fn(cohort.sc.astype(np.float64), cohort.fc.astype(np.float64),
             cohort.iu, cohort.iv, cohort.n_nodes)
    res = _coerce(out, cohort)
    if isinstance(res, tuple):
        kind, arr = res
    else:
        kind, arr = "edge", res
    if not np.all(np.isfinite(arr)) and np.isnan(arr).all():
        raise CustomCodeError("custom function returned all-NaN values")
    return kind, arr


def register_custom(name, source, description="", family="custom"):
    """Validate, compile and register a custom definition under ``name``."""
    from .coords import register_metric

    name = str(name).strip()
    if not name:
        raise CustomCodeError("a name is required")
    if not name.replace("_", "").isalnum():
        raise CustomCodeError("name may contain letters, digits and underscores only")
    if name in ("direct", "reliance", "cosine", "matching", "communicability",
                "mfpt", "path_ease", "search_info"):
        raise CustomCodeError("%r is a built-in definition; pick another name" % name)
    fn = compile_source(source, name=name)
    CUSTOM_METRICS[name] = dict(fn=fn, source=source, description=description, family=family)
    register_metric(name, family=family, label=name, description=description, needs="edges")
    return name


def unregister_custom(name):
    from .coords import unregister_metric
    CUSTOM_METRICS.pop(name, None)
    unregister_metric(name)


def list_custom():
    return [dict(name=k, description=v["description"], family=v["family"])
            for k, v in CUSTOM_METRICS.items()]


def dry_run(source, n_subjects=8, n_nodes=24, seed=0):
    """Compile a source and run it on a small synthetic cohort."""
    from .demo import make_cohort
    fn = compile_source(source, name="custom")
    coh = make_cohort(n_subjects=n_subjects, n_nodes=n_nodes, seed=seed, behaviour=False)
    kind, arr = evaluate(fn, coh)
    return dict(kind=kind, shape=list(arr.shape),
                finite=bool(np.all(np.isfinite(arr))),
                mean=float(np.nanmean(arr)), sd=float(np.nanstd(arr)),
                min=float(np.nanmin(arr)), max=float(np.nanmax(arr)))
