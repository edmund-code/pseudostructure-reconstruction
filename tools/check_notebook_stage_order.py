#!/usr/bin/env python3
"""Static ordering check for the workflow notebooks.

Two failure modes this catches, both of which have broken a real run:

1. **use before definition** - a cell reads a name that no earlier cell binds. The notebooks share
   one namespace and are executed top to bottom, so a name that only appears later is a NameError.
2. **maybe-assigned** - a name that an earlier cell binds *only inside a conditional branch* (an
   ``if`` body, a ``try`` block, a loop) and a later cell then reads unconditionally. That is the
   subtler case: the notebook works when the branch was taken and fails otherwise, e.g. a cache flag
   initialised only on a cache hit, which crashes on the first run when the cache is empty.

Notebooks are canonical and the ``.py`` files are generated mirrors, so this reads the ``.ipynb``
files and ignores magics (which are not valid Python). Run it before committing a notebook change:

    python tools/check_notebook_stage_order.py analysis/notebooks/*.ipynb
"""
from __future__ import annotations

import argparse
import ast
import builtins
import json
import re
import sys
from pathlib import Path

BUILTIN_NAMES = set(dir(builtins)) | {'display', 'get_ipython', '__file__', 'Image', 'In', 'Out'}

# Stage-cache entry points take a zero-argument compute callback and store what it returns, so a
# callback with no `return` silently caches None. `cached_neighbor_graph` takes it third (adata,
# key, compute); the others take it second.
CACHE_CALLBACKS = {
    'cached_anndata': 1, 'cached_payload': 1, 'cached_frame': 1, 'cached_neighbor_graph': 2,
}
MAGIC = re.compile(r'^\s*[%!]')
CONTINUATION = re.compile(r'^\s*[%!]\s*=')          # a line starting with '!=' is not a magic


def _strip_magics(source: str) -> str:
    lines = []
    for line in source.splitlines():
        lines.append('' if (MAGIC.match(line) and not CONTINUATION.match(line)) else line)
    return '\n'.join(lines)


def _bindings(tree: ast.AST) -> tuple[set[str], set[str]]:
    """Return (certainly bound, maybe bound) top-level names for one cell.

    A name is *certain* when a top-level statement binds it, or when both branches of a top-level
    ``if``/``else`` bind it. It is *maybe* when only one conditional branch binds it.
    """
    certain: set[str] = set()
    maybe: set[str] = set()

    def bound_names(node: ast.AST) -> set[str]:
        names: set[str] = set()
        for child in ast.walk(node):
            if isinstance(child, ast.Name) and isinstance(child.ctx, (ast.Store, ast.Del)):
                names.add(child.id)
            elif isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                names.add(child.name)
            elif isinstance(child, ast.arg):
                names.add(child.arg)
            elif isinstance(child, (ast.Import, ast.ImportFrom)):
                for alias in child.names:
                    names.add((alias.asname or alias.name).split('.')[0])
            elif isinstance(child, ast.ExceptHandler) and child.name:
                names.add(child.name)
        return names

    for statement in tree.body:
        if isinstance(statement, ast.If):
            body = set().union(*(bound_names(part) for part in statement.body)) if statement.body else set()
            other = (set().union(*(bound_names(part) for part in statement.orelse))
                     if statement.orelse else set())
            certain |= body & other
            maybe |= (body | other) - (body & other)
        elif isinstance(statement, (ast.Try, ast.For, ast.While, ast.With)):
            # Bound when the construct's body runs; a later read of a loop target is the author's
            # business, and flagging every such name drowns the check in noise.
            certain |= bound_names(statement)
        else:
            certain |= bound_names(statement)
    return certain, maybe


def _returns(tree: ast.AST) -> bool:
    return any(isinstance(node, ast.Return) for node in ast.walk(tree))


def _compute_callback_problems(tree: ast.AST, index: int, known: dict[str, bool]) -> list[str]:
    problems = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Name):
            continue
        position = CACHE_CALLBACKS.get(node.func.id)
        if position is None:
            continue
        callback = node.args[position] if len(node.args) > position else None
        if callback is None:
            for keyword in node.keywords:
                if keyword.arg in ('compute', 'callback'):
                    callback = keyword.value
        if isinstance(callback, ast.Name) and known.get(callback.id) is False:
            problems.append(
                f'cell {index}: {callback.id!r} is passed as the compute callback to '
                f'{node.func.id} but has no return statement, so None would be cached')
    return problems


def check_notebook(path: Path) -> list[str]:
    notebook = json.loads(path.read_text())
    function_returns: dict[str, bool] = {}
    certain: set[str] = set()
    maybe: dict[str, int] = {}
    problems: list[str] = []
    for index, cell in enumerate(notebook['cells']):
        if cell['cell_type'] != 'code':
            continue
        try:
            tree = ast.parse(_strip_magics(''.join(cell['source'])))
        except SyntaxError as exc:
            problems.append(f'cell {index}: cannot parse ({exc.msg}, line {exc.lineno})')
            continue
        # A cell's own bindings count for its own reads (function bodies legitimately reference names
        # bound further down), so they are collected before the cell's uses are examined.
        cell_certain, cell_maybe = _bindings(tree)
        visible = certain | cell_certain | cell_maybe
        problems.extend(_compute_callback_problems(tree, index, function_returns))
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef):
                function_returns[node.name] = _returns(node)
        for statement in tree.body:
            for node in ast.walk(statement):
                if not isinstance(node, ast.Name) or not isinstance(node.ctx, ast.Load):
                    continue
                if node.id in BUILTIN_NAMES or node.id in visible:
                    continue
                if node.id in maybe:
                    origin = maybe.get(node.id, index)
                    problems.append(
                        f'cell {index}: {node.id!r} may be unbound - it is assigned only inside a '
                        f'conditional branch of cell {origin}')
                else:
                    problems.append(f'cell {index}: {node.id!r} is used before any cell defines it')
        certain |= cell_certain
        new_maybe = cell_maybe
        new_certain = cell_certain
        for name in new_maybe - new_certain:
            maybe.setdefault(name, index)
        for name in new_certain:
            maybe.pop(name, None)

    # de-duplicate: one report per (cell, name)
    seen: set[tuple[str, str]] = set()
    unique = []
    for problem in problems:
        key = problem.split(':', 2)[:2] if problem.count(':') >= 2 else (problem, '')
        key = (str(key), problem.split("'")[1] if "'" in problem else problem)
        if key in seen:
            continue
        seen.add(key)
        unique.append(problem)
    return unique


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('notebooks', nargs='+', type=Path)
    args = parser.parse_args(argv)
    failed = False
    for path in args.notebooks:
        problems = check_notebook(path)
        if problems:
            failed = True
            print(f'{path}: {len(problems)} ordering problem(s)')
            for problem in problems:
                print(f'  {problem}')
        else:
            print(f'{path}: name order clean')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
