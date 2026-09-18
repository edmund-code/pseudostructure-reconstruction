#!/usr/bin/env python3
"""Static audit of 06: every frame column read by literal name must exist in that frame.

Three ways a column can disappear without a syntax error, all checked here:

1. a column read from a metric frame (``gene_metrics["..."]``) that no cell ever defines;
2. a column written into the atlas (``atlas_columns``) that is not in the metric frame - a failure
   that would otherwise surface only at the end of a long run, after every figure had been drawn;
3. a column present in *both* the metric frame and ``curve_table``, which ``merge`` silently renames
   to ``_x``/``_y`` so the name read afterwards does not exist.

The extraction is guarded against its own vacuity: the metric universe must contain the columns the
notebook's own section-2 markdown lists, so a refactor that renames the metric table cannot quietly
turn this gate into a no-op.
"""
import ast
import json
import re
import sys

# Columns section 2 documents. If a refactor renames the metric table the extraction below returns
# nothing at all, so it is cross-checked against this list rather than trusted.
REQUIRED_METRIC_COLUMNS = {
    "level_effect_human_minus_mouse", "mouse_amplitude", "human_amplitude",
    "amplitude_log2_ratio_human_over_mouse", "shape_corr", "pattern_rms_z",
    "best_shift_human_minus_mouse", "post_shift_shape_corr", "shift_improvement",
    "residual_rms_after_shift", "shift_at_search_boundary", "mouse_position_centroid",
    "human_position_centroid", "centroid_shift_human_minus_mouse", "mouse_early_to_late",
    "human_early_to_late", "n_shared_grid_points", "enough_shared_support",
}

path = sys.argv[1]
nb = json.load(open(path))
cells = ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]
code = "\n".join(cells)


def cell_with(marker):
    hits = [text for text in cells if marker in text]
    assert len(hits) == 1, (marker, len(hits))
    return hits[0]


def _is_dataframe(value):
    """Is this call `pd.DataFrame(...)` / `DataFrame(...)`?"""
    function = getattr(value, "func", None)
    return ((isinstance(function, ast.Attribute) and function.attr == "DataFrame")
            or (isinstance(function, ast.Name) and function.id == "DataFrame"))


def literal_strings(cell_text, name, node_type):
    """The string elements of the dict/list literal assigned to `name` in `cell_text`."""
    found = set()
    for node in ast.walk(ast.parse(cell_text)):
        if not isinstance(node, ast.Assign) or not isinstance(node.value, node_type):
            continue
        if not any(isinstance(target, ast.Name) and target.id == name for target in node.targets):
            continue
        elements = node.value.keys if node_type is ast.Dict else node.value.elts
        found |= {element.value for element in elements
                  if isinstance(element, ast.Constant) and isinstance(element.value, str)}
    return found


def dataframe_dict_keys(cell_text, name):
    """The string keys of a `name = pd.DataFrame({...})` (or `name = {...}`) constructor."""
    keys = set()
    for node in ast.walk(ast.parse(cell_text)):
        if not isinstance(node, ast.Assign):
            continue
        if not any(isinstance(target, ast.Name) and target.id == name for target in node.targets):
            continue
        value = node.value
        if isinstance(value, ast.Call) and value.args and _is_dataframe(value):
            value = value.args[0]
        if isinstance(value, ast.Dict):
            keys |= {key.value for key in value.keys
                     if isinstance(key, ast.Constant) and isinstance(key.value, str)}
    return keys


# --- the metric frame: the columns the _curve_metrics constructor returns ----------------
metric_columns = set()
for node in ast.walk(ast.parse(cell_with("def _curve_metrics"))):
    if (isinstance(node, ast.Return) and isinstance(node.value, ast.Call)
            and _is_dataframe(node.value) and node.value.args
            and isinstance(node.value.args[0], ast.Dict)):
        metric_columns |= {key.value for key in node.value.args[0].keys
                           if isinstance(key, ast.Constant) and isinstance(key.value, str)}
missing_required = sorted(REQUIRED_METRIC_COLUMNS - metric_columns)
assert not missing_required, f"not extracted from _curve_metrics: {missing_required}"

# --- the atlas: the columns written to tables/gene_spatial_rewiring_atlas.csv -------------
atlas_columns = literal_strings(cell_with("atlas_columns = ["), "atlas_columns", ast.List)
assert len(atlas_columns) >= 20, sorted(atlas_columns)

# --- curve_table: parsed from its constructor rather than regexed ------------------------
# The same cell also holds a mapping dict of 05's atlas names onto ours, whose keys are not columns
# of curve_table, so the assignment target is checked.
curve_cell = cell_with("curve_table = pd.DataFrame({")
curve_columns = dataframe_dict_keys(curve_cell, "curve_table")
for node in ast.walk(ast.parse(curve_cell)):
    if (isinstance(node, ast.Assign) and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Subscript)
            and isinstance(node.targets[0].value, ast.Name)
            and node.targets[0].value.id == "curve_table"
            and isinstance(node.targets[0].slice, ast.Constant)):
        curve_columns.add(node.targets[0].slice.value)
assert curve_columns, "curve_table columns not extracted"

flags_cell = cell_with("def _phenotype_flags")
flag_columns = set(re.findall(r'flags\["([A-Za-z0-9_]+)"\]', flags_cell))
assert len(flag_columns) >= 15, sorted(flag_columns)

assigned = set(re.findall(r'gene_metrics\["([A-Za-z0-9_]+)"\]\s*=', code))
merged = re.search(r'gene_metrics\.merge\(\s*robustness\[\[(.*?)\]\]', code, re.S)
merged_columns = set(re.findall(r'"([A-Za-z0-9_]+)"', merged.group(1))) if merged else set()

known = ({"gene", "spatial_phenotype", "broad_phenotype"} | metric_columns | atlas_columns
         | curve_columns | flag_columns | assigned | merged_columns)

# --- every read of a gene_metrics-derived frame -----------------------------------------
frames = ("gene_metrics", "variant_metrics", "quadrant_frame", "displacement_frame", "landscape",
          "block", "robustness", "pathway_enrichment", "membership", "pathway_tested")
reads = {}
for text in cells:
    for frame in frames:
        # both quote styles: the closing summary uses single quotes inside f-strings
        found = set(re.findall(rf"{frame}\[\"([A-Za-z0-9_]+)\"\]", text))
        found |= set(re.findall(rf"{frame}\['([A-Za-z0-9_]+)'\]", text))
        if found:
            reads.setdefault(frame, set()).update(found)

# Columns assigned onto any of these frames count as defined, not only gene_metrics assignments.
frame_assigned = set()
for frame in frames:
    frame_assigned |= set(re.findall(rf"{frame}\[\"([A-Za-z0-9_]+)\"\]\s*=", code))
    frame_assigned |= set(re.findall(rf"{frame}\['([A-Za-z0-9_]+)'\]\s*=", code))
known |= frame_assigned

metric_reads = (reads.get("gene_metrics", set()) | reads.get("variant_metrics", set())
                | reads.get("quadrant_frame", set()) | reads.get("displacement_frame", set())
                | reads.get("landscape", set()))
missing = sorted(metric_reads - known -
                 {"conventional_abs", "conventional_percentile", "spatial_discovery_score",
                  "spatial_phenotype", "broad_phenotype"})

# A column present in both the metric frames and curve_table would be silently renamed to _x/_y by
# the merge, so the name read later would not exist. Check the overlap explicitly: this is the failure
# mode a union-of-all-names audit cannot see.
collisions = sorted(metric_columns & curve_columns)
atlas_bad = sorted(atlas_columns - known)

print(f"merge collisions between the metric frames and curve_table: {collisions or 'none'}")
print(f"metric columns: {len(metric_columns)} | atlas columns: {len(atlas_columns)} "
      f"| curve_table columns: {len(curve_columns)} | flags: {len(flag_columns)}")
print(f"gene_metrics column universe: {len(known)} | literal columns read: {len(metric_reads)}")
print("atlas columns that do not exist in the metric frames:", atlas_bad or "none")
print("UNRESOLVED READS:", missing or "none")
sys.exit(1 if (missing or atlas_bad or collisions) else 0)
