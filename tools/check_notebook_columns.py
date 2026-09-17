"""Static audit of 06: every frame column read by literal name must exist in that frame."""
import ast
import json, re, sys

path = sys.argv[1]
nb = json.load(open(path))
cells = ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]
code = "\n".join(cells)

def cell_with(marker):
    hits = [text for text in cells if marker in text]
    assert len(hits) == 1, (marker, len(hits))
    return hits[0]

# --- gene_metrics column universe ------------------------------------------------------
templates = re.findall(r'f"([A-Za-z0-9_]+)\{suffix\}"', cell_with("def _curve_metrics"))
alt_suffix = ("_landmark_registered" if '"primary_coordinate": "shared_dpt"' in code
              else "_shared_dpt")
metric_columns = {name for template in templates
                  for name in (template, template + alt_suffix)}

alt_bases = re.search(r"ALT_COLUMNS = \[f\"\{name\}\{ALT_SUFFIX\}\" for name in \[(.*?)\]\]", code, re.S)
curated = {f"{name}{alt_suffix}" for name in re.findall(r'"([A-Za-z0-9_]+)"', alt_bases.group(1))}
not_carried = sorted(curated - metric_columns)
curated_bad = sorted(column for column in curated if column not in metric_columns)

# Parse the constructor rather than regexing the cell: the same cell also holds a mapping dict of
# 05's atlas names onto ours, whose keys are not columns of curve_table.
curve_cell = cell_with("curve_table = pd.DataFrame({")
curve_columns = set()
for node in ast.walk(ast.parse(curve_cell)):
    if (isinstance(node, ast.Assign) and isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Attribute) and node.value.func.attr == "DataFrame"
            and node.value.args and isinstance(node.value.args[0], ast.Dict)
            and isinstance(node.targets[0], ast.Name) and node.targets[0].id == "curve_table"):
        curve_columns |= {key.value for key in node.value.args[0].keys
                          if isinstance(key, ast.Constant) and isinstance(key.value, str)}
    if (isinstance(node, ast.Assign) and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Subscript)
            and isinstance(node.targets[0].value, ast.Name)
            and node.targets[0].value.id == "curve_table"
            and isinstance(node.targets[0].slice, ast.Constant)):
        curve_columns.add(node.targets[0].slice.value)
assert curve_columns, "curve_table columns not extracted"
print(f"curve_table columns parsed: {len(curve_columns)}")
print(f"curve_table columns parsed: {len(curve_columns)}")

flags_cell = cell_with("def _phenotype_flags")
flag_columns = set(re.findall(r'flags\["([A-Za-z0-9_]+)"\]', flags_cell))

assigned = set(re.findall(r'gene_metrics\["([A-Za-z0-9_]+)"\]\s*=', code))
merged = re.search(r'gene_metrics\.merge\(\s*robustness\[\[(.*?)\]\]', code, re.S)
merged_columns = set(re.findall(r'"([A-Za-z0-9_]+)"', merged.group(1))) if merged else set()

known = ({"gene", "spatial_phenotype", "broad_phenotype"} | set(templates) | curated
         | curve_columns | flag_columns | assigned | merged_columns)

# --- every read of a gene_metrics-derived frame -----------------------------------------
frames = ("gene_metrics", "variant_metrics", "quadrant_frame", "phase_frame", "landscape",
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

metric_reads = reads.get("gene_metrics", set()) | reads.get("variant_metrics", set()) \
    | reads.get("quadrant_frame", set()) | reads.get("phase_frame", set()) | reads.get("landscape", set())
missing = sorted(metric_reads - known -
                 {"conventional_abs", "conventional_percentile", "spatial_discovery_score",
                  "spatial_phenotype", "broad_phenotype"})

# A column present in both the metric frames and curve_table would be silently renamed to _x/_y by
# the merge, so the name read later would not exist. Check the overlap explicitly: this is the failure
# mode a union-of-all-names audit cannot see.
collisions = sorted((set(templates) | curated) & set(curve_columns))
print(f"merge collisions between the metric frames and curve_table: {collisions or 'none'}")

print(f"metric templates: {len(templates)} | curated unregistered columns: {len(curated)}")
print(f"gene_metrics column universe: {len(known)} | literal columns read: {len(metric_reads)}")
print("curated columns that do not exist in the metric frames:", curated_bad or "none")
print("UNRESOLVED READS:", missing or "none")
sys.exit(1 if (missing or curated_bad or collisions) else 0)
