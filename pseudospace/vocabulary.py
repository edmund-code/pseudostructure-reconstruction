"""The single nephron segment vocabulary shared by every analysis workflow.

The mouse-only and the cross-species workflows cluster the same tissue and label the same kind of
Leiden cluster, so they have to agree on what a label means. This module is the one definition of
that vocabulary; the notebooks import it instead of each keeping its own copy.

Two levels, as in the mouse-only workflow:

* **fine segments** -- the segments the marker panel can distinguish, in anatomical order;
* **coarse families** -- the nephron continuum ``PT, DTL, AL, DCT, CNT_CD``.

A coarse family is also a valid label in its own right, not merely a rollup target: at the
resolutions these pipelines cluster at, most clusters carry coarse identity only, and forcing a
fine label onto such a cluster invents a distinction the data does not support. Label each cluster
at the granularity its markers support, and no finer.

``AL`` is the ascending limb family: ``ATL``, ``mTAL``, ``cTAL`` and ``Macula-densa`` are not
separated at this resolution, so they share one family. A cluster whose markers do resolve the
thick limb can be labelled ``mTAL``/``cTAL`` instead -- a fine label rolls up to ``AL`` either way.
"""

from __future__ import annotations

# --- fine segments, in anatomical order: glomerulus -> PT -> thin limb -> ascending limb -> DCT -> CNT/CD
FINE_SEGMENTS: tuple[str, ...] = (
    "Podocyte",
    "PT-S1", "PT-S2", "PT-S3",
    "DTL1", "DTL2", "DTL3",
    "ATL", "mTAL", "cTAL", "Macula-densa",
    "DCT1", "DCT2",
    "CNT", "CCD", "OMCD", "IMCD",
)

# --- the nephron continuum, in anatomical order (the reported coarse classes)
COARSE_FAMILIES: tuple[str, ...] = ("PT", "DTL", "AL", "DCT", "CNT_CD")

# --- labels outside the tubular continuum; removed before the pass-2 nephron analysis
REMOVE_CLASSES: tuple[str, ...] = (
    "Glomerulus", "Vessel", "Stroma", "SmoothMuscle", "Immune", "Unassigned",
)

# --- fine segment -> coarse family
SEGMENT_TO_COARSE: dict[str, str] = {
    "Podocyte": "Glomerulus",
    "PT-S1": "PT", "PT-S2": "PT", "PT-S3": "PT",
    "DTL1": "DTL", "DTL2": "DTL", "DTL3": "DTL",
    "ATL": "AL", "mTAL": "AL", "cTAL": "AL", "Macula-densa": "AL",
    "DCT1": "DCT", "DCT2": "DCT",
    "CNT": "CNT_CD", "CCD": "CNT_CD", "OMCD": "CNT_CD", "IMCD": "CNT_CD",
}
# A family is a label in its own right, so each maps to itself, as does Glomerulus.
SEGMENT_TO_COARSE.update({family: family for family in COARSE_FAMILIES})
SEGMENT_TO_COARSE["Glomerulus"] = "Glomerulus"

# --- every label that may be typed at a manual cluster-review checkpoint
LABEL_VOCABULARY: tuple[str, ...] = (
    FINE_SEGMENTS
    + COARSE_FAMILIES
    + ("Glomerulus",)
    + tuple(c for c in REMOVE_CLASSES if c != "Glomerulus")
)

# --- labels that survive into the pass-2 tubular-nephron analysis
KEEP_TUBULE_CLASSES: tuple[str, ...] = COARSE_FAMILIES


def _build_display_order() -> tuple[str, ...]:
    """Each coarse family immediately before its own first fine member.

    ``segment_class`` may hold a mix of fine and coarse labels, so ordering by families keeps the
    anatomical left-to-right reading intact: ... PT, PT-S1, PT-S2, PT-S3, DTL, DTL1, ...
    """
    order: list[str] = []
    seen: set[str] = set()
    for segment in FINE_SEGMENTS:
        family = SEGMENT_TO_COARSE.get(segment)
        if family in COARSE_FAMILIES and family not in seen:
            order.append(family)
            seen.add(family)
        order.append(segment)
    for extra in ("Glomerulus",) + tuple(c for c in REMOVE_CLASSES if c != "Glomerulus"):
        if extra not in order:
            order.append(extra)
    return tuple(order)


# --- display order for the fine `segment_class` column
SEGMENT_DISPLAY_ORDER: tuple[str, ...] = _build_display_order()


def coarse_for(label: str) -> str:
    """Coarse family for any label. Unmapped labels become ``Unassigned`` rather than raising."""
    return SEGMENT_TO_COARSE.get(label, "Unassigned")


def invalid_labels(labels) -> list[str]:
    """Labels that cannot be typed at a checkpoint, sorted. Empty means the set is legal."""
    return sorted(set(map(str, labels)) - set(LABEL_VOCABULARY))


def unlabelled_clusters(cluster_ids, labels) -> list[str]:
    """Cluster IDs with no entry in the label map, sorted. Empty means every cluster is reviewed."""
    return sorted(set(map(str, cluster_ids)) - set(map(str, labels)))


def stale_cluster_ids(cluster_ids, labels) -> list[str]:
    """Label-map keys that are absent from this run, sorted. Empty means no stale entries."""
    return sorted(set(map(str, labels)) - set(map(str, cluster_ids)))


def assert_consistent() -> None:
    """Guard the definitions against each other. Raises ``ValueError`` on any inconsistency."""
    if len(set(FINE_SEGMENTS)) != len(FINE_SEGMENTS):
        raise ValueError("FINE_SEGMENTS contains duplicates")
    if len(set(COARSE_FAMILIES)) != len(COARSE_FAMILIES):
        raise ValueError("COARSE_FAMILIES contains duplicates")
    unmapped = [s for s in FINE_SEGMENTS if s not in SEGMENT_TO_COARSE]
    if unmapped:
        raise ValueError(f"fine segments with no coarse family: {unmapped}")
    to_unknown = sorted({v for v in SEGMENT_TO_COARSE.values() if v not in COARSE_FAMILIES + ("Glomerulus",)})
    if to_unknown:
        raise ValueError(f"rollup targets outside COARSE_FAMILIES/Glomerulus: {to_unknown}")
    duplicate_labels = sorted({x for x in LABEL_VOCABULARY if LABEL_VOCABULARY.count(x) > 1})
    if duplicate_labels:
        raise ValueError(f"labels listed twice in LABEL_VOCABULARY: {duplicate_labels}")
    missing_from_display = sorted(set(LABEL_VOCABULARY) - set(SEGMENT_DISPLAY_ORDER))
    if missing_from_display:
        raise ValueError(f"labels missing from SEGMENT_DISPLAY_ORDER: {missing_from_display}")


assert_consistent()
