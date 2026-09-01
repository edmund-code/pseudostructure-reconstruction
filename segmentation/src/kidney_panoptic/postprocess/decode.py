"""Panoptic instance decode: local fields -> instance ids -> per-class cleanup.

The decode reads only LOCAL fields (semantic probability, interior/boundary probability,
center heatmaps). Nothing here depends on a non-local quantity such as an object's centroid
offset, which is what lets the same function run identically on a 512 px patch and on a
2048 px window of a whole slide, and is why a tiling seam cannot manufacture a split
(docs/adr/0001, docs/adr/0005).

Pipeline
--------
1. foreground   = 1 - p_background  (from the boundary head), gated by tissue if supplied
2. seeds        = connected components of p_interior above threshold, specks dropped
3. instances    = marker-controlled watershed flooding the BOUNDARY probability inside the
                  foreground mask — basins meet on the ridge between touching objects
4. class vote   = majority semantic argmax inside each instance, boundary pixels excluded
                  (they are half-and-half by construction and would bias the vote)
5. cleanup      = class-heterogeneous, because a glomerulus and a tubule fail differently:
                  compact classes get center-peak splitting; elongated classes get fragment
                  rejection and an optional, deliberately conservative re-merge; every class
                  gets its own minimum area
6. tissue gate  = drop instances lying mostly off tissue

Note on NMS: a label map assigns each pixel exactly one id, so decoded instances are
disjoint by construction and IoU between any two of them is zero. Overlap deduplication is
therefore meaningless HERE and lives where overlapping candidates genuinely exist — the
sliding-window stitcher (postprocess/stitch.py), which matches instances across window
seams at ``infer.stitch_iou``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.measure import regionprops
from skimage.segmentation import watershed

from ..data.classes import ClassSpec


@dataclass
class InstanceResult:
    id_map: np.ndarray                       # (H, W) int32, 0 = background
    class_of_id: dict[int, int] = field(default_factory=dict)   # id -> class index
    score_of_id: dict[int, float] = field(default_factory=dict)

    @property
    def n_instances(self) -> int:
        return len(self.class_of_id)

    def class_name_of_id(self, spec: ClassSpec) -> dict[int, str]:
        return {i: spec.names[c] for i, c in self.class_of_id.items()}


# ---------------------------------------------------------------------- decode
def decode_panoptic(
    semantic_prob: np.ndarray,      # (C, H, W) float32
    boundary_prob: np.ndarray,      # (3|4, H, W) — bg / interior / boundary [/ touching]
    center_heatmaps: Optional[np.ndarray],   # (K, H, W) float32 or None
    cfg: dict,
    spec: ClassSpec,
    tissue_mask: Optional[np.ndarray] = None,
) -> InstanceResult:
    dcfg = dict(cfg.get("decode", {}) or {})
    # A 4-channel boundary head splits the wall into object-to-background (2) and
    # object-to-INSTANCE (3). Both are ridges for the flooding, so the elevation sums them.
    p_bg, p_int = boundary_prob[0], boundary_prob[1]
    has_touch = boundary_prob.shape[0] > 3
    p_touch = boundary_prob[3] if has_touch else None
    p_bnd = boundary_prob[2:].sum(axis=0) if has_touch else boundary_prob[2]

    fg = (1.0 - p_bg) >= float(dcfg.get("fg_threshold", 0.5))
    if tissue_mask is not None and cfg.get("infer", {}).get("tissue_gate", True):
        fg &= tissue_mask.astype(bool)
    if not fg.any():
        return InstanceResult(np.zeros(p_bg.shape, np.int32))

    seeds = (p_int >= float(dcfg.get("interior_threshold", 0.5))) & fg

    # Carve the instance-instance wall OUT of the seed region.
    #
    # This is where the touching channel has to act. Merging is a SEEDING failure: when two
    # objects' interiors stay connected they yield a single marker, and a marker-controlled
    # watershed can never split one marker no matter how high the ridge between them is
    # raised. Vetoing seeds on predicted wall pixels disconnects the two interiors, which
    # turns one marker into two and lets the watershed separate them.
    #
    # The threshold is low on purpose. Measured on validation patches, mean p(touching) is
    # 0.142 on true instance-instance walls versus 0.014 on object-to-background contours
    # and 0.0014 in interiors — 10x discriminative but poorly calibrated, because the class
    # is a tiny fraction of pixels. An absolute cut near 0.05 sits in that gap; reading it
    # as a calibrated probability and using 0.5 would make the veto never fire.
    veto = float(dcfg.get("touching_seed_veto", 0.0))
    if has_touch and veto > 0:
        seeds &= p_touch < veto
    markers, n = ndi.label(seeds)
    if n == 0:
        return InstanceResult(np.zeros(p_bg.shape, np.int32))
    markers = _drop_small_labels(markers, int(dcfg.get("min_seed_area", 24)))
    if markers.max() == 0:
        return InstanceResult(np.zeros(p_bg.shape, np.int32))

    labels = watershed(
        p_bnd.astype(np.float32), markers=markers, mask=fg,
        watershed_line=bool(dcfg.get("watershed_line", False)),
    ).astype(np.int32)

    # Split BEFORE class assignment, so each resulting fragment gets its own honest
    # majority vote and score from `assign_classes` rather than inheriting the parent's.
    if bool(dcfg.get("semantic_split", {}).get("enabled", False)):
        labels = split_by_semantic(labels, semantic_prob, p_int, p_bnd, spec, dcfg)

    class_of_id, score_of_id = assign_classes(labels, semantic_prob, p_bnd, dcfg)

    labels, class_of_id, score_of_id = class_specific_cleanup(
        labels, class_of_id, score_of_id, semantic_prob, boundary_prob,
        center_heatmaps, cfg, spec,
    )

    if tissue_mask is not None and cfg.get("infer", {}).get("tissue_gate", True):
        labels, class_of_id, score_of_id = drop_off_tissue(
            labels, class_of_id, score_of_id, tissue_mask,
            float(cfg.get("infer", {}).get("tissue_min_instance_frac", 0.5)),
        )

    labels, class_of_id, score_of_id = relabel(labels, class_of_id, score_of_id)
    return InstanceResult(labels, class_of_id, score_of_id)


# ------------------------------------------------------------------ class vote
def assign_classes(
    labels: np.ndarray,
    semantic_prob: np.ndarray,
    p_bnd: np.ndarray,
    dcfg: dict,
) -> tuple[dict[int, int], dict[int, float]]:
    """Majority semantic vote per instance, boundary pixels excluded.

    Background never wins: an instance exists because the instance head said so, so the
    vote is taken over thing-classes only. The reported score is the mean probability of
    the winning class inside the instance.
    """
    C = semantic_prob.shape[0]
    sem_arg = np.argmax(semantic_prob, axis=0)
    core = p_bnd < 0.5   # exclude the boundary band from the vote

    class_of_id: dict[int, int] = {}
    score_of_id: dict[int, float] = {}
    objects = ndi.find_objects(labels)
    for idx, sl in enumerate(objects, start=1):
        if sl is None:
            continue
        sub = labels[sl] == idx
        if not sub.any():
            continue
        vote_mask = sub & core[sl]
        if not vote_mask.any():
            vote_mask = sub
        votes = np.bincount(sem_arg[sl][vote_mask], minlength=C)
        votes[0] = 0                       # background is not a thing class
        if votes.max() == 0:
            # Every pixel argmaxed to background. Fall back to the strongest thing-class
            # by mean probability rather than dropping a detection the instance head made.
            means = [semantic_prob[c][sl][sub].mean() if c > 0 else -1.0 for c in range(C)]
            cls = int(np.argmax(means))
        else:
            cls = int(np.argmax(votes))
        class_of_id[idx] = cls
        score_of_id[idx] = float(semantic_prob[cls][sl][sub].mean())
    return class_of_id, score_of_id


# ----------------------------------------------------------------- semantic split
def split_by_semantic(
    labels: np.ndarray,          # (H, W) int32 watershed output
    semantic_prob: np.ndarray,   # (C, H, W) float32
    p_int: np.ndarray,           # (H, W) float32 — the confident-interior probability
    p_bnd: np.ndarray,           # (H, W) float32 — boundary ridge
    spec: ClassSpec,
    dcfg: dict,
) -> np.ndarray:
    """Split an instance whose confident interior holds >= 2 sizeable regions of
    DIFFERENT semantic classes.

    This is a second separation signal, independent of the boundary head: two touching
    tubules of different anatomical segments are cut even where the boundary head missed
    the wall. It CANNOT help same-class neighbours, which are the majority of contacts —
    that is the boundary head's job (docs/adr/0006).

    Returns only a label map. ``assign_classes`` runs afterwards and rebuilds the class and
    score dicts over the finer partition, so nothing is bookkept here.

    Over-splitting is the risk this function creates, so it is guarded five ways: a
    confidence floor, restriction to the confident interior, morphological opening, an
    absolute minimum component area, and a minority-share floor. All operate inside one
    instance's bounding box, so the operation stays local and translation-equivariant.
    """
    scfg = dict(dcfg.get("semantic_split", {}) or {})
    min_prob = float(scfg.get("min_class_prob", 0.50))
    open_px = int(scfg.get("open_px", 2))
    min_comp = int(scfg.get("min_component_area", 600))
    min_frac = float(scfg.get("min_minor_frac", 0.15))
    max_parts = int(scfg.get("max_parts", 4))
    mode = str(scfg.get("elevation", "boundary_plus_ambiguity"))
    amb_w = float(scfg.get("ambiguity_weight", 1.0))

    allowed = _pair_set(scfg.get("pairs"), spec)
    core_all = p_int >= float(dcfg.get("interior_threshold", 0.5))
    eligible = tuple(range(1, semantic_prob.shape[0]))   # background never seeds

    out = labels.copy()
    next_id = int(labels.max()) + 1
    for idx, sl in enumerate(ndi.find_objects(labels), start=1):
        if sl is None:
            continue
        sub = labels[sl] == idx
        area = int(sub.sum())
        if area == 0:
            continue
        sem_sub = semantic_prob[(slice(None),) + sl]

        markers, area_of_class = _semantic_markers(
            sub, sub & core_all[sl], sem_sub, eligible, min_prob, open_px, min_comp
        )
        if len(area_of_class) < 2:
            continue
        ranked = sorted(area_of_class.items(), key=lambda kv: -kv[1])
        if ranked[1][1] < min_frac * area:
            continue
        if allowed is not None and tuple(sorted((ranked[0][0], ranked[1][0]))) not in allowed:
            continue
        n_parts = int(markers.max())
        if n_parts < 2 or n_parts > max_parts:
            continue

        elev = _split_elevation(mode, p_bnd[sl], sem_sub, [c for c, _ in ranked], amb_w)
        parts = watershed(elev.astype(np.float32), markers=markers, mask=sub)

        # Every pixel of the parent must land in exactly one child. `watershed` with
        # mask=sub guarantees that; anything it left at 0 keeps the parent id so the split
        # can never punch a hole in the instance.
        region = out[sl]
        for j in range(2, n_parts + 1):
            m = (parts == j) & sub
            if m.any():
                region[m] = next_id
                next_id += 1
        out[sl] = region
    return out


def _semantic_markers(sub, core, sem_sub, eligible, min_prob, open_px, min_comp):
    """Cleaned per-class seed components inside one instance.

    Returns ``(marker map 1..M, {class index: surviving area})``. Opening is the speckle
    guard and matters most: without it, a handful of stray argmax pixels of a neighbouring
    class would seed a spurious basin.
    """
    marker = np.zeros(sub.shape, np.int32)
    area_of_class: dict[int, int] = {}
    m = 0
    fp = _disk(open_px)
    for c in eligible:
        cm = core & (sem_sub[c] >= min_prob)
        if not cm.any():
            continue
        if open_px > 0:
            cm = ndi.binary_opening(cm, structure=fp)
            if not cm.any():
                continue
        lab, n = ndi.label(cm)
        if n == 0:
            continue
        total = 0
        for k in range(1, n + 1):
            comp = lab == k
            a = int(comp.sum())
            if a < min_comp:
                continue
            m += 1
            marker[comp] = m
            total += a
        if total:
            area_of_class[c] = total
    return marker, area_of_class


def _split_elevation(mode: str, p_bnd_sub, sem_sub, candidates, amb_w: float):
    """Surface the split watershed floods.

    ``boundary`` alone is wrong here by construction: if the boundary head had fired, the
    blob would not have merged in the first place, so p_bnd is flat exactly where the cut
    is needed. The ambiguity term ``1 - max(class prob)`` peaks where the two classes'
    probabilities cross, which is where the segment border actually is. Summing them uses
    real boundary evidence when it exists and falls back on the class crossover when it
    does not.
    """
    if mode == "boundary":
        return p_bnd_sub
    amb = 1.0 - np.max(sem_sub[list(candidates)], axis=0)
    if mode == "ambiguity":
        return amb
    if mode == "boundary_plus_ambiguity":
        return p_bnd_sub + amb_w * amb
    raise ValueError(
        f"decode.semantic_split.elevation must be boundary | ambiguity | "
        f"boundary_plus_ambiguity, got '{mode}'"
    )


def _pair_set(pairs, spec: ClassSpec):
    """Config class-name pairs -> a set of sorted index tuples. None means 'any pair'."""
    if not pairs:
        return None
    out = set()
    for a, b in pairs:
        if a not in spec.names or b not in spec.names:
            raise ValueError(
                f"decode.semantic_split.pairs references unknown class in ({a}, {b}). "
                f"Known: {list(spec.names)}"
            )
        out.add(tuple(sorted((spec.index(a), spec.index(b)))))
    return out


def _disk(radius: int) -> np.ndarray:
    """Boolean disk footprint (radius 0 -> single pixel)."""
    r = int(radius)
    if r <= 0:
        return np.ones((1, 1), dtype=bool)
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return (x * x + y * y) <= r * r


# --------------------------------------------------------------------- cleanup
def class_specific_cleanup(
    labels: np.ndarray,
    class_of_id: dict[int, int],
    score_of_id: dict[int, float],
    semantic_prob: np.ndarray,
    boundary_prob: np.ndarray,
    center_heatmaps: Optional[np.ndarray],
    cfg: dict,
    spec: ClassSpec,
):
    dcfg = dict(cfg.get("decode", {}) or {})

    if center_heatmaps is not None and len(spec.compact_indices):
        labels, class_of_id, score_of_id = split_compact_by_centers(
            labels, class_of_id, score_of_id, center_heatmaps, spec, dcfg
        )

    if bool(dcfg.get("merge_fragments", False)):
        labels, class_of_id, score_of_id = merge_elongated_fragments(
            labels, class_of_id, score_of_id, boundary_prob[2], spec, dcfg
        )

    if bool(dcfg.get("fill_holes", False)):
        labels = _fill_instance_holes(labels)

    return filter_min_area(labels, class_of_id, score_of_id, spec, dcfg)


def split_compact_by_centers(
    labels: np.ndarray,
    class_of_id: dict[int, int],
    score_of_id: dict[int, float],
    center_heatmaps: np.ndarray,
    spec: ClassSpec,
    dcfg: dict,
):
    """Split a compact-class instance that contains several center peaks.

    Two adjacent glomeruli often have no boundary evidence strong enough to separate
    them — they are separated in the image by a thin capsule the boundary head can miss —
    but they do produce two distinct center peaks. Splitting is geometric (flood the
    negative distance transform from the peaks) rather than boundary-driven, which is the
    right prior for round objects.

    An instance with NO peak is kept, not dropped: the center head is an auxiliary signal
    and treating its silence as a veto would cost recall for nothing.
    """
    chan_of_class = {c: k for k, c in enumerate(spec.compact_indices)}
    thr = float(dcfg.get("center_threshold", 0.35))
    min_dist = int(dcfg.get("center_min_distance", 20))
    out = labels.copy()
    classes = dict(class_of_id)
    scores = dict(score_of_id)
    next_id = int(labels.max()) + 1

    objects = ndi.find_objects(labels)
    for idx, sl in enumerate(objects, start=1):
        if sl is None or idx not in classes:
            continue
        k = chan_of_class.get(classes[idx])
        if k is None:
            continue
        sub = labels[sl] == idx
        heat = center_heatmaps[k][sl] * sub
        if heat.max() < thr:
            continue
        peaks = peak_local_max(heat, min_distance=min_dist, threshold_abs=thr, labels=sub)
        if len(peaks) < 2:
            continue

        seeds = np.zeros(sub.shape, np.int32)
        for j, (py, px) in enumerate(peaks, start=1):
            seeds[py, px] = j
        dist = ndi.distance_transform_edt(sub)
        parts = watershed(-dist, markers=seeds, mask=sub)

        region = out[sl]
        for j in range(1, len(peaks) + 1):
            m = (parts == j) & sub
            if not m.any():
                continue
            if j == 1:
                region[m] = idx
            else:
                region[m] = next_id
                classes[next_id] = classes[idx]
                scores[next_id] = scores.get(idx, 0.0)
                next_id += 1
        out[sl] = region
    return out, classes, scores


def merge_elongated_fragments(
    labels: np.ndarray,
    class_of_id: dict[int, int],
    score_of_id: dict[int, float],
    p_bnd: np.ndarray,
    spec: ClassSpec,
    dcfg: dict,
):
    """Conservative re-merge of adjacent same-class elongated fragments.

    BOTH conditions must hold:
      * the shared interface carries little boundary probability (the split had no
        evidence behind it), and
      * the union's minor axis barely grows.

    The second is a hard reject for "parallel glue": two tubules lying side by side also
    have a weak interface, but fusing them roughly doubles the minor axis. Without that
    test an interface-only rule fuses genuinely distinct neighbours, which is a worse
    error than the split it repairs — hence merging is OFF by default.
    """
    elong = {spec.index(c) for c in spec.elongated if c in spec.names}
    max_bnd = float(dcfg.get("merge_max_boundary_prob", 0.30))
    max_growth = float(dcfg.get("merge_max_minor_axis_growth", 1.15))
    max_ratio = float(dcfg.get("merge_max_area_ratio", 4.0))

    props = {p.label: p for p in regionprops(labels)}
    pairs = _adjacent_pairs(labels)
    parent = {i: i for i in class_of_id}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for (a, b), interface in pairs.items():
        if a not in class_of_id or b not in class_of_id:
            continue
        if class_of_id[a] != class_of_id[b] or class_of_id[a] not in elong:
            continue
        if len(interface[0]) < 5:
            continue
        if float(p_bnd[interface].mean()) > max_bnd:
            continue
        pa, pb = props.get(a), props.get(b)
        if pa is None or pb is None:
            continue
        if max(pa.area, pb.area) / max(min(pa.area, pb.area), 1) > max_ratio:
            continue
        union = (labels == a) | (labels == b)
        minor_u = _minor_axis(union)
        if minor_u > max_growth * max(pa.minor_axis_length, pb.minor_axis_length, 1e-6):
            continue     # parallel glue — reject
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    out = labels.copy()
    classes, scores = {}, {}
    for i in class_of_id:
        root = find(i)
        if root != i:
            out[labels == i] = root
    for i in np.unique(out):
        if i == 0:
            continue
        classes[int(i)] = class_of_id[int(i)]
        scores[int(i)] = score_of_id.get(int(i), 0.0)
    return out, classes, scores


def filter_min_area(labels, class_of_id, score_of_id, spec: ClassSpec, dcfg: dict):
    """Drop instances below their class's minimum area."""
    min_area = dict(dcfg.get("min_area", {}) or {})
    default = int(min_area.get("default", 0))
    out = labels.copy()
    keep_c, keep_s = {}, {}
    objects = ndi.find_objects(labels)
    for idx, sl in enumerate(objects, start=1):
        if sl is None or idx not in class_of_id:
            continue
        area = int((labels[sl] == idx).sum())
        name = spec.names[class_of_id[idx]]
        if area < int(min_area.get(name, default)):
            sub = out[sl]
            sub[sub == idx] = 0
            out[sl] = sub
            continue
        keep_c[idx] = class_of_id[idx]
        keep_s[idx] = score_of_id.get(idx, 0.0)
    return out, keep_c, keep_s


def drop_off_tissue(labels, class_of_id, score_of_id, tissue_mask, min_frac: float):
    out = labels.copy()
    keep_c, keep_s = {}, {}
    tissue = tissue_mask.astype(bool)
    objects = ndi.find_objects(labels)
    for idx, sl in enumerate(objects, start=1):
        if sl is None or idx not in class_of_id:
            continue
        m = labels[sl] == idx
        if not m.any():
            continue
        if (tissue[sl] & m).sum() / m.sum() < min_frac:
            sub = out[sl]
            sub[sub == idx] = 0
            out[sl] = sub
            continue
        keep_c[idx] = class_of_id[idx]
        keep_s[idx] = score_of_id.get(idx, 0.0)
    return out, keep_c, keep_s


# ----------------------------------------------------------------- primitives
def relabel(labels: np.ndarray, class_of_id: dict, score_of_id: dict):
    """Compact ids to 1..N, carrying the class/score dicts along.

    A label with no entry in ``class_of_id`` is DROPPED, not defaulted to class 0 — an
    instance silently labelled "background" would be exported as a background polygon and
    counted as a detection, which is worse than not emitting it at all.
    """
    present = [int(i) for i in np.unique(labels) if i != 0 and int(i) in class_of_id]
    out = np.zeros_like(labels, dtype=np.int32)
    classes, scores = {}, {}
    for new_id, old in enumerate(sorted(present), start=1):
        out[labels == old] = new_id
        classes[new_id] = int(class_of_id[old])
        scores[new_id] = float(score_of_id.get(old, 0.0))
    return out, classes, scores


def _drop_small_labels(labels: np.ndarray, min_area: int) -> np.ndarray:
    if min_area <= 0:
        return labels
    counts = np.bincount(labels.ravel())
    small = np.where(counts < min_area)[0]
    if len(small) == 0:
        return labels
    lut = np.arange(len(counts), dtype=np.int32)
    lut[small] = 0
    lut[0] = 0
    out = lut[labels]
    return ndi.label(out > 0)[0].astype(np.int32)


def _minor_axis(mask: np.ndarray) -> float:
    props = regionprops(mask.astype(np.int32))
    return float(props[0].minor_axis_length) if props else 0.0


def _adjacent_pairs(labels: np.ndarray) -> dict[tuple[int, int], tuple]:
    """Map each adjacent label pair to the index tuple of their shared interface pixels."""
    pairs: dict[tuple[int, int], list] = {}
    for axis in (0, 1):
        a = labels
        b = np.roll(labels, -1, axis=axis)
        m = (a > 0) & (b > 0) & (a != b)
        if axis == 0:
            m[-1, :] = False
        else:
            m[:, -1] = False
        ys, xs = np.nonzero(m)
        for y, x in zip(ys, xs):
            key = (int(min(a[y, x], b[y, x])), int(max(a[y, x], b[y, x])))
            pairs.setdefault(key, []).append((y, x))
    return {k: (np.array([p[0] for p in v]), np.array([p[1] for p in v]))
            for k, v in pairs.items()}


def _fill_instance_holes(labels: np.ndarray) -> np.ndarray:
    out = labels.copy()
    for idx, sl in enumerate(ndi.find_objects(labels), start=1):
        if sl is None:
            continue
        m = labels[sl] == idx
        filled = ndi.binary_fill_holes(m)
        new = filled & (out[sl] == 0)
        sub = out[sl]
        sub[new] = idx
        out[sl] = sub
    return out
