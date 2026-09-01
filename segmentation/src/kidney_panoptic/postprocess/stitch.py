"""Slide-wide instance assembly from overlapping decode windows.

Instances are decoded on windows of ``infer.decode_window`` stepping by
``decode_window - decode_overlap``. An object straddling a window edge is decoded twice,
once per window, and the two copies must become one global instance.

Matching rule
-------------
When a window writes into territory an earlier window already covered, each local instance
is compared with whatever global instance already occupies those pixels. The two are joined
when

    |local ∩ global| / min(|local ∩ covered|, |global ∩ window|)  >=  infer.stitch_iou

Overlap over the SMALLER of the two is used rather than plain IoU because the two copies
are usually very unequal — one window may hold 10% of a long tubule and the other 100% —
and plain IoU would refuse the obviously correct match at 0.1.

Joins are recorded in a union-find and resolved once at the end, so a chain of three
windows crossing one vessel collapses to a single id without any ordering dependence.
Pixels already written are never overwritten; a later window only fills what is still empty.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
from scipy import ndimage as ndi


class UnionFind:
    def __init__(self):
        self.parent: dict[int, int] = {}

    def add(self, x: int) -> None:
        self.parent.setdefault(x, x)

    def find(self, x: int) -> int:
        self.add(x)
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:      # path compression
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a: int, b: int) -> int:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return ra
        lo, hi = (ra, rb) if ra < rb else (rb, ra)
        self.parent[hi] = lo               # smallest id wins — deterministic
        return lo


class InstanceStitcher:
    """Accumulates window decodes into one slide-wide instance map."""

    def __init__(self, height: int, width: int, stitch_iou: float = 0.5,
                 dtype=np.int32):
        self.labels = np.zeros((height, width), dtype=dtype)
        self.covered = np.zeros((height, width), dtype=bool)
        self.stitch_iou = float(stitch_iou)
        self.uf = UnionFind()
        self.class_of_id: dict[int, int] = {}
        self.score_of_id: dict[int, float] = {}
        self.area_of_id: dict[int, int] = {}
        self._next_id = 1

    # ------------------------------------------------------------------- add
    def add_window(
        self,
        x0: int,
        y0: int,
        local_labels: np.ndarray,
        class_of_id: dict[int, int],
        score_of_id: Optional[dict[int, float]] = None,
    ) -> None:
        h, w = local_labels.shape
        y1, x1 = y0 + h, x0 + w
        g_view = self.labels[y0:y1, x0:x1]
        cov_view = self.covered[y0:y1, x0:x1]
        prev = cov_view & (g_view > 0)
        score_of_id = score_of_id or {}

        for lid, sl in enumerate(ndi.find_objects(local_labels), start=1):
            if sl is None or lid not in class_of_id:
                continue
            m = local_labels[sl] == lid
            area = int(m.sum())
            if area == 0:
                continue

            gid = self._next_id
            self._next_id += 1
            self.uf.add(gid)
            self.class_of_id[gid] = int(class_of_id[lid])
            self.score_of_id[gid] = float(score_of_id.get(lid, 0.0))
            self.area_of_id[gid] = area

            g_sub = g_view[sl]
            prev_sub = prev[sl]
            hit = m & prev_sub
            if hit.any():
                match = self._best_match(g_sub, m, hit)
                if match is not None:
                    self.uf.union(match, gid)

            # Fill only what is still empty; earlier windows keep their pixels.
            write = m & (g_sub == 0)
            g_sub[write] = gid
            g_view[sl] = g_sub

        cov_view[:] = True   # cov_view is a view into self.covered

    def _best_match(self, g_sub: np.ndarray, m: np.ndarray, hit: np.ndarray) -> Optional[int]:
        vals = g_sub[hit]
        vals = vals[vals > 0]
        if vals.size == 0:
            return None
        counts = np.bincount(vals)
        cand = int(np.argmax(counts))
        inter = int(counts[cand])
        local_in_covered = int(hit.sum())
        global_in_window = int((g_sub == cand).sum())
        denom = max(min(local_in_covered, global_in_window), 1)
        return cand if inter / denom >= self.stitch_iou else None

    # --------------------------------------------------------------- finalize
    def finalize(self) -> tuple[np.ndarray, dict[int, int], dict[int, float]]:
        """Resolve the union-find and relabel to contiguous ids starting at 1.

        A merged instance takes the class of its LARGEST contributing fragment: the window
        that saw most of the object had the most context for the vote.
        """
        if self._next_id == 1:
            return self.labels, {}, {}

        roots = np.arange(self._next_id, dtype=np.int64)
        for gid in range(1, self._next_id):
            roots[gid] = self.uf.find(gid)

        # Aggregate class by largest-fragment vote and score by area-weighted mean.
        best_area: dict[int, int] = {}
        cls: dict[int, int] = {}
        score_num: dict[int, float] = {}
        score_den: dict[int, float] = {}
        for gid in range(1, self._next_id):
            r = int(roots[gid])
            a = self.area_of_id.get(gid, 0)
            if a > best_area.get(r, -1):
                best_area[r] = a
                cls[r] = self.class_of_id.get(gid, 0)
            score_num[r] = score_num.get(r, 0.0) + self.score_of_id.get(gid, 0.0) * a
            score_den[r] = score_den.get(r, 0.0) + a

        # Compact only the roots that actually survive in the map: an instance whose every
        # pixel was claimed by an earlier window has no pixels left and must not get an id.
        used = (np.unique(roots[np.unique(self.labels[self.labels > 0])])
                if (self.labels > 0).any() else np.array([], dtype=np.int64))
        remap = np.zeros(self._next_id, dtype=np.int32)
        out_cls: dict[int, int] = {}
        out_score: dict[int, float] = {}
        for new_id, r in enumerate(sorted(int(u) for u in used), start=1):
            out_cls[new_id] = int(cls.get(r, 0))
            den = score_den.get(r, 0.0)
            out_score[new_id] = float(score_num.get(r, 0.0) / den) if den else 0.0
            remap[r] = new_id

        lut = remap[roots]                       # gid -> final id
        lut[0] = 0
        self.labels = lut[self.labels]
        return self.labels, out_cls, out_score
