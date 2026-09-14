"""No-private-data coverage for the shared nephron segment vocabulary."""
from __future__ import annotations

import pytest

from pseudospace.vocabulary import (
    COARSE_FAMILIES,
    FINE_SEGMENTS,
    KEEP_TUBULE_CLASSES,
    LABEL_VOCABULARY,
    REMOVE_CLASSES,
    SEGMENT_DISPLAY_ORDER,
    SEGMENT_TO_COARSE,
    assert_consistent,
    coarse_for,
    invalid_labels,
    stale_cluster_ids,
    unlabelled_clusters,
)

# The four segments the AL family deliberately merges. If this changes, the reported coarse classes
# change, so it is asserted rather than left to be rediscovered.
ASCENDING_LIMB_MEMBERS = ("ATL", "mTAL", "cTAL", "Macula-densa")


def test_definitions_are_internally_consistent():
    assert_consistent()


def test_segments_and_families_have_no_duplicates():
    assert len(set(FINE_SEGMENTS)) == len(FINE_SEGMENTS)
    assert len(set(COARSE_FAMILIES)) == len(COARSE_FAMILIES)


def test_every_fine_segment_rolls_up():
    for segment in FINE_SEGMENTS:
        assert segment in SEGMENT_TO_COARSE, segment


def test_rollup_targets_are_families_or_glomerulus():
    allowed = set(COARSE_FAMILIES) | {"Glomerulus"}
    assert set(SEGMENT_TO_COARSE.values()) <= allowed


def test_each_family_is_itself_a_label():
    for family in COARSE_FAMILIES:
        assert SEGMENT_TO_COARSE[family] == family


def test_al_is_the_merged_ascending_limb():
    for segment in ASCENDING_LIMB_MEMBERS:
        assert SEGMENT_TO_COARSE[segment] == "AL", segment


def test_thick_ascending_limb_is_not_reported_as_its_own_family():
    # The reported set is the five-family one, so AL is the only ascending-limb class.
    assert "TAL" not in COARSE_FAMILIES
    assert "TAL" not in LABEL_VOCABULARY


def test_label_vocabulary_covers_every_label_exactly_once():
    for label in FINE_SEGMENTS + COARSE_FAMILIES + REMOVE_CLASSES:
        assert LABEL_VOCABULARY.count(label) == 1, label


def test_keep_tubule_classes_are_the_coarse_families():
    assert KEEP_TUBULE_CLASSES == COARSE_FAMILIES


def test_display_order_places_each_family_before_its_members():
    position = {label: i for i, label in enumerate(SEGMENT_DISPLAY_ORDER)}
    for segment in FINE_SEGMENTS:
        family = SEGMENT_TO_COARSE[segment]
        if family in COARSE_FAMILIES:
            assert position[family] < position[segment], segment


def test_display_order_contains_every_label_once():
    assert len(set(SEGMENT_DISPLAY_ORDER)) == len(SEGMENT_DISPLAY_ORDER)
    assert set(LABEL_VOCABULARY) <= set(SEGMENT_DISPLAY_ORDER)


@pytest.mark.parametrize(
    "label, family",
    [
        ("PT-S2", "PT"),
        ("DTL3", "DTL"),
        ("ATL", "AL"),
        ("mTAL", "AL"),
        ("cTAL", "AL"),
        ("Macula-densa", "AL"),
        ("DCT1", "DCT"),
        ("IMCD", "CNT_CD"),
        ("Podocyte", "Glomerulus"),
        ("PT", "PT"),
        ("CNT_CD", "CNT_CD"),
        ("Glomerulus", "Glomerulus"),
    ],
)
def test_coarse_for_maps_known_labels(label, family):
    assert coarse_for(label) == family


def test_coarse_for_falls_back_to_unassigned():
    # Unmapped labels are reported, not raised: a cluster may carry a label this module does not
    # know and still needs to appear in the pass-1 object.
    assert coarse_for("NotASegment") == "Unassigned"


def test_invalid_labels_flags_only_unknown_labels():
    assert invalid_labels(["PT-S1", "AL", "Unassigned"]) == []
    assert invalid_labels(["PT-S1", "cTAL"]) == []
    assert invalid_labels(["PT-S1", "Nope", "AlsoNope"]) == ["AlsoNope", "Nope"]


def test_unlabelled_and_stale_cluster_helpers():
    clusters = ["0", "1", "2"]
    assert unlabelled_clusters(clusters, {"0": "PT", "1": "AL", "2": "DCT"}) == []
    assert unlabelled_clusters(clusters, {"0": "PT"}) == ["1", "2"]
    assert stale_cluster_ids(clusters, {"0": "PT", "9": "DCT"}) == ["9"]
    assert stale_cluster_ids(clusters, {"0": "PT", "1": "AL", "2": "DCT"}) == []
