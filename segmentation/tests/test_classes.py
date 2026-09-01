import pytest

from kidney_panoptic.data.classes import ClassSpec


def test_tubule_variants_map_to_distinct_anatomical_segment_classes(spec):
    """They are anatomical segments, kept separate so the decode can split on the border."""
    assert spec.canonical("Tubules") == "tubule_proximal"
    assert spec.canonical("Tubules2") == "tubule_distal"
    assert spec.canonical("Tubules3") == "tubule_collecting"
    assert len({spec.canonical_index(n)
                for n in ("Tubules", "Tubules2", "Tubules3")}) == 3


def test_weak_classes_are_ignored_not_trained(spec):
    """v3 folded Vessel/Blood Cells/Nerve into one `other_structure` class and it collapsed
    onto its smallest member: predicted areas capped at p95 1847 px against annotations
    running to 26022 px, i.e. it never predicted a vessel, while emitting 14% of the slide
    output as objects that are not per-object interpretable. They are now IGNORED, so those
    pixels contribute no gradient in either direction."""
    for raw in ("Vessel", "Blood Cells", "Nerve"):
        assert spec.role(raw) == "ignore"
        assert spec.canonical(raw) is None
    for gone in ("other_structure", "vessel", "blood_cell", "nerve"):
        assert gone not in spec.names


def test_all_labels_present_in_the_real_geojson_resolve(spec):
    """Every class name that appears in the cohort's annotations must have a role."""
    for label in ("Tubules", "Tubules2", "Tubules3", "Glomeruli", "Vessel",
                  "Blood Cells", "Nerve", "Background"):
        assert spec.role(label) in {"instance", "background", "dense_roi", "ignore"}
        if spec.role(label) == "instance":
            assert spec.canonical(label) in spec.names


def test_background_is_not_an_instance_class(spec):
    assert spec.role("Background") == "background"
    assert spec.canonical("Background") is None
    assert "Background" not in spec.names


def test_dense_roi_patterns_match_by_glob(spec):
    assert spec.role("ROI") == "dense_roi"
    assert spec.role("Region 3") == "dense_roi"


def test_unmapped_instance_label_fails_loudly(spec):
    with pytest.raises(KeyError, match="no entry in class_map"):
        spec.canonical("Podocytes")


def test_background_must_be_index_zero():
    with pytest.raises(ValueError, match="classes\\[0\\]"):
        ClassSpec.from_config({"classes": ["tubule", "background"]})


def test_class_map_target_must_exist():
    with pytest.raises(ValueError, match="not in `classes`"):
        ClassSpec.from_config({"classes": ["background", "tubule"],
                               "class_map": {"Glomeruli": "glomerulus"}})


def test_compact_and_elongated_must_be_known_classes():
    with pytest.raises(ValueError, match="compact_classes"):
        ClassSpec.from_config({"classes": ["background", "tubule"],
                               "compact_classes": ["glomerulus"]})


def test_compact_indices_align_with_center_channels(spec):
    assert spec.compact_indices == (spec.index("glomerulus"),)
    assert len(spec.compact_indices) == len(spec.compact)


def test_only_glomerulus_is_compact(spec):
    """Compact membership enables a round-object distance-transform split. Applying that to
    an elongated tubule cross-section would shatter it (ADR 0005), so the tubule segments
    must stay out of it."""
    assert spec.compact == ("glomerulus",)
    for t in ("tubule_proximal", "tubule_distal", "tubule_collecting"):
        assert t not in spec.compact
        assert spec.index(t) not in spec.compact_indices


def test_thing_names_exclude_background(spec):
    assert "background" not in spec.thing_names
    assert len(spec.thing_names) == spec.n_classes - 1


def test_canonical_classes_map_to_themselves_without_an_explicit_entry():
    s = ClassSpec.from_config({"classes": ["background", "tubule"]})
    assert s.canonical("tubule") == "tubule"
