import numpy as np
import pytest

from kidney_panoptic.data.augment import (
    Sample, augment_sample, copy_paste, photometric, random_affine, random_flip_rot90,
)
from kidney_panoptic.data.targets import build_targets


def make_sample(h=64, w=64):
    rng = np.random.default_rng(0)
    img = rng.integers(80, 220, (h, w, 3), dtype=np.uint8)
    lbl = np.zeros((h, w), np.int32)
    lbl[10:30, 10:30] = 1
    lbl[35:55, 35:55] = 2
    return Sample(img, lbl, np.ones((h, w), bool), np.zeros((h, w), bool),
                  np.zeros((h, w), bool), np.ones((h, w), bool), {1: 1, 2: 2})


def test_geometry_is_applied_identically_to_image_and_every_mask():
    s = make_sample()
    rng = np.random.default_rng(3)
    out = random_affine(s, rng, {"rotate_deg": 90, "scale": [1.0, 1.0], "translate_px": 0})
    assert out.image.shape == s.image.shape
    for m in (out.id_map, out.tissue, out.bg_poly, out.dense_roi, out.inbounds):
        assert m.shape == s.id_map.shape


def test_masks_stay_integral_under_warping():
    """Nearest interpolation only — a blended label id would be a new, fake instance."""
    s = make_sample()
    out = random_affine(s, np.random.default_rng(1),
                        {"rotate_deg": 37, "scale": [0.9, 1.1], "translate_px": 5})
    assert set(np.unique(out.id_map)) <= {0, 1, 2}
    assert out.id_map.dtype == np.int32


def test_rotation_marks_out_of_frame_pixels_as_not_inbounds():
    """Corners pulled in from outside the patch are unlabelled and must be ignorable."""
    s = make_sample()
    out = random_affine(s, np.random.default_rng(0),
                        {"rotate_deg": 45, "scale": [1.0, 1.0], "translate_px": 0})
    assert not out.inbounds.all(), "a 45-degree rotation must expose out-of-frame corners"
    assert out.inbounds[32, 32], "the centre is always in bounds"


def test_no_geometry_leaves_inbounds_untouched():
    s = make_sample()
    out = random_affine(s, np.random.default_rng(0),
                        {"rotate_deg": 0, "scale": [1.0, 1.0], "translate_px": 0})
    assert out.inbounds.all()


def test_flip_is_exact_and_reversible():
    s = make_sample()
    out = random_flip_rot90(s, np.random.default_rng(0),
                            {"hflip": 1.0, "vflip": 0.0, "rot90": 0.0})
    assert np.array_equal(out.id_map, s.id_map[:, ::-1])
    assert np.array_equal(out.image, s.image[:, ::-1])


def test_targets_regenerated_after_geometry_stay_consistent_with_the_image():
    """The invariant the whole design rests on: targets are a function of the warped map."""
    s = make_sample()
    out = random_flip_rot90(s, np.random.default_rng(0),
                            {"hflip": 1.0, "vflip": 0.0, "rot90": 0.0})
    t = build_targets(out.id_map, out.class_idx_of_id, [2])
    assert np.array_equal(t.semantic > 0, out.id_map > 0)
    assert ((t.boundary > 0) == (out.id_map > 0)).all()


def test_photometric_changes_pixels_only():
    s = make_sample()
    img = photometric(s.image, np.random.default_rng(0),
                      {"hed_sigma": 0.05, "hsv_p": 1.0, "brightness_contrast_p": 1.0,
                       "blur_p": 1.0})
    assert img.shape == s.image.shape and img.dtype == np.uint8
    assert not np.array_equal(img, s.image)


def test_copy_paste_adds_a_new_instance_with_the_donor_class():
    s = make_sample()
    donor = make_sample()
    out = copy_paste(s, donor.image, donor.id_map, 1, donor_class_idx=3,
                     rng=np.random.default_rng(2), blend_px=1)
    if out.id_map.max() > s.id_map.max():           # a paste site was found
        new_id = int(out.id_map.max())
        assert out.class_idx_of_id[new_id] == 3
        assert (out.id_map == new_id).sum() > 0


def test_copy_paste_does_not_bury_an_existing_annotated_instance():
    """A paste that would destroy a real labelled object must be rejected."""
    s = make_sample(32, 32)
    s.id_map[:] = 1                                  # every pixel already annotated
    donor = make_sample(32, 32)
    out = copy_paste(s, donor.image, donor.id_map, 1, 3, np.random.default_rng(0),
                     max_overlap_frac=0.15)
    assert out.id_map.max() == 1, "no paste should have happened"


def test_full_pipeline_preserves_shapes_and_dtypes():
    s = make_sample()
    cfg = {"hflip": 0.5, "vflip": 0.5, "rot90": 0.5, "rotate_deg": 180,
           "translate_px": 8, "scale": [0.85, 1.15], "hed_sigma": 0.05, "hsv_p": 0.5,
           "brightness_contrast_p": 0.5, "blur_p": 0.2, "elastic_p": 1.0,
           "elastic_alpha": 8.0, "elastic_sigma": 4.0}
    for seed in range(5):
        out = augment_sample(make_sample(), np.random.default_rng(seed), cfg)
        assert out.image.shape == s.image.shape and out.image.dtype == np.uint8
        assert out.id_map.shape == s.id_map.shape and out.id_map.dtype == np.int32
        assert set(np.unique(out.id_map)) <= {0, 1, 2}


# ------------------------------------------------------ cutout / stain / shear
def test_cutout_marks_the_hole_ignored_and_never_as_background():
    """The whole point: a flat erased rectangle must NOT be supervised as background."""
    from kidney_panoptic.data.augment import random_cutout

    s = make_sample(64, 64)
    out = random_cutout(s, np.random.default_rng(0),
                        {"cutout_n": [2, 2], "cutout_frac": [0.05, 0.05]})
    assert not out.inbounds.all(), "erased pixels must be marked out of bounds"
    assert np.array_equal(out.id_map, s.id_map), "cutout must not touch the instance map"
    assert np.array_equal(out.tissue, s.tissue)
    erased = ~out.inbounds
    assert not np.array_equal(out.image[erased], s.image[erased])


def test_cutout_hole_is_stripped_from_the_valid_mask():
    """End-to-end contract with the dataset: inbounds is ANDed into `valid`."""
    from kidney_panoptic.data.augment import random_cutout
    from kidney_panoptic.data.ignore import build_valid_mask

    s = make_sample(64, 64)
    out = random_cutout(s, np.random.default_rng(1),
                        {"cutout_n": [1, 1], "cutout_frac": [0.10, 0.10]})
    valid = build_valid_mask(out.id_map, strategy="none") & out.inbounds
    assert not valid[~out.inbounds].any()


def test_hed_beta_is_interpreted_in_optical_density_units():
    """Regression guard on a real trap: beta is NOT dimensionless.

    hed2rgb is exponential, so an OD offset that looks tiny next to the channel range
    produces a large RGB shift. beta ~ 0.004 is a sane stain jitter; beta ~ 0.05 destroys
    the image. Measured mean |dRGB|: 0.004 -> a few units, 0.05 -> tens of units.
    """
    from kidney_panoptic.data.augment import _HAVE_HED, hed_jitter

    if not _HAVE_HED:
        pytest.skip("skimage HED unavailable")

    # An H&E-like image, not uniform noise: eosinophilic pink background with
    # haematoxylin-dark nuclei. Optical density depends strongly on stain content, so
    # random noise would not exercise the property under test.
    img = np.full((64, 64, 3), (235, 190, 215), np.uint8)
    img[20:40, 20:40] = (110, 70, 150)

    def delta(beta):
        return float(np.abs(hed_jitter(img, np.random.default_rng(0), 0.0, beta).astype(float)
                            - img.astype(float)).mean())

    d_small, d_large = delta(0.004), delta(0.05)
    assert delta(0.0) == 0.0, "beta=0 must be a no-op"
    assert d_large > 5 * d_small, "beta is a scale parameter, and 0.05 is not 'a bit more'"
    assert d_small < 0.25 * 255, "beta=0.004 must remain a recognisable stain shift"


def test_legacy_hed_sigma_key_keeps_its_original_meaning():
    """An old config must not silently become ~20x stronger."""
    from kidney_panoptic.data.augment import _hed_params

    assert _hed_params({"hed_sigma": 0.05}) == (0.05, 0.05 * 0.05)
    assert _hed_params({"hed_alpha_sigma": 0.1, "hed_beta_od": 0.004}) == (0.1, 0.004)


def test_shear_and_anisotropic_scale_keep_image_and_masks_aligned():
    """A square instance must still coincide with its image footprint after shearing."""
    s = make_sample(96, 96)
    s.image[:] = 0
    s.image[s.id_map == 1] = 255
    out = random_affine(s, np.random.default_rng(4),
                        {"rotate_deg": 0, "scale": [1.0, 1.0], "translate_px": 0,
                         "aspect": [0.9, 1.11], "shear_deg": 6})
    bright = out.image[..., 0] > 128
    inst = out.id_map == 1
    overlap = (bright & inst).sum() / max(inst.sum(), 1)
    assert overlap > 0.9, f"mask drifted from the image under shear (overlap {overlap:.2f})"


def test_noise_and_jpeg_only_touch_the_image():
    img = make_sample(64, 64).image
    out = photometric(img, np.random.default_rng(0),
                      {"noise_p": 1.0, "noise_sigma": [5.0, 5.0], "jpeg_p": 1.0,
                       "jpeg_quality": [60, 60]})
    assert out.shape == img.shape and out.dtype == np.uint8
    assert not np.array_equal(out, img)
