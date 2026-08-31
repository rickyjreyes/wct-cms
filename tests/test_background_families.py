import numpy as np

from cms_wct.background_families import fit_background_family


def synthetic_spectrum():
    centers = np.geomspace(2.0, 120.0, 180)
    x = np.log(centers)
    model = np.exp(6.0 - 0.7 * x + 0.04 * x * x)
    counts = model.copy()
    return centers, counts


def test_background_families_return_positive_models():
    centers, counts = synthetic_spectrum()
    for family in ("chebyshev", "bernstein", "spline"):
        model, mask = fit_background_family(
            centers,
            counts,
            family=family,
            degree=7,
            iterations=3,
            clip_sigma=4.0,
            excluded_windows=[(2.9, 3.3), (8.5, 11.5)],
            spline_smoothing_factor=1.0,
        )
        assert model.shape == counts.shape
        assert mask.shape == counts.shape
        assert np.all(np.isfinite(model))
        assert np.all(model > 0.0)
        assert np.count_nonzero(mask) > 100


def test_unknown_background_family_fails():
    centers, counts = synthetic_spectrum()
    try:
        fit_background_family(centers, counts, family="not-a-model")
    except ValueError as exc:
        assert "Unknown background family" in str(exc)
    else:
        raise AssertionError("unknown family should fail")
