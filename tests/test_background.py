import numpy as np

from cms_wct.background import make_histogram, robust_smooth_background


def test_histogram_count_conservation():
    masses = np.array([1.0, 2.0, 2.5, 3.0, 4.0, 8.0])
    counts, _, _ = make_histogram(masses, 2.0, 4.0, 4, False)
    assert int(counts.sum()) == 4


def test_background_returns_positive_model():
    centers = np.geomspace(2.0, 100.0, 100)
    counts = 10000.0 * centers ** -1.5 + 5.0
    model, mask = robust_smooth_background(centers, counts, degree=4)
    assert np.all(model > 0)
    assert mask.sum() > 50
