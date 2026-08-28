import numpy as np

from cms_wct.nulls import parametric_background_null


def test_parametric_background_null_is_reproducible_and_bounded():
    centers = np.geomspace(2.0, 120.0, 80)
    x = np.log(centers)
    null_model = 180.0 * np.exp(-0.35 * (x - x.min())) + 20.0
    mask = np.ones_like(centers, dtype=bool)
    omegas = np.linspace(3.1, 15.0, 80)

    kwargs = dict(
        centers=centers,
        null_model=null_model,
        analysis_mask=mask,
        x=x,
        omegas=omegas,
        observed_best_score=25.0,
        frozen_omega=7.0,
        observed_frozen_score=10.0,
        degree=3,
        iterations=3,
        clip_sigma=4.0,
        excluded_windows=(),
        n_bootstrap=12,
        seed=12345,
    )

    out1 = parametric_background_null(**kwargs)
    out2 = parametric_background_null(**kwargs)

    p_global, p_frozen, max_scores, frozen_scores = out1
    assert 0.0 < p_global <= 1.0
    assert p_frozen is not None and 0.0 < p_frozen <= 1.0
    assert max_scores.shape == (12,)
    assert frozen_scores is not None and frozen_scores.shape == (12,)
    assert np.all(np.isfinite(max_scores))
    assert np.all(np.isfinite(frozen_scores))
    assert np.allclose(max_scores, out2[2])
    assert np.allclose(frozen_scores, out2[3])
