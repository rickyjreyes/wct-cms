import math

from cms_wct.significance import (
    add_one_monte_carlo_p,
    clopper_pearson_upper_bound,
    five_sigma_requirements,
    minimum_trials_for_add_one_resolution,
    minimum_zero_exceedance_trials,
    one_sided_gaussian_p,
)


def test_five_sigma_one_sided_tail():
    p5 = one_sided_gaussian_p(5.0)
    assert math.isclose(p5, 2.866515718791933e-7, rel_tol=1e-12)


def test_add_one_resolution_count_for_five_sigma():
    p5 = one_sided_gaussian_p(5.0)
    n = minimum_trials_for_add_one_resolution(p5)
    assert n == 3_488_555
    assert add_one_monte_carlo_p(0, n) <= p5
    assert add_one_monte_carlo_p(0, n - 1) > p5


def test_zero_exceedance_confidence_requirements():
    p5 = one_sided_gaussian_p(5.0)
    n95 = minimum_zero_exceedance_trials(p5, 0.95)
    n99 = minimum_zero_exceedance_trials(p5, 0.99)
    assert n95 == 10_450_778
    assert n99 == 16_065_391
    assert clopper_pearson_upper_bound(0, n95, 0.95) <= p5
    assert clopper_pearson_upper_bound(0, n95 - 1, 0.95) > p5


def test_requirement_snapshot():
    req = five_sigma_requirements()
    assert req["minimum_trials_add_one_resolution"] == 3_488_555
    assert req["minimum_zero_exceedance_trials_95pct"] == 10_450_778
    assert req["minimum_zero_exceedance_trials_99pct"] == 16_065_391
