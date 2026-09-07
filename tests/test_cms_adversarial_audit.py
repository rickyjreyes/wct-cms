"""Regression checks for independent score derivation and smooth-null witness."""
import pathlib,sys
import numpy as np
from scipy.stats import chi2
sys.path.insert(0,str(pathlib.Path(__file__).parents[1]/'experiments/cms_audit'))
from models import locked,poisson_fit,OMEGA,PHASE,MASKS,fit_model,eligible
from cms_wct.locked import fit_phase_locked_waveform
from cms_wct.background_kill import evaluate_pair_pipeline,PairKillConfig
from cms_wct.background_families import fit_background_family
ROOT=pathlib.Path(__file__).parents[1]

def test_independent_score_matches_reference():
    x=np.linspace(.7,4.8,287);rng=np.random.default_rng(22)
    for amp in (-1.,0.,1.):
        y=amp*np.cos(OMEGA*x-PHASE)+rng.normal(size=len(x));r=locked(x,y);ref=fit_phase_locked_waveform(x,y,OMEGA,PHASE)
        np.testing.assert_allclose([r['amplitude'],r['q']],[ref.signed_amplitude,ref.delta_chi2],atol=1e-11)

def test_poisson_fit_includes_zero_counts_and_bin_width():
    y=np.array([0.,2.,8.,0.]);width=np.array([1.,2.,3.,4.]);r=poisson_fit(y,np.ones((4,1)),np.log(width),np.ones(4,bool))
    np.testing.assert_allclose(r.mean,width*y.sum()/width.sum(),rtol=1e-9)
    assert r.converged

def test_bernstein_is_same_log_polynomial_space():
    m=np.geomspace(2,120,350);y=np.exp(5-.2*np.log(m)+.01*np.log(m)**2)
    a,_=fit_background_family(m,y,family='chebyshev',degree=7);b,_=fit_background_family(m,y,family='bernstein',degree=7)
    np.testing.assert_allclose(a,b,rtol=1e-10)

def test_phase_selection_invalidates_fixed_phase_reference():
    assert .008<chi2.sf(chi2.isf(.002,1),2)<.009

def test_smooth_background_seeded_failure_witness():
    p=ROOT/'audit/2026-09-06';s=np.genfromtxt(p/'baseline/H2/spectrum.csv',delimiter=',',names=True);means=np.load(p/'attacks/background_means.npz');rng=np.random.default_rng(20461607)
    a=rng.poisson(means['H2__pois_poly12']).astype(float);b=rng.poisson(means['G1__pois_poly12']).astype(float)
    r,_,_=evaluate_pair_pipeline(s['mass_center'],a,s['mass_center'],b,None,PairKillConfig(OMEGA,PHASE,MASKS))
    assert r['primary_pair_score_min_locked_delta_chi2']>120
    assert r['selected_background']['name']=='spline_s2'

def test_poisson_degree12_removes_positive_g2_component():
    s=np.genfromtxt(ROOT/'audit/2026-09-06/baseline/G2/spectrum.csv',delimiter=',',names=True);m=s['mass_center'];f=fit_model('pois_poly12',m,s['count'],s['mass_high']-s['mass_low']);g=eligible(m)
    assert f.converged
    assert locked(np.log(m[g]),((s['count']-f.mean)/np.sqrt(f.mean))[g])['amplitude']<0
