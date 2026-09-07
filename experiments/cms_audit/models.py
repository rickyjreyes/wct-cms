"""Independent audit models. No historical fit or acceptance rule is modified."""
from dataclasses import dataclass
import numpy as np
from scipy.interpolate import BSpline
from scipy.special import xlogy
from scipy.stats import norm

OMEGA = 7.025825825825827
PHASE = -0.2313916852932179
MASKS = ((2.9, 3.3), (3.55, 3.85), (8.5, 11.5), (80., 100.))


def eligible(m, windows=MASKS):
    good = np.isfinite(m) & (m > 0)
    for lo, hi in windows:
        good &= ~((m >= lo) & (m <= hi))
    return good


def deviance(y, mu):
    mu = np.maximum(mu, 1e-12)
    return float(2 * np.sum(mu-y+xlogy(y, y/mu)))


def locked(x, r, omega=OMEGA, phase=PHASE):
    """Independent scalar reconstruction of the historical statistic."""
    g = np.cos(omega*x-phase)
    g = g-g.mean()
    den = g@g
    amp = float(g@(r-r.mean())/den)
    q = max(amp, 0)**2*den
    return dict(amplitude=amp, q=float(q), naive_p=float(norm.sf(np.sqrt(q))) if amp>0 else 1., naive_se=float(den**-.5))


@dataclass
class Fit:
    mean: np.ndarray
    edf: float
    design: np.ndarray | None = None
    penalty: np.ndarray | None = None
    beta: np.ndarray | None = None
    converged: bool = True


def poisson_fit(y, X, offset, train, penalty=None, initial=None):
    """Convex Poisson log-link fit, with explicit Newton convergence check.

    Fixed penalty = Gaussian prior on spline second differences. No clipping,
    log-count response, exclusion of zero counts, or WCT template in the fit.
    """
    y = np.asarray(y, float)
    P = np.zeros((X.shape[1], X.shape[1])) if penalty is None else penalty
    beta = np.linalg.lstsq(X[train], np.log(y[train]+.5)-offset[train], rcond=1e-11)[0] if initial is None else initial.copy()
    Xt = X[train]; yt = y[train]; off = offset[train]
    def objective(b):
        eta = off+Xt@b
        if np.max(eta)>60: return np.inf
        return float(np.sum(np.exp(eta)-yt*eta)+.5*b@P@b)
    converged = False
    for _ in range(100):
        mu = np.exp(np.clip(off+Xt@beta,-50,60))
        grad = Xt.T@(mu-yt)+P@beta
        H = Xt.T@(mu[:,None]*Xt)+P
        step = np.linalg.lstsq(H, grad, rcond=1e-12)[0]
        old = objective(beta)
        scale = 1.
        while objective(beta-scale*step)>old+1e-9 and scale>2**-30:
            scale *= .5
        beta -= scale*step
        if np.max(np.abs(Xt@(scale*step)))<1e-8:
            converged = True
            break
    mean = np.exp(np.clip(offset+X@beta,-50,60))
    H0 = Xt.T@(mean[train,None]*Xt)
    edf = np.trace(np.linalg.pinv(H0+P,rcond=1e-11)@H0)
    return Fit(mean, float(edf), X, P, beta, converged)


def design_for(name, centers):
    x = np.log(centers); t = 2*(x-x.min())/np.ptp(x)-1
    if name.startswith('pois_poly'):
        d = int(name.split('poly')[1]); X = np.polynomial.chebyshev.chebvander(t,d)
        return X, np.zeros((d+1,d+1))
    if name.startswith('pspline'):
        n, lam = name.split('_')[1:]; n=int(n); lam=float(lam)
        knots = np.r_[[x.min()]*4, np.linspace(x.min(),x.max(),n+2)[1:-1], [x.max()]*4]
        X = BSpline.design_matrix(x,knots,3).toarray()
        D = np.diff(np.eye(X.shape[1]),2,axis=0)
        return X,lam*D.T@D
    if name=='cms_dijet_form':
        z=centers/13000.; l=np.log(z)
        X=np.column_stack([np.ones(len(x)),np.log1p(-z),l,l*l,l**3])
        # Basis preconditioning preserves the function space.
        X=np.linalg.qr(X)[0]
        return X,np.zeros((5,5))
    raise ValueError(name)


def fit_model(name, centers, counts, widths, train=None):
    train=eligible(centers) if train is None else np.asarray(train,bool)
    if name.startswith('local_'):
        bandwidth=float(name.split('_')[1]); x=np.log(centers)
        logdensity=np.log(counts+.5)-np.log(widths)
        result=[]; leverage=[]
        for i,xi in enumerate(x):
            X=np.column_stack([np.ones(len(x)),x-xi,(x-xi)**2])
            w=np.exp(-.5*((x-xi)/bandwidth)**2)*np.maximum(counts,1)*train
            h=np.linalg.pinv(X.T@(w[:,None]*X),rcond=1e-11)[0]@X.T*w
            result.append(np.exp(np.clip(h@logdensity,-50,50))*widths[i]);leverage.append(h[i])
        return Fit(np.array(result),float(np.sum(leverage)))
    X,P=design_for(name,centers)
    return poisson_fit(counts,X,np.log(widths),train,P)


def candidates():
    return ([f'pois_poly{d}' for d in (2,3,5,7,9,12,16,20)]
            +[f'pspline_{n}_{s}' for n in (8,16,24,40) for s in (.1,1,10,100)]
            +['cms_dijet_form']+[f'local_{h}' for h in (.15,.3,.5)])


def profile_locked(name, centers, counts, widths):
    """Joint Poisson profile diagnostic, not the historical residual score.

    Alternative log mu = log width + X beta + A cos(...)/sqrt(B0).
    Template scale is anchored at the null fit. A is approximately in Pearson
    units near zero. Both signs are fitted; discovery statistic is positive-only.
    Frozen penalties are retained, so penalized results are NOT Wilks p-values.
    """
    good=eligible(centers)
    null=fit_model(name,centers,counts,widths,good)
    g=np.cos(OMEGA*np.log(centers)-PHASE)/np.sqrt(null.mean)
    X=np.column_stack([null.design,g]); P=np.zeros((X.shape[1],X.shape[1]));P[:-1,:-1]=null.penalty
    alt=poisson_fit(counts,X,np.log(widths),good,P,np.r_[null.beta,0.])
    improvement=deviance(counts[good],null.mean[good])-deviance(counts[good],alt.mean[good])
    info=X[good].T@(alt.mean[good,None]*X[good])+P
    se=float(np.sqrt(np.linalg.pinv(info,rcond=1e-11)[-1,-1]))
    return dict(amplitude=float(alt.beta[-1]),se=se,ci95=[float(alt.beta[-1]-1.96*se),float(alt.beta[-1]+1.96*se)],deviance_improvement=float(improvement),positive_q=float(max(0,improvement)) if alt.beta[-1]>0 else 0.,converged=alt.converged and null.converged,penalized=bool(np.any(P)))
