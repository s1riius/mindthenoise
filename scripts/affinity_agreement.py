"""
affinity_agreement.py - statistical analysis of agreement between affinity metrics (Supplementary S7).
Reproduces all numbers of Tables S3, S5, and S6-S10 from two input files (Figures S15, S18 and S19: figures.py):
  vegfr.csv     columns: compound, Ki (nM), IC50 (nM, mean of reports), stdev (nM; 0 = single report)      [Table S11]
  lindberg.csv  columns: compound, Kd_D1A, Kd_CLK1, Kd_GSK3B (nM), R_D1A, R_CLK1, R_GSK3B (radiometric IC50, uM),
                P_D1A, P_CLK1 (ADP-Glo IC50, uM); empty = not available                                      [Table S12]
Usage:  python affinity_agreement.py vegfr  |  lindberg  |  ratios  |  extras
Conventions: natural logarithms; affinities in uM. Literature uncertainties in pK (log10) units are multiplied by ln 10.
Requires Python >= 3.9, numpy, scipy, pandas, matplotlib, scikit-learn. Random seeds are fixed.
"""
import sys
import numpy as np, pandas as pd
from scipy import stats, optimize
LN10 = np.log(10)
S_REP = 0.2 * LN10                      # within-assay repeatability, 0.2 pK units -> 0.46 ln units
S_KI, S_IC = 0.54 * LN10, 0.68 * LN10   # inter-laboratory reproducibility of pKi (Kramer) and pIC50 (Kalliokoski)
NB, NNULL = 2000, 3000                  # bootstrap replicates, null-model datasets
# ---------------------------------------------------------------- basic estimators
def deming(x, y, lam=1.0):
    mx, my = x.mean(), y.mean(); sxx = ((x-mx)**2).mean(); syy = ((y-my)**2).mean(); sxy = ((x-mx)*(y-my)).mean()
    return (syy - lam*sxx + np.sqrt((syy - lam*sxx)**2 + 4*lam*sxy**2)) / (2*sxy)
def passing_bablok(x, y):
    i, j = np.triu_indices(len(x), 1); dx, dy = x[j]-x[i], y[j]-y[i]; ok = dx != 0
    s = np.sort(dy[ok]/dx[ok]); s = s[s != -1]; n = len(s); k = int((s < -1).sum())
    return s[(n+1)//2+k-1] if n % 2 else 0.5*(s[n//2+k-1] + s[n//2+k])
def ccc(x, y):
    return 2*np.cov(x, y, bias=True)[0, 1] / (x.var() + y.var() + (x.mean()-y.mean())**2)
def r2(x, y):
    return stats.pearsonr(x, y)[0]**2 if len(x) > 2 else np.nan
def cooks(x, y):
    X = np.c_[np.ones(len(x)), x]; h = np.diag(X @ np.linalg.inv(X.T @ X) @ X.T)
    s, i = stats.linregress(x, y)[:2]; res = y - (i + s*x); mse = (res**2).sum()/(len(x)-2)
    return res**2/(2*mse)*h/(1-h)**2
def window_r2(x, y, w, step=0.2, minpts=5):
    v = [r2(x[m], y[m]) for c in np.arange(x.min(), x.max(), step) if (m := (x >= c) & (x < c+w)).sum() >= minpts]
    return np.mean(v) if v else np.nan
def inversion_beta(x, y):
    """P(inversion | separation s) = 1/(1+exp(beta*s)); returns beta (maximum likelihood over all pairs)."""
    i, j = np.triu_indices(len(x), 1); s = np.abs(x[j]-x[i]); inv = (np.sign(x[j]-x[i]) != np.sign(y[j]-y[i])).astype(float)
    nll = lambda b: -np.sum(inv*np.log(1/(1+np.exp(b*s))+1e-12) + (1-inv)*np.log(1-1/(1+np.exp(b*s))+1e-12))
    return optimize.minimize_scalar(nll, bounds=(1e-3, 30), method='bounded').x
def threshold(beta, p):           # separation (ln units) at which the inversion probability equals p
    return np.log(1/p - 1)/beta
def inversion_bins(x, y, edges=(0, 1, 2, 3, 4, 5, 7, 30)):
    i, j = np.triu_indices(len(x), 1); s = np.abs(x[j]-x[i]); inv = np.sign(x[j]-x[i]) != np.sign(y[j]-y[i])
    return np.array([inv[(s >= a) & (s < b)].mean() if ((s >= a) & (s < b)).any() else np.nan for a, b in zip(edges[:-1], edges[1:])])
def boot(f, x, y, rng, nb=NB):
    v = []
    for _ in range(nb):
        k = rng.integers(0, len(x), len(x))
        if len(set(k)) < 4: continue
        try: v.append(f(x[k], y[k]))
        except Exception: pass
    return np.nanpercentile(v, [2.5, 97.5])
# ---------------------------------------------------------------- null model (position-preserving)
def null_sets(x, y, sx, sy, rng, n=NNULL):
    """Datasets in which both metrics measure the same quantity (up to the observed mean offset) and differ only by
    Gaussian measurement error of SD sx and sy; compound positions T are taken from the observed data."""
    b = (y-x).mean(); T = (x+y-b)/2
    for _ in range(n):
        yield T + rng.normal(0, sx, len(x)), T + b + rng.normal(0, sy, len(x))
def q95(v):
    v = np.asarray(v, float); v = v[~np.isnan(v)]; return np.percentile(v, [2.5, 50, 97.5])
# ---------------------------------------------------------------- full analysis of one pair
def analyse(x, y, label, lam, noise, rng, split=None, windows=(2, 3, 4), scenarios=None):
    n = len(x); D = y-x; out = dict(label=label, n=n)
    s, i, r = stats.linregress(x, y)[:3]
    out.update(R2=r*r, p=stats.pearsonr(x, y)[1], OLS=s, OLS_ci=boot(lambda a, b: stats.linregress(a, b)[0], x, y, rng),
               intercept=i, Dem=deming(x, y, lam), Dem_ci=boot(lambda a, b: deming(a, b, lam), x, y, rng),
               Dem_l1=deming(x, y, 1.0), Dem_l25=deming(x, y, 2.5), PB=passing_bablok(x, y),
               PB_ci=boot(passing_bablok, x, y, rng), CCC=ccc(x, y), CCC_ci=boot(ccc, x, y, rng),
               rho=stats.spearmanr(x, y)[0], tau=stats.kendalltau(x, y)[0],
               bias=D.mean(), bias_ci=boot(lambda a, b: (b-a).mean(), x, y, rng), sdD=D.std(ddof=1),
               LoA=(np.exp(D.mean()-1.96*D.std(ddof=1)), np.exp(D.mean()+1.96*D.std(ddof=1))),
               BA_trend=stats.linregress((x+y)/2, D)[0::3], neg=int((D < 0).sum()),
               var=D.var(ddof=1), between=D.var(ddof=1)-2*S_REP**2,
               between_ci=boot(lambda a, b: (b-a).var(ddof=1)-2*S_REP**2, x, y, rng),
               cook_max=cooks(x, y).max(), windows={w: window_r2(x, y, w) for w in windows})
    beta = inversion_beta(x, y); bb = [inversion_beta(x[k], y[k]) for k in (rng.integers(0, n, n) for _ in range(500)) if len(set(k)) > 3]
    out['thr'] = {p: (threshold(beta, p), np.percentile([threshold(v, p) for v in bb], [2.5, 97.5])) for p in (0.10, 0.05, 0.01)}
    out['inv_bins'] = inversion_bins(x, y)
    if split is not None:
        out['clusters'] = {lab: (int(k.sum()), r2(x[k], y[k]), stats.pearsonr(x[k], y[k])[1]) for lab, k in split.items()}
    out['null'] = {}
    for name, (sx, sy) in noise.items():
        acc = {k: [] for k in ['R2', 'Dem', 'sdD', 'neg', 'CCC', 'thr05', 'thr10', 'inv'] + [f'w{w}' for w in windows] + (list(split) if split else [])}
        for it, (xs, ys) in enumerate(null_sets(x, y, sx, sy, rng)):
            acc['R2'].append(r2(xs, ys)); acc['Dem'].append(deming(xs, ys, lam)); acc['sdD'].append((ys-xs).std(ddof=1))
            acc['neg'].append(((ys-xs) < 0).mean()); acc['CCC'].append(ccc(xs, ys))
            if split:
                for lab, k in split.items(): acc[lab].append(r2(xs[k], ys[k]))
            if it < 300:
                bt = inversion_beta(xs, ys); acc['thr05'].append(threshold(bt, 0.05)); acc['thr10'].append(threshold(bt, 0.10))
                acc['inv'].append(inversion_bins(xs, ys))
                for w in windows: acc[f'w{w}'].append(window_r2(xs, ys, w))
        out['null'][name] = {k: (np.nanpercentile(np.array(v), [2.5, 50, 97.5], axis=0) if k == 'inv' else q95(v)) for k, v in acc.items() if len(v)}
    if scenarios:
        out['scen'] = {}
        for sname, f in scenarios.items():
            xs_, ys_ = f(x, y); Ds = ys_-xs_
            o = dict(n=len(xs_), R2=r2(xs_, ys_), OLS=stats.linregress(xs_, ys_)[0], Dem=deming(xs_, ys_, lam), bias=Ds.mean(), sdD=Ds.std(ddof=1), CCC=ccc(xs_, ys_))
            nul = {}
            for name, (sx, sy) in noise.items():
                acc = {k: [] for k in ['R2', 'Dem', 'bias', 'sdD', 'CCC']}
                for xa, ya in null_sets(x, y, sx, sy, rng):
                    xb, yb = f(xa, ya)
                    if len(xb) < 5: continue
                    acc['R2'].append(r2(xb, yb)); acc['Dem'].append(deming(xb, yb, lam)); acc['bias'].append((yb-xb).mean())
                    acc['sdD'].append((yb-xb).std(ddof=1)); acc['CCC'].append(ccc(xb, yb))
                nul[name] = {k: q95(v) for k, v in acc.items()}
            out['scen'][sname] = (o, nul)
    return out
def report(o):
    f = lambda v, c: f'{v:.2f} [{c[0]:.2f}, {c[1]:.2f}]'
    print(f"\n######## {o['label']}  n = {o['n']}")
    print(f"R2 {o['R2']:.2f} (p = {o['p']:.2g}); OLS {f(o['OLS'], o['OLS_ci'])}, intercept {o['intercept']:.2f}; Deming {f(o['Dem'], o['Dem_ci'])} "
          f"(lambda=1: {o['Dem_l1']:.2f}, 2.5: {o['Dem_l25']:.2f}); Passing-Bablok {f(o['PB'], o['PB_ci'])}")
    print(f"CCC {f(o['CCC'], o['CCC_ci'])}; Spearman {o['rho']:.2f}; Kendall {o['tau']:.2f}; max Cook {o['cook_max']:.3f} (4/n = {4/o['n']:.3f})")
    print(f"bias ln(y/x) {f(o['bias'], o['bias_ci'])} = {np.exp(o['bias']):.1f}-fold [{np.exp(o['bias_ci'][0]):.1f}, {np.exp(o['bias_ci'][1]):.1f}]; "
          f"SD {o['sdD']:.2f}; limits of agreement {o['LoA'][0]:.2g}- to {o['LoA'][1]:.3g}-fold; BA trend {o['BA_trend'][0]:.2f} (p = {o['BA_trend'][1]:.2g}); n(y<x) {o['neg']}")
    print(f"variance {o['var']:.2f}; between-assay {o['between']:.2f} [{o['between_ci'][0]:.2f}, {o['between_ci'][1]:.2f}] "
          f"({100*o['between']/o['var']:.0f}%); between-assay SD per metric {np.sqrt(max(o['between'], 0)/2):.2f}")
    print('inversion thresholds: ' + '; '.join(f"{p:.0%}: {t:.1f} ln = {np.exp(t):.0f}-fold [{np.exp(c[0]):.0f}, {np.exp(c[1]):.0f}]" for p, (t, c) in o['thr'].items()))
    print('window R2: ' + ', '.join(f'{w}: {v:.2f}' for w, v in o['windows'].items()))
    for lab, (k, rr, pp) in o.get('clusters', {}).items(): print(f'cluster {lab}: n = {k}, R2 = {rr:.2f} (p = {pp:.2f})')
    for name, d in o['null'].items():
        g = lambda k: f"[{d[k][0]:.2f}, {d[k][2]:.2f}]"
        line = f"null {name}: R2 {g('R2')} Deming {g('Dem')} SD {g('sdD')} frac(y<x) {g('neg')} CCC {g('CCC')} "
        line += f"thr10 {np.exp(d['thr10'][1]):.0f}-fold thr05 {np.exp(d['thr05'][1]):.0f}-fold; windows " + ', '.join(f"{w[1:]}: {d[w][1]:.2f} {g(w)}" for w in d if w.startswith('w'))
        for lab in o.get('clusters', {}): line += f"; cluster {lab} {d[lab][1]:.2f} {g(lab)}"
        print(line)
    for sname, (ob, nul) in o.get('scen', {}).items():
        print(f"scenario {sname}: " + ', '.join(f'{k} {v:.2f}' if isinstance(v, float) else f'{k} {v}' for k, v in ob.items()))
        for name, d in nul.items(): print(f"    null {name}: " + ', '.join(f"{k} [{v[0]:.2f}, {v[2]:.2f}]" for k, v in d.items()))
# ---------------------------------------------------------------- datasets
def vegfr(path='vegfr.csv'):
    d = pd.read_csv(path); rng = np.random.default_rng(2026)
    x, y = np.log(d.Ki.values/1000), np.log(d.IC50.values/1000)
    lam = (S_IC/S_KI)**2; noise = {'within-assay': (S_REP, S_REP), 'inter-lab': (S_KI, S_IC)}
    line = y > 3 + (-10/10.6)*(x + 9)                      # dashed blue line of Figure 4, through (-9, 3) and (1.6, -7)
    onl = (d.compound == 'OSI-632').values                  # compound lying on the line, excluded from both clusters
    split = {'upper': line & ~onl, 'lower': ~line & ~onl}
    sdw = np.sqrt(2)*S_REP
    def sc(kind):
        def f(a, b):
            Dv = b-a; sig = Dv < -2*sdw; comp = (Dv < 0) & ~sig; bb = b.copy(); bb[comp] = a[comp]
            if kind == 1: return a[~sig], b[~sig]
            if kind == 2: return a[Dv >= 0], b[Dv >= 0]
            if kind == 3: return a[~sig], bb[~sig]
            return a, bb
        return f
    scen = {'S1 exclude IC50<KI beyond 2 SD': sc(1), 'S2 exclude all IC50<KI': sc(2),
            'S3 S1 + set compatible IC50=KI': sc(3), 'S4 keep all, set compatible IC50=KI': sc(4)}
    report(analyse(x, y, 'VEGFR: ln KI (x) vs ln IC50 (y), all compounds', lam, noise, rng, split, windows=(2, 3, 4, 5, 6, 7), scenarios=scen))
    k = ~d.compound.isin(['DIGITONIN', '2378-TCDD']).values
    report(analyse(x[k], y[k], 'VEGFR without digitonin and TCDD', lam, noise, rng))
    rep = (d.stdev > 0).values
    for lab, m in [('IC50 reported once', ~rep), ('IC50 replicated', rep)]:
        print(f"\nsubset {lab}: n = {m.sum()}, R2 = {r2(x[m], y[m]):.2f}, OLS = {stats.linregress(x[m], y[m])[0]:.2f}, SD ln KI = {x[m].std(ddof=1):.2f}")
        for name, (sx, sy) in noise.items():
            v = [r2(xs[m], ys[m]) for xs, ys in null_sets(x, y, sx, sy, rng)]; print(f"    null {name}: R2 {q95(v)[0]:.2f}-{q95(v)[2]:.2f}")
    cv = d.stdev[rep]/d.IC50[rep]; print(f"\nln-SD between IC50 reports (n = {rep.sum()}): median {np.median(np.sqrt(np.log(1+cv**2))):.2f}")
    print('compounds with IC50 < KI: ' + ', '.join(f"{c} ({v:.2f}, z = {v/sdw:.1f})" for c, v in zip(d.compound[(y-x) < 0], (y-x)[(y-x) < 0])))
    from sklearn.mixture import GaussianMixture as GM
    M = ((x+y)/2).reshape(-1, 1); print('Gaussian mixture BIC (1 vs 2 components): ' + ', '.join(f"{GM(k, random_state=0, n_init=10).fit(M).bic(M):.1f}" for k in (1, 2)))
    for u, lab in [(1, 'uM'), (1e-3, 'nM')]:
        print(f'OLS intercept with affinities in {lab}: {stats.linregress(np.log(d.Ki/1000/u), np.log(d.IC50/1000/u))[1]:.2f}')
LIND_PAIRS = [('Kd_CLK1', 'R_CLK1'), ('Kd_CLK1', 'P_CLK1'), ('R_CLK1', 'P_CLK1'),
              ('Kd_D1A', 'R_D1A'), ('Kd_D1A', 'P_D1A'), ('R_D1A', 'P_D1A'), ('Kd_GSK3B', 'R_GSK3B')]
CAP = {'Kd_D1A': 1, 'Kd_CLK1': 1, 'Kd_GSK3B': 10, 'R_D1A': 10, 'R_CLK1': 10, 'R_GSK3B': 10, 'P_D1A': 10, 'P_CLK1': 10}   # censoring limits (uM)
RANGE = {'Kd': (0.001, 1), 'R': (0.03, 10), 'P_D1A': (0.001, 10), 'P_CLK1': (0.005, 10)}   # tested range / enzyme-limited floor (uM)
def lindberg_pair(d, xc, yc):
    x, y = d[xc].values.astype(float), d[yc].values.astype(float)
    if xc.startswith('Kd'): x = x/1000
    m = ~(np.isnan(x) | np.isnan(y)); cens = m & ((x >= CAP[xc]) | (y >= CAP[yc])); m &= ~cens
    rx, ry = RANGE.get(xc, RANGE.get(xc.split('_')[0])), RANGE.get(yc, RANGE.get(yc.split('_')[0]))
    inr = (x >= rx[0]) & (x <= rx[1]) & (y >= ry[0]) & (y <= ry[1])
    return np.log(x[m]), np.log(y[m]), inr[m], int(cens.sum()), d.compound.values[m]
def lindberg(path='lindberg.csv'):
    d = pd.read_csv(path); rng = np.random.default_rng(2027)
    for xc, yc in LIND_PAIRS:
        x, y, inr, nc, _ = lindberg_pair(d, xc, yc)
        both_ic = xc[0] != 'K'
        noise = {'within-assay': (S_REP, S_REP), 'inter-lab': (S_IC, S_IC) if both_ic else (S_KI, S_IC)}
        split = {'lower': x < -5.6, 'upper': x >= -5.6} if xc == 'Kd_CLK1' else None
        print(f"\n==== {xc} vs {yc}: {nc} censored values excluded; {int((~inr).sum())} of {len(x)} usable values outside the tested range")
        report(analyse(x, y, f'Lindberg {xc} vs {yc} (all usable)', 1.0, noise, rng, split))
        if inr.sum() >= 7: report(analyse(x[inr], y[inr], f'Lindberg {xc} vs {yc} (within tested range)', 1.0, noise, rng))
    x, y = np.log(d.Kd_CLK1/1000), np.log(d.R_CLK1)
    print(f"\nFigure 3 as originally plotted (all {len(x)} compounds, censored values included): R2 = {r2(x.values, y.values):.2f}")
# ---------------------------------------------------------------- all-pairs log ratio-of-ratios (Tables S3 and S5)
def pair_r(x, y, idx=None):
    """r = ln[(x_j/x_i)/(y_j/y_i)] for all compound pairs, oriented so that x_j > x_i; pairs of the same compound skipped."""
    idx = np.arange(len(x)) if idx is None else idx
    i, j = np.triu_indices(len(x), 1); ok = idx[i] != idx[j]; i, j = i[ok], j[ok]
    dx, dy = x[j]-x[i], y[j]-y[i]; sg = np.where(dx >= 0, 1, -1)
    return sg*(dx-dy)
def ratio_stats(r):
    a = np.abs(r)
    return dict(med=np.median(r), mean=r.mean(), sd=r.std(ddof=1), f2=(a > np.log(2)).mean(), f10=(a > np.log(10)).mean())
def ratios(vpath='vegfr.csv', lpath='lindberg.csv'):
    rng = np.random.default_rng(2029); L = pd.read_csv(lpath); V = pd.read_csv(vpath); sets = []
    for lab, xc, yc in [('CLK1 KD vs radiometric', 'Kd_CLK1', 'R_CLK1'), ('DYRK1A KD vs radiometric', 'Kd_D1A', 'R_D1A'),
                        ('GSK-3b KD vs radiometric', 'Kd_GSK3B', 'R_GSK3B'), ('CLK1 radiometric vs ADP-Glo', 'R_CLK1', 'P_CLK1'),
                        ('DYRK1A radiometric vs ADP-Glo', 'R_D1A', 'P_D1A')]:
        x, y, *_ = lindberg_pair(L, xc, yc); sets.append((lab, x, y, (S_IC, S_IC) if xc[0] == 'R' else (S_KI, S_IC)))
    sets.insert(3, ('VEGFR KI vs IC50', np.log(V.Ki.values/1000), np.log(V.IC50.values/1000), (S_KI, S_IC)))
    for lab, x, y, (sx, sy) in sets:
        n = len(x); o = ratio_stats(pair_r(x, y)); B = []
        for _ in range(NB):
            k = rng.integers(0, n, n)
            if len(set(k)) > 3: B.append(ratio_stats(pair_r(x[k], y[k], k)))
        ci = {q: np.percentile([b[q] for b in B], [2.5, 97.5]) for q in o}
        nul = {nm: [ratio_stats(pair_r(a, b)) for a, b in null_sets(x, y, s1, s2, rng, 1000)] for nm, (s1, s2) in [('within-assay', (S_REP, S_REP)), ('inter-lab', (sx, sy))]}
        print(f"\n{lab}: n = {n}, pairs = {n*(n-1)//2}")
        print(f"  median ratio-of-ratios {np.exp(o['med']):.2f}-fold [{np.exp(ci['med'][0]):.2f}, {np.exp(ci['med'][1]):.2f}]; mean r {o['mean']:.2f} [{ci['mean'][0]:.2f}, {ci['mean'][1]:.2f}]")
        print(f"  SD r {o['sd']:.2f} [{ci['sd'][0]:.2f}, {ci['sd'][1]:.2f}] (sqrt2 x SD of D = {np.sqrt(2)*(y-x).std(ddof=1):.2f}); "
              f"pairs >2-fold {o['f2']:.2f} [{ci['f2'][0]:.2f}, {ci['f2'][1]:.2f}]; >10-fold {o['f10']:.2f} [{ci['f10'][0]:.2f}, {ci['f10'][1]:.2f}]")
        for nm, v in nul.items():
            print(f"  null {nm}: SD r {q95([w['sd'] for w in v])[0]:.2f}-{q95([w['sd'] for w in v])[2]:.2f}; >2-fold {q95([w['f2'] for w in v])[0]:.2f}-{q95([w['f2'] for w in v])[2]:.2f}; "
                  f">10-fold {q95([w['f10'] for w in v])[0]:.2f}-{q95([w['f10'] for w in v])[2]:.2f}")
# ---------------------------------------------------------------- additional checks (rank correlations, cluster rule, IC50 averaging)
def rank_line(x, y, k, rng, noise):
    from scipy.stats import spearmanr, kendalltau
    bs = [(spearmanr(x[k][b], y[k][b])[0], kendalltau(x[k][b], y[k][b])[0]) for b in (rng.integers(0, k.sum(), k.sum()) for _ in range(NB)) if len(set(b)) > 3]
    s = f"n = {k.sum()}, R2 = {r2(x[k], y[k]):.2f}, Spearman {spearmanr(x[k], y[k])[0]:.2f} [{np.nanpercentile([b[0] for b in bs], 2.5):.2f}, {np.nanpercentile([b[0] for b in bs], 97.5):.2f}], " \
        f"Kendall {kendalltau(x[k], y[k])[0]:.2f} [{np.nanpercentile([b[1] for b in bs], 2.5):.2f}, {np.nanpercentile([b[1] for b in bs], 97.5):.2f}]"
    for nm, (sx, sy) in noise.items():
        v = q95([spearmanr(a[k], b[k])[0] for a, b in null_sets(x, y, sx, sy, rng, 1000)]); s += f"; null {nm} Spearman {v[1]:.2f} [{v[0]:.2f}, {v[2]:.2f}]"
    return s
def extras(vpath='vegfr.csv', lpath='lindberg.csv'):
    from sklearn.cluster import KMeans
    from sklearn.mixture import GaussianMixture as GM
    rng = np.random.default_rng(2028)
    d = pd.read_csv(vpath); x, y = np.log(d.Ki.values/1000), np.log(d.IC50.values/1000)
    noise = {'within-assay': (S_REP, S_REP), 'inter-lab': (S_KI, S_IC)}
    line = y > 3 + (-10/10.6)*(x + 9); onl = (d.compound == 'OSI-632').values
    km = KMeans(2, n_init=20, random_state=0).fit_predict(np.c_[x, y]); up = km == km[np.argmax(line)]
    print(f"VEGFR: k-means (k = 2) on (ln KI, ln IC50) reproduces the visual clusters for {int((up == line)[~onl].sum())} of {int((~onl).sum())} compounds; OSI-632 assigned to the {'upper' if up[onl][0] else 'lower'} cluster")
    for lab, k in [('upper cluster', line & ~onl), ('lower cluster', ~line & ~onl), ('lower cluster incl. OSI-632 (k-means)', ~line),
                   ('IC50 reported once', (d.stdev == 0).values), ('IC50 replicated', (d.stdev > 0).values)]:
        print(f"  {lab}: " + rank_line(x, y, k, rng, noise))
    rep = (d.stdev > 0).values; gm = d.IC50.values/np.sqrt(1 + (d.stdev.values/d.IC50.values)**2)   # log-normal geometric mean from arithmetic mean and SD
    yg = np.where(rep, np.log(gm/1000), y); Dg = yg-x
    print(f"VEGFR with geometric instead of arithmetic means for the {rep.sum()} replicated IC50 values: R2 {r2(x, yg):.2f}, Deming {deming(x, yg, (S_IC/S_KI)**2):.2f}, "
          f"bias {Dg.mean():.2f} ({np.exp(Dg.mean()):.1f}-fold), SD {Dg.std(ddof=1):.2f}, CCC {ccc(x, yg):.2f}, Spearman {stats.spearmanr(x, yg)[0]:.2f}, n(IC50<KI) {(Dg < 0).sum()}; "
          f"largest shift {np.max(y-yg):.2f} ln units")
    L = pd.read_csv(lpath)
    for xc, yc in [('Kd_CLK1', 'R_CLK1'), ('Kd_CLK1', 'P_CLK1')]:
        x, y, *_ = lindberg_pair(L, xc, yc); vis = x < -5.6
        km = KMeans(2, n_init=20, random_state=0).fit_predict(np.c_[x, y]); lo = km == km[np.argmin(x)]
        bic = [GM(c, random_state=0, n_init=10).fit(np.c_[x, y]).bic(np.c_[x, y]) for c in (1, 2)]
        print(f"{xc} vs {yc}: k-means (k = 2) reproduces the visual split for {int((lo == vis).sum())} of {len(x)} compounds; 2D Gaussian mixture BIC 1 vs 2 components: {bic[0]:.1f}, {bic[1]:.1f}")
        nz = {'within-assay': (S_REP, S_REP), 'inter-lab': (S_KI, S_IC)}
        for lab, k in [('lower', vis), ('upper', ~vis)]: print(f"  {lab}: " + rank_line(x, y, k, rng, nz))
if __name__ == '__main__':
    {'vegfr': vegfr, 'lindberg': lindberg, 'ratios': ratios, 'extras': extras}.get(sys.argv[1] if len(sys.argv) > 1 else 'vegfr', vegfr)()
