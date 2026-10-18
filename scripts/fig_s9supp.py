"""fig_s9supp.py - reproduces the four figures of Supplement S9: the agreement
between K_I and IC50 for the current VEGFR1 and VEGFR2 extractions, and the
distribution of ln(IC50/K_I) for each.  Panels, null models, colours, bin edges
and random seed are those of figures.py (Figure S17), so that the two targets and
the main-text dataset are presented identically.  Requires affinity_agreement.py
in the same folder.

Usage:  python fig_s9supp.py vegfr1.csv vegfr2.csv
        ->  FigS19_vegfr1_agreement.png   FigS20_vegfr1_lnratio.png
            FigS21_vegfr2_agreement.png   FigS22_vegfr2_lnratio.png
"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.ticker import MaxNLocator
import affinity_agreement as A
B, G, O, T, MUT = '#2a78d6', '#52514e', '#eb6834', '#1baf7a', '#8a8985'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': G})
p1, p2 = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else ('vegfr1.csv', 'vegfr2.csv')
WID = [2, 3, 4, 5, 6, 7]; EDG = (0, 1, 2, 3, 4, 5, 7, 30); CTR = [0.5, 1.5, 2.5, 3.5, 4.5, 6, 8]
NSIM, NBOOT, BINW = 300, 1000, 0.4      # null datasets per model, bootstrap resamples, histogram bin width


def agreement(x, y, out, rng):
    """Four-panel agreement figure, identical in construction to Figure S17."""
    lam = (A.S_IC / A.S_KI) ** 2

    def curves(sx, sy):
        W, I = [], []
        for xs, ys in A.null_sets(x, y, sx, sy, rng, NSIM):
            W.append([A.window_r2(xs, ys, w) for w in WID]); I.append(A.inversion_bins(xs, ys, EDG))
        return np.array(W), np.array(I)
    Wi, Ii = curves(A.S_KI, A.S_IC); Wa, Ia = curves(A.S_REP, A.S_REP)
    obs_w = [A.window_r2(x, y, w) for w in WID]; obs_i = A.inversion_bins(x, y, EDG)
    bi = np.array([A.inversion_bins(x[k], y[k], EDG) for k in (rng.integers(0, len(x), len(x)) for _ in range(NBOOT))])
    ci = np.nanpercentile(bi, [2.5, 97.5], 0)
    fig, ax = plt.subplots(2, 2, figsize=(9, 8))
    a = ax[0, 0]; a.scatter(x, y, s=36, color=B, edgecolor='white', lw=1, zorder=3)
    xx = np.linspace(x.min() - 0.5, x.max() + 0.5, 10)
    s, i = stats.linregress(x, y)[:2]; bd = A.deming(x, y, lam); ad = y.mean() - bd * x.mean()
    a.plot(xx, xx, color=MUT, lw=1, ls=':', label='identity'); a.plot(xx, s * xx + i, color=G, lw=2, label=f'OLS, slope {s:.2f}')
    a.plot(xx, bd * xx + ad, color=O, lw=2, label=f'Deming (λ = {lam:.2f}), slope {bd:.2f}')
    a.set_xlabel('ln(K$_I$/µM)'); a.set_ylabel('ln(IC$_{50}$/µM)'); a.legend(frameon=False, fontsize=8)
    a.set_title('a  Errors-in-variables regression', loc='left', fontweight='bold')
    a = ax[0, 1]; M = (x + y) / 2; Dv = y - x; md, sd = Dv.mean(), Dv.std(ddof=1)
    a.scatter(M, Dv, s=36, color=B, edgecolor='white', lw=1, zorder=3)
    for v, l, c in [(md, f'bias {md:.2f} ({np.exp(md):.1f}-fold)', O),
                    (md + 1.96 * sd, f'+1.96 SD ({np.exp(md + 1.96 * sd):.0f}-fold)', G),
                    (md - 1.96 * sd, f'−1.96 SD ({np.exp(md - 1.96 * sd):.2f}-fold)', G)]:
        a.axhline(v, color=c, lw=1.5 if c == O else 1, ls='-' if c == O else '--')
        a.text(M.max() + 0.2, v + 0.1, l, va='bottom', fontsize=7.5, bbox=dict(fc='white', ec='none', pad=0.5))
    a.axhline(0, color=MUT, lw=1, ls=':'); a.set_xlim(None, M.max() + 3.6)
    a.set_xlabel('mean of ln(K$_I$/µM) and ln(IC$_{50}$/µM)'); a.set_ylabel('ln(IC$_{50}$/K$_I$)')
    a.set_title('b  Bland–Altman', loc='left', fontweight='bold')
    a = ax[1, 0]
    for W, c, l in [(Wi, O, 'null, inter-laboratory reproducibility'), (Wa, T, 'null, within-assay repeatability')]:
        lo, me, hi = np.nanpercentile(W, [2.5, 50, 97.5], 0)
        a.fill_between(WID, lo, hi, color=c, alpha=0.18, lw=0); a.plot(WID, me, color=c, lw=2, label=l)
    a.plot(WID, obs_w, 'o-', color=B, lw=2, ms=6, label='observed', zorder=3); a.set_ylim(0, 1)
    a.set_xlabel('window width (ln units of K$_I$)'); a.set_ylabel('mean R$^2$ within sliding window')
    a.legend(frameon=False, fontsize=8, loc='upper left'); a.set_title('c  Local correlation vs. null models', loc='left', fontweight='bold')
    a = ax[1, 1]
    for I, c, l in [(Ii, O, 'null, inter-laboratory reproducibility'), (Ia, T, 'null, within-assay repeatability')]:
        lo, me, hi = np.nanpercentile(I, [2.5, 50, 97.5], 0)
        a.fill_between(CTR, lo, hi, color=c, alpha=0.18, lw=0); a.plot(CTR, me, color=c, lw=2, label=l)
    a.errorbar(CTR, obs_i, yerr=[obs_i - ci[0], ci[1] - obs_i], fmt='o-', color=B, lw=2, ms=6, capsize=3,
               label='observed (bootstrap 95% CI)', zorder=3)
    a.axhline(0.05, color=MUT, lw=1, ls=':'); a.text(7.6, 0.065, '5%', fontsize=7.5, color=G)
    a.set_xlabel('|Δ ln K$_I$| between compound pairs'); a.set_ylabel('fraction of pairs with reversed IC$_{50}$ order')
    a.legend(frameon=False, fontsize=8); a.set_title('d  Rank inversions vs. separation', loc='left', fontweight='bold')
    fig.tight_layout(); fig.savefig(out, dpi=300); plt.close(fig)
    print(f'{out}: n = {len(x)}, OLS {s:.2f}, Deming {bd:.2f}, bias {md:.2f} ({np.exp(md):.1f}-fold), SD {sd:.2f}')


def lnratio(D, out):
    """Histogram of ln(IC50/K_I), as in Figure S14."""
    lo = np.floor(D.min() / BINW) * BINW; hi = np.ceil(D.max() / BINW) * BINW + BINW / 2
    fig, ax = plt.subplots(figsize=(5, 3.2))
    n, _, _ = ax.hist(D, bins=np.arange(lo, hi, BINW), color=B, edgecolor='white', linewidth=2)
    ax.axvline(0, color=G, ls=':', lw=1); ax.axvline(D.mean(), color=O, lw=2)
    top = n.max() * 1.25; ax.set_ylim(0, top)
    ax.text(lo, top * 0.95, f'mean = {D.mean():.2f} (IC$_{{50}}$/K$_I$ = {np.exp(D.mean()):.1f})\n'
                            f'median = {np.median(D):.2f}\nn = {len(D)}; IC$_{{50}}$ < K$_I$ for {(D < 0).sum()}',
            fontsize=8, va='top')
    ax.set_xlabel('ln(IC$_{50}$/K$_I$)'); ax.set_ylabel('number of compounds')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    fig.tight_layout(); fig.savefig(out, dpi=300); plt.close(fig)
    print(f'{out}: n = {len(D)}, mean {D.mean():.2f}, median {np.median(D):.2f}, SD {D.std(ddof=1):.2f}, '
          f'IC50 < KI for {(D < 0).sum()}')


for path, ag, lr in ((p1, 'FigS19_vegfr1_agreement.png', 'FigS20_vegfr1_lnratio.png'),
                     (p2, 'FigS21_vegfr2_agreement.png', 'FigS22_vegfr2_lnratio.png')):
    d = pd.read_csv(path)
    x, y = np.log(d.Ki.values / 1000), np.log(d.IC50.values / 1000)
    agreement(x, y, ag, np.random.default_rng(11))
    lnratio(np.log(d.IC50.values / d.Ki.values), lr)
