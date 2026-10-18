"""figures.py - reproduces Figures S15, S18 and S19 (requires affinity_agreement.py in the same folder).
Usage: python figures.py vegfr.csv lindberg.csv"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import stats
import affinity_agreement as A
B, G, O, T, MUT = '#2a78d6', '#52514e', '#eb6834', '#1baf7a', '#8a8985'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': G})
vpath, lpath = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else ('vegfr.csv', 'lindberg.csv')
rng = np.random.default_rng(11)
# ------------------------------------------------ Figure S14
d = pd.read_csv(vpath); D = np.log(d.IC50/d.Ki)
fig, ax = plt.subplots(figsize=(5, 3.2)); ax.hist(D, bins=np.arange(-3.25, 3.76, 0.5), color=B, edgecolor='white', linewidth=2)
ax.axvline(0, color=G, ls=':', lw=1); ax.axvline(D.mean(), color=O, lw=2)
ax.text(-3.2, 5.9, f'mean = {D.mean():.2f} (IC$_{{50}}$/K$_I$ = {np.exp(D.mean()):.1f})\nmedian = {np.median(D):.2f}\nn = {len(D)}; IC$_{{50}}$ < K$_I$ for {(D < 0).sum()}', fontsize=8, va='top')
ax.set_xlabel('ln(IC$_{50}$/K$_I$)'); ax.set_ylabel('number of compounds'); ax.set_ylim(0, 6.5); fig.tight_layout(); fig.savefig('FigS14.png', dpi=300)
# ------------------------------------------------ Figure S17 (VEGFR)
x, y = np.log(d.Ki.values/1000), np.log(d.IC50.values/1000); lam = (A.S_IC/A.S_KI)**2
WID = [2, 3, 4, 5, 6, 7]; EDG = (0, 1, 2, 3, 4, 5, 7, 30); CTR = [0.5, 1.5, 2.5, 3.5, 4.5, 6, 8]
def curves(sx, sy, n=300):
    W, I = [], []
    for k, (xs, ys) in enumerate(A.null_sets(x, y, sx, sy, rng, n)):
        W.append([A.window_r2(xs, ys, w) for w in WID]); I.append(A.inversion_bins(xs, ys, EDG))
    return np.array(W), np.array(I)
Wi, Ii = curves(A.S_KI, A.S_IC); Wa, Ia = curves(A.S_REP, A.S_REP)
obs_w = [A.window_r2(x, y, w) for w in WID]; obs_i = A.inversion_bins(x, y, EDG)
bi = np.array([A.inversion_bins(x[k], y[k], EDG) for k in (rng.integers(0, len(x), len(x)) for _ in range(1000))]); ci = np.nanpercentile(bi, [2.5, 97.5], 0)
fig, ax = plt.subplots(2, 2, figsize=(9, 8))
a = ax[0, 0]; a.scatter(x, y, s=36, color=B, edgecolor='white', lw=1, zorder=3); xx = np.linspace(x.min()-0.5, x.max()+0.5, 10)
s, i = stats.linregress(x, y)[:2]; bd = A.deming(x, y, lam); ad = y.mean()-bd*x.mean()
a.plot(xx, xx, color=MUT, lw=1, ls=':', label='identity'); a.plot(xx, s*xx+i, color=G, lw=2, label=f'OLS, slope {s:.2f}')
a.plot(xx, bd*xx+ad, color=O, lw=2, label=f'Deming (λ = {lam:.2f}), slope {bd:.2f}')
a.set_xlabel('ln(K$_I$/µM)'); a.set_ylabel('ln(IC$_{50}$/µM)'); a.legend(frameon=False, fontsize=8); a.set_title('a  Errors-in-variables regression', loc='left', fontweight='bold')
a = ax[0, 1]; M = (x+y)/2; Dv = y-x; md, sd = Dv.mean(), Dv.std(ddof=1); a.scatter(M, Dv, s=36, color=B, edgecolor='white', lw=1, zorder=3)
for v, l, c in [(md, f'bias {md:.2f} ({np.exp(md):.1f}-fold)', O), (md+1.96*sd, f'+1.96 SD ({np.exp(md+1.96*sd):.0f}-fold)', G), (md-1.96*sd, f'−1.96 SD ({np.exp(md-1.96*sd):.2f}-fold)', G)]:
    a.axhline(v, color=c, lw=1.5 if c == O else 1, ls='-' if c == O else '--'); a.text(M.max()+0.2, v+0.1, l, va='bottom', fontsize=7.5, bbox=dict(fc='white', ec='none', pad=0.5))
a.axhline(0, color=MUT, lw=1, ls=':'); a.set_xlim(None, M.max()+3.6)
a.set_xlabel('mean of ln(K$_I$/µM) and ln(IC$_{50}$/µM)'); a.set_ylabel('ln(IC$_{50}$/K$_I$)'); a.set_title('b  Bland–Altman', loc='left', fontweight='bold')
a = ax[1, 0]
for W, c, l in [(Wi, O, 'null, inter-laboratory reproducibility'), (Wa, T, 'null, within-assay repeatability')]:
    lo, me, hi = np.nanpercentile(W, [2.5, 50, 97.5], 0); a.fill_between(WID, lo, hi, color=c, alpha=0.18, lw=0); a.plot(WID, me, color=c, lw=2, label=l)
a.plot(WID, obs_w, 'o-', color=B, lw=2, ms=6, label='observed', zorder=3); a.set_ylim(0, 1)
a.set_xlabel('window width (ln units of K$_I$)'); a.set_ylabel('mean R$^2$ within sliding window'); a.legend(frameon=False, fontsize=8, loc='upper left'); a.set_title('c  Local correlation vs. null models', loc='left', fontweight='bold')
a = ax[1, 1]
for I, c, l in [(Ii, O, 'null, inter-laboratory reproducibility'), (Ia, T, 'null, within-assay repeatability')]:
    lo, me, hi = np.nanpercentile(I, [2.5, 50, 97.5], 0); a.fill_between(CTR, lo, hi, color=c, alpha=0.18, lw=0); a.plot(CTR, me, color=c, lw=2, label=l)
a.errorbar(CTR, obs_i, yerr=[obs_i-ci[0], ci[1]-obs_i], fmt='o-', color=B, lw=2, ms=6, capsize=3, label='observed (bootstrap 95% CI)', zorder=3)
a.axhline(0.05, color=MUT, lw=1, ls=':'); a.text(7.6, 0.065, '5%', fontsize=7.5, color=G)
a.set_xlabel('|Δ ln K$_I$| between compound pairs'); a.set_ylabel('fraction of pairs with reversed IC$_{50}$ order'); a.legend(frameon=False, fontsize=8); a.set_title('d  Rank inversions vs. separation', loc='left', fontweight='bold')
fig.tight_layout(); fig.savefig('FigS17.png', dpi=300)
# ------------------------------------------------ Figure S18 (Lindberg)
L = pd.read_csv(lpath)
NAME = {'Kd': 'K$_D$', 'R': 'IC$_{50}$ radiometric', 'P': 'IC$_{50}$ ADP-Glo'}
ROWS, COLS = 4, 2                                       # Figure S18 layout: 4 rows x 2 columns
fig, ax = plt.subplots(ROWS, COLS, figsize=(4.6*COLS, 2.9*ROWS)); ax = ax.ravel(); sds = []
for k, (xc, yc) in enumerate(A.LIND_PAIRS):
    xx_, yy_, inr, nc, _ = A.lindberg_pair(L, xc, yc); a = ax[k]
    a.scatter(xx_[inr], yy_[inr], s=26, color=B, edgecolor='white', lw=0.8, zorder=3, label='within tested range')
    a.scatter(xx_[~inr], yy_[~inr], s=26, facecolor='white', edgecolor=B, lw=1.2, zorder=3, label='outside tested range')
    g = np.linspace(min(xx_.min(), yy_.min())-0.5, max(xx_.max(), yy_.max())+0.5, 10)
    bd = A.deming(xx_, yy_, 1.0); a.plot(g, bd*g+yy_.mean()-bd*xx_.mean(), color=O, lw=1.8, label='Deming fit')
    a.plot(g, g+(yy_-xx_).mean(), color=G, lw=1, ls='--', label='slope 1 + bias'); a.plot(g, g, color=MUT, lw=1, ls=':', label='identity')
    kin = xc.split('_')[1].replace('D1A', 'DYRK1A').replace('GSK3B', 'GSK-3β')
    a.set_title(f'{kin}: {NAME[yc.split("_")[0]]} vs {NAME[xc.split("_")[0]]}', loc='left', fontsize=8.5, fontweight='bold')
    a.set_xlabel(f'ln({NAME[xc.split("_")[0]]}/µM)', fontsize=8); a.set_ylabel(f'ln({NAME[yc.split("_")[0]]}/µM)', fontsize=8)
    a.text(0.03, 0.97, f'n = {len(xx_)} ({nc} censored excl.)\nR$^2$ = {A.r2(xx_, yy_):.2f}; Deming {bd:.2f}\nbias {np.exp((yy_-xx_).mean()):.1f}-fold; SD {np.std(yy_-xx_, ddof=1):.2f}', transform=a.transAxes, va='top', fontsize=7.5)
    sds.append(np.std(yy_-xx_, ddof=1))
    if k == 0: a.legend(frameon=False, fontsize=6.5, loc='lower right')
a = ax[7]; labs = ['CLK1 K$_D$–R', 'CLK1 K$_D$–A', 'CLK1 R–A', 'DYRK1A K$_D$–R', 'DYRK1A K$_D$–A', 'DYRK1A R–A', 'GSK-3β K$_D$–R']
yy = np.arange(len(sds))[::-1]; a.axvspan(0.41, 0.92, color=T, alpha=0.18, lw=0, label='within-assay null (95%)'); a.axvspan(1.3, 2.9, color=O, alpha=0.15, lw=0, label='inter-laboratory null (95%)')
a.scatter(sds, yy, color=B, s=36, zorder=3, label='observed'); a.set_yticks(yy); a.set_yticklabels(labs, fontsize=7.5); a.set_xlabel('SD of ln(y/x)'); a.set_xlim(0, 3.2)
a.set_title('Scatter vs. noise scales', loc='left', fontsize=8.5, fontweight='bold'); a.legend(frameon=False, fontsize=6.5, loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=2)
fig.tight_layout(); fig.savefig('FigS18.png', dpi=300)
print('figures written; SDs', np.round(sds, 2))
