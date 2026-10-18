"""bsa_gamma.py - BSA activity coefficients from second virial coefficients (Supplementary S2, Figure S9).
ln(gamma) = 2*C2*m + ...,  C2 = rho_w*(B2 - nu*M),  B2 = B22*M**2   (Agena et al.; equations 11-12)
B22 (mL mol g^-2) from Ma et al. (CG-MALS), Table 1; standard error 1.8e-6 mL mol g^-2."""
import numpy as np
M, NU, RHO = 66000.0, 0.733, 0.997e-3      # g/mol, mL/g, kg/mL (as in BSAAgena.xlsx)
SE = 1.8e-6
B22 = {0.025: {3.76: 15.30, 4.64: 1.13, 4.89: 1.92, 5.11: 3.18, 5.47: 3.34, 6.39: 6.30, 7.37: 8.60},
       0.05:  {3.73: 10.80, 4.61: 1.62, 4.86: 2.16, 5.08: 3.62, 5.44: 4.25, 6.36: 5.19, 7.34: 6.93},
       0.075: {3.70: 9.94, 4.58: 1.63, 4.83: 3.11, 5.05: 3.02, 5.41: 3.60, 6.33: 4.78, 7.31: 6.86},
       0.15:  {3.62: 6.97, 4.50: 3.16, 4.75: 3.27, 4.97: 3.47, 5.33: 3.39, 6.25: 6.16, 7.23: 6.64}}   # x 1e-5
def C2(b22): return RHO*(b22*M**2 - NU*M)          # kg/mol
def lng(b22, c): return 2*C2(b22)*c                 # c in mol/L ~ mol/kg in the dilute limit
if __name__ == '__main__':
    print('NaCl  pH    B22(1e-5)  C2(kg/mol)  ln g(1 nM)  ln g(1 uM)  +-SE(1 uM)  ln g(10 uM)  ln g(100 uM)  ln g(1 mM)  c(5%) uM')
    for s, d in B22.items():
        for ph, b in d.items():
            b *= 1e-5; c5 = np.log(1.05)/(2*C2(b))*1e6 if C2(b) > 0 else np.inf
            print(f"{s:<5} {ph:<5} {b*1e5:6.2f}    {C2(b):8.1f}   {lng(b,1e-9):.1e}     {lng(b,1e-6):.2e}    {2*RHO*SE*M**2*1e-6:.1e}     {lng(b,1e-5):.4f}      {lng(b,1e-4):.3f}        {lng(b,1e-3):.3f}      {c5:.0f}")

def figures():
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    G = '#52514e'; plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': G})
    c = np.logspace(-9, -3, 200); cm = plt.get_cmap('viridis')
    fig, ax = plt.subplots(2, 2, figsize=(9, 7), sharey=True); ax = ax.ravel()
    for a, (s, d) in zip(ax, B22.items()):
        a.axvspan(1e-3, 1, color='#d9d8d4', alpha=0.6, lw=0, label='typical assay range (1 nM–1 µM)')
        for k, (ph, b) in enumerate(d.items()):
            b *= 1e-5; col = cm(k/(len(d)-1))
            a.fill_between(c*1e6, np.exp(lng(b-SE, c)), np.exp(lng(b+SE, c)), color=col, alpha=0.25, lw=0)
            a.plot(c*1e6, np.exp(lng(b, c)), color=col, lw=1.8, label=f'pH {ph}')
        a.axhline(1.05, color=G, ls=':', lw=1); a.set_xscale('log'); a.set_xlim(1e-3, 1e3)
        a.set_title(f'[NaCl] = {int(s*1000)} mM', loc='left', fontweight='bold'); a.set_xlabel('[BSA] (µM)'); a.set_ylabel('γ (virial, molar scale)')
    ax[0].legend(frameon=False, fontsize=7, loc='upper left'); fig.tight_layout(); fig.savefig('FigS9_corrected.png', dpi=300)
    fig, a = plt.subplots(figsize=(4.8, 3.6)); cols = ['#2a78d6', '#eb6834', '#1baf7a', '#c89211']
    a.axvspan(1e-3, 1, color='#d9d8d4', alpha=0.6, lw=0, label='typical assay range')
    for (s, d), col in zip(B22.items(), cols):
        ph, b = list(d.items())[-1]; b *= 1e-5
        a.plot(c*1e6, np.exp(lng(b, c)), color=col, lw=2, label=f'[NaCl] = {int(s*1000)} mM (pH {ph})')
    a.axhline(1.05, color=G, ls=':', lw=1); a.set_xscale('log'); a.set_xlim(1e-3, 1e3); a.set_xlabel('[BSA] (µM)'); a.set_ylabel('γ$^{*,V}$')
    a.legend(frameon=False, fontsize=7, loc='upper left'); fig.tight_layout(); fig.savefig('FigS3D_corrected.png', dpi=300)

if __name__ == '__main__' and len(__import__('sys').argv) > 1: figures()
