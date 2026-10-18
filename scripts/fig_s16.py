"""fig_s16.py - Figure S16: dependence of the apparent K_I on K_M(ATP), computed
with the Cheng-Prusoff equation from the IC50 values and ATP concentrations of
Lindberg et al. (25 uM ATP) and Grygier et al. (100 uM ATP).  No input files.

Usage:  python fig_s16.py       ->  FigS16.png
"""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
S={'Grygier':100.0,'Lindberg':25.0}   # ATP, uM
ser=[('DYRK1A','harmine','Grygier',97,'#2f6bbf'),('DYRK1A','AZ191','Grygier',188,'#e0352b'),
     ('DYRK1A','harmine','Lindberg',61.7,'#4aa85a'),('DYRK1A','AZ191','Lindberg',19.1,'#5fb2e6'),
     ('DYRK1B','harmine','Grygier',115,'#8e3a99'),('DYRK1B','AZ191','Grygier',66,'#f2c13d'),
     ('DYRK1B','harmine','Lindberg',94,'#8a8a8a'),('DYRK1B','AZ191','Lindberg',8.3,'#b5562a')]
ki=lambda ic,s,km: ic/(1+s/km)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(2,2,figsize=(10,7.4),dpi=200); ax=ax.ravel()
km=np.logspace(-2,3,400)
for k,c,src,ic,col in ser:
    lab=f'{k}–{c} ({src})'
    ax[0].plot(km,ki(ic,S[src],km),color=col,lw=1.8,label=lab)
    ax[3].loglog(np.logspace(-4,1,200),ki(ic,S[src],np.logspace(-4,1,200)),color=col,lw=1.6)
    if k=='DYRK1A':
        x=np.linspace(30,120,200); ax[1].plot(x,ki(ic,S[src],x),color=col,lw=1.8)
    else:
        x=np.linspace(80,105,200); ax[2].plot(x,ki(ic,S[src],x),color=col,lw=1.8)
ax[0].set_xscale('log'); ax[0].set_ylim(0,105)
ax[0].set_title('A) Full range',loc='left',fontsize=10)
ax[1].set_title('B) DYRK1A, reported K$_M$ range',loc='left',fontsize=10)
ax[2].set_title('C) DYRK1B, reported K$_M$ range',loc='left',fontsize=10)
ax[3].set_title('D) Low K$_M$ (log scale)',loc='left',fontsize=10)
for a in ax:
    a.set_xlabel('K$_M$(ATP) (µM)'); a.set_ylabel('K$_I$ (nM)'); a.grid(color='#e6e6e6',lw=0.7); a.set_axisbelow(True)
ax[0].legend(frameon=False,fontsize=7.5,loc='upper left')
fig.tight_layout(); fig.savefig('FigS16.png')
