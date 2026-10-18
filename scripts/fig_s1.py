"""fig_s1.py - Figure S1: sensitivity of the non-ideality correction factor
K_D^gamma = gamma_P gamma_L / gamma_PL to deviations of the activity coefficients
from unity (Supplementary S1).  Computed from the equations; no input data.

Usage:  python fig_s1.py        ->  FigS1.png
"""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10.5,'axes.spines.top':False,'axes.spines.right':False,'mathtext.default':'it'})
fig=plt.figure(figsize=(10.2,5.6),dpi=220)
ax=fig.add_axes([0.10,0.12,0.54,0.82])
g=np.logspace(-2,0,300)
cols={'lin':'#2f5fd0','quad':'#d9472b','cub':'#3a9a4a'}
ax.axhline(1,color='#222222',lw=1.4,ls='-',zorder=1)
ax.plot(g,g,color=cols['lin'],lw=2,label=r'linear, $K_D^\gamma=\gamma$')
ax.plot(g,g**2,color=cols['quad'],lw=2,label=r'quadratic, $K_D^\gamma=\gamma^2$')
ax.plot(g,g**3,color=cols['cub'],lw=2,label=r'cubic (worst case), $K_D^\gamma=\gamma^3$')
ax.plot([],[],color='#222222',lw=1.4,label=r'complete cancellation, $K_D^\gamma=1$')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.008,1.05); ax.set_ylim(1e-6,2)
for y,lab in [(0.1,'×10'),(0.02,'×50'),(0.01,'×100')]:
    ax.axhline(y,color='#e08a2c',lw=0.9,ls=(0,(3,3)),zorder=0); ax.text(0.0085,y*1.12,lab,va='bottom',fontsize=9,color='#8a5a1a')
pts=[('A',0.1,1),('B',0.02,1),('C',0.01,1),('D',0.1**0.5,2),('E',0.02**0.5,2),('F',0.01**0.5,2),('G',0.1**(1/3),3),('H',0.02**(1/3),3),('I',0.01**(1/3),3)]
for n,x,p in pts:
    ax.plot([x],[x**p],'o',ms=6.5,mfc='white',mec='#222222',mew=1.2,zorder=5)
    ax.annotate(n,(x,x**p),xytext=(4,-13) if n=='C' else (-11,4),textcoords='offset points',fontsize=10,fontweight='bold')
from matplotlib.ticker import FuncFormatter
import math
fmt=FuncFormatter(lambda v,p: ('%g'%v) if v>=0.001 else r'$10^{%d}$'%round(math.log10(v)))
ax.xaxis.set_major_formatter(fmt); ax.yaxis.set_major_formatter(fmt)
ax.set_xlabel(r'activity coefficient, $\gamma$'); ax.set_ylabel(r'$K_D^\gamma=\gamma_P\gamma_L/\gamma_{PL}$')
ax.legend(frameon=False,loc='lower right',fontsize=9)
# tables
t1=fig.add_axes([0.70,0.58,0.29,0.34]); t1.axis('off')
r1=[['linear',r'1 or $\gamma_{PL}$',r'$\gamma$',r'1 or $\gamma_P$',r'$\gamma$'],['quadratic',r'$\gamma$',r'$\gamma$','1',r'$\gamma^2$'],['cubic',r'$\gamma$',r'$\gamma$',r'$1/\gamma$',r'$\gamma^3$']]
tb=t1.table(cellText=r1,colLabels=['scenario',r'$\gamma_P$',r'$\gamma_L$',r'$\gamma_{PL}$',r'$K_D^\gamma$'],loc='center',cellLoc='center',edges='horizontal')
tb.auto_set_font_size(False); tb.set_fontsize(9.5); tb.scale(1,1.5)
t2=fig.add_axes([0.73,0.10,0.23,0.42]); t2.axis('off')
r2=[[n,f'{x**p:.2f}',f'{x:.2f}'] for n,x,p in pts]
tb2=t2.table(cellText=r2,colLabels=['point',r'$K_D^\gamma$',r'$\gamma$'],loc='center',cellLoc='center',edges='horizontal')
tb2.auto_set_font_size(False); tb2.set_fontsize(9.5); tb2.scale(1,1.25)
fig.savefig('FigS1.png')
