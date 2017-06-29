# codes to work with the NGC 5548 photoionization grids supplied by Otho Adam Ulrich
import math as m
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
from itertools import cycle
from astropy.io import ascii

plot_fontsize = 20
plot_ticksize = 19
plot_lwidth = 2.5
xtick_padding = 15
tickwidth=2.5
markersize=14

#import seaborn as sns

lightday = 2.59e15 # cm - lazy/evil duplication...
    
def p_cf(r_ld,y,title='something vs r(lightdays)',vlines='none',
         log='yes',axes='none',hlines='none',label='a label',legends=['1','2','3','4'],
         xtitle='$R$ [lightdays]',ytitle='y',xlog='no',label_xy=[0.3,0.3]):
    """plot something vs radius in lightdays"""
    #sns.set_palette(sns.color_palette("hls", 20))
    markersymbols = ['b','r--','k:','m-.','co','gv']
    markers = cycle(markersymbols)
    l = iter(legends)
    try:
        for yi in y:
            plt.plot(r_ld,yi,next(markers),label=next(l),lw=2.5)
        if log=='yes':
            plt.yscale('log')
        if xlog=='yes':
            plt.xscale('log')
        plt.legend(fontsize=plot_fontsize-5,loc='lower right')
    except:
        markers = cycle(markersymbols)
        plt.plot(r_ld,y,next(markers),lw=2.5)
        if log=='yes':
            plt.yscale('log')
        if xlog=='yes':
            plt.xscale('log')
    #plt.title(title)
    if hlines != 'none':
        for hline in hlines:
            plt.axhline(y=hline,color='black')
    if vlines != 'none':
        for vline in vlines:
            plt.axvline(x=vline,color='black')
    if axes != 'none':
        plt.axis(axes)
    plt.xlabel(xtitle,fontsize=plot_fontsize)
    plt.ylabel(ytitle,rotation=90,fontsize=plot_fontsize)
    plt.figtext(label_xy[0],label_xy[1],label,backgroundcolor='white',fontsize=plot_fontsize)
    plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
    plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
    plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
    plt.tight_layout()
    plt.savefig(title+'.eps')
    plt.close()

data1 = np.array(ascii.read('lya_s0_cumdC.dat',data_start=2))
data2 = np.array(ascii.read('lya_s2_cumdC.dat',data_start=2))
r1 = [d[0] for d in data1]
r2 = [d[0] for d in data2]
c1 = [d[1]/(4*m.pi) for d in data1]
c2 = [d[1]/(4*m.pi) for d in data2]

p_cf(r1,[c1,c2],title='cf',log='no',legends=['s=0','s=2'],label='',ytitle=r'BLR coverage / 4$\pi$ steradian',axes=[0,140,0,1.1])
