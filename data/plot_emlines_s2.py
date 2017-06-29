# =======================================================================================
# get BEL luminosities + weighted radii for a grid of nH, Ncol, for s=0 models
# - basically calls do_line() for each line, nH, Ncol combination.
# =======================================================================================

import math as m
import matplotlib.pyplot as plt
from matplotlib import ticker
from ewgrids import *
import json

plot_fontsize = 20
plot_ticksize = 19
plot_lwidth = 2.5
xtick_padding = 15
tickwidth=2.5
markersize=14

def load_line_s2(line):
    linefile = 'lineresults_s2/s2_linedata_'+line+'.js'
    with open(linefile) as json_data:
        linedata = json.load(json_data)
    data = []
    for logU in linedata['U_arr']:
        data.append(linedata['data_logU_'+logU])
    return [linedata['U_arr'],data]

def plot_linelumin(linelist):
    L_lya,eL_lya = observed_lumin('lya')
    luminblock = []
    markersymbols = ['b','r','k','m']
    markers = cycle(markersymbols)
    for i in range(len(linelist)):
        U_arr,data = load_line_s2(linelist[i])
        lumins = []
        for j in range(len(U_arr)):
            lumins.append((data[j])['L'])
        plt.plot(U_arr,lumins,next(markers),label=typeline(linelist[i]))
    L_lya,eL_lya = observed_lumin('lya')
    L_civ,eL_civ = observed_lumin('civ')
    L_hbeta,eL_hbeta = observed_lumin('hbeta')
    L_he4686,eL_he4686 = observed_lumin('he4686')
    plt.axhline(y=L_lya,color='blue',ls='--')
    plt.axhline(y=L_civ,color='red',ls='--')
    plt.axhline(y=L_hbeta,color='black',ls='--')
    plt.axhline(y=L_he4686,color='magenta',ls='--')
    plt.axvline(x=-1.23,color='black')
    plt.xlabel(r'$\log$(U)',fontsize=plot_fontsize)
    plt.ylabel(r'$L_{\mathrm{line}}$ (erg s$^{-1}$)',rotation=90,fontsize=plot_fontsize)
    plt.yscale('log')
    plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
    plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
    plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
    plt.tight_layout()
    plt.legend(loc='upper right',numpoints=1,fontsize=plot_fontsize)
    plt.savefig('lineresults_s2/Llines.eps')
    plt.close()

plot_linelumin(linelist=['lya','civ','hbeta','he4686'])

