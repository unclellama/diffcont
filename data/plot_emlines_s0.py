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

def load_line(line):
    linefile = 'lineresults_s0/s0_linedata_'+line+'.js'
    with open(linefile) as json_data:
        linedata = json.load(json_data)
    data = []
    for Ncol in linedata['Ncol_arr']:
        data.append(linedata['data_Ncol_'+Ncol])
    return [linedata['Ncol_arr'],data]

def plot_lineratios(line1,line2,nHprint='none'):
    Ncol_arr,data_line1 = load_line(line1)
    Ncol_arr,data_line2 = load_line(line2)
    nH_arr = (data_line1[0])['nHarr']
    print(nH_arr)
    if nHprint == 'none':
        nHprint = [str(n) for n in nH_arr]
    try:
        L_line1_obs,eL_line1_obs = observed_lumin(line1)
        L_line2_obs,eL_line2_obs = observed_lumin(line2)
        hasData = 1
    except:
        hasData = 0
    L_line1_Ncol = []
    L_line2_Ncol = []
    for i in range(len(Ncol_arr)):
        L_line1_Ncol.append((data_line1[i])['L'])
        L_line2_Ncol.append((data_line2[i])['L'])

    for i in range(len(Ncol_arr)):
        plt.plot(L_line1_Ncol[i],L_line2_Ncol[i],
                label=r'$\log(N_{\mathrm{col}})=$'+str(Ncol_arr[i]))
        for j in range(len(nH_arr)):
            plt.text((L_line1_Ncol[0])[j],(L_line2_Ncol[0])[j],nHprint[j],
                     fontsize=plot_fontsize-5)
    plt.xlabel(r'L('+typeline(line1)+') (erg s$^{-1}$)',fontsize=plot_fontsize)
    plt.ylabel(r'L('+typeline(line2)+') (erg s$^{-1}$)',rotation=90,fontsize=plot_fontsize)
    plt.yscale('log')
    plt.xscale('log')
    plt.plot([1,1e50],[1,1e50],color='gray',linestyle='--')
    plt.axis([1e39,3e43,1e37,4e43])
    if hasData == 1:
        plt.plot([L_line1_obs,L_line1_obs],[L_line2_obs,L_line2_obs],
            color='black',lw=0,marker='*',markersize=20,label='Observed ratio')
    plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
    plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
    plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
    plt.tight_layout()
    plt.legend(loc='lower right',numpoints=1)
    plt.savefig('lineresults_s0/L_'+line1+'_'+line2+'.eps')
    plt.close()

plotstring_nH = ['     7','','','',r'     $\log(n_\mathrm{H})$=8','','','',
          '       9','','','','','','','','      11','','','',
          '','','','','','','','','14']
plot_lineratios('lya','civ',plotstring_nH)
plot_lineratios('lya','hbeta',plotstring_nH)
plot_lineratios('lya','he4686',plotstring_nH)
plot_lineratios('lya','he1640',plotstring_nH)
plot_lineratios('lya','mgii',plotstring_nH)
plot_lineratios('lya','halpha',plotstring_nH)
plot_lineratios('halpha','hbeta',plotstring_nH)
