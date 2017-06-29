# =======================================================================================
# get BEL luminosities + weighted radii for a grid of nH, Ncol, for s=0 models
# - basically calls do_line() for each line, nH, Ncol combination.
# =======================================================================================

import math as m
import matplotlib.pyplot as plt
from matplotlib import ticker
from ewgrids import *
import json

plotNcol = 22.5

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

# plot lumins

Ncol_arr,data_lya = load_line('lya')
Ncol_arr,data_civ = load_line('civ')
Ncol_arr,data_halpha = load_line('halpha')
Ncol_arr,data_hbeta = load_line('hbeta')
Ncol_arr,data_he4686 = load_line('he4686')
Ncol_arr,data_he1640 = load_line('he1640')
floatNc = [float(N) for N in Ncol_arr]
print(floatNc)
from bisect import bisect_left
i = bisect_left(floatNc, plotNcol+0.01)
markersymbols = ['b','r','k','m']
markers = cycle(markersymbols)
plt.plot((data_lya[i])['nHarr'],(data_lya[i])['L'],next(markers),label=typeline('lya'),lw=2.5)
plt.plot((data_civ[i])['nHarr'],(data_civ[i])['L'],next(markers),label=typeline('civ'),lw=2.5)
plt.plot((data_hbeta[i])['nHarr'],(data_hbeta[i])['L'],next(markers),label=typeline('hbeta'),lw=2.5)
plt.plot((data_he4686[i])['nHarr'],(data_he4686[i])['L'],next(markers),label=typeline('he4686'),lw=2.5)
L_lya,eL_lya = observed_lumin('lya')
L_civ,eL_civ = observed_lumin('civ')
L_hbeta,eL_hbeta = observed_lumin('hbeta')
L_he4686,eL_he4686 = observed_lumin('he4686')
plt.axhline(y=L_lya,color='blue',ls='--',lw=2.5)
plt.axhline(y=L_civ,color='red',ls='--',lw=2.5)
plt.axhline(y=L_hbeta,color='black',ls='--',lw=2.5)
plt.axhline(y=L_he4686,color='magenta',ls='--',lw=2.5)
plt.xlabel(r'$\log(n_\mathrm{H}/1 \mathrm{cm}^{-3})$',fontsize=plot_fontsize)
plt.ylabel(r'$L_{\mathrm{line}}$ [erg cm$^{-1}$ s$^{-1}$]',fontsize=plot_fontsize)
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.figtext(0.2,0.25,r's=0, $\log(n_\mathrm{col})=\mathrm{'+str(plotNcol)+'}$',fontsize=plot_fontsize)
plt.yscale('log')
plt.legend(loc='lower right',numpoints=1,fontsize=plot_fontsize)
plt.axvline(x=10.75,color='black',lw=1.5)
plt.axis([7,14,5e37,5e43])
#plt.tight_layout()
plt.savefig('lineresults_s0/L_nH.eps',bbox_inches='tight')
plt.close()

markersymbols = ['b','r--','k:','m-.','co','gs']
markers = cycle(markersymbols)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_lya[i])['T_aniso']],next(markers),label=typeline('lya'),lw=2.5)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_civ[i])['T_aniso']],next(markers),label=typeline('civ'),lw=2.5)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_halpha[i])['T_aniso']],next(markers),label=typeline('halpha'),lw=2.5)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_hbeta[i])['T_aniso']],next(markers),label=typeline('hbeta'),lw=2.5)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_he4686[i])['T_aniso']],next(markers),label=typeline('he4686'),lw=2.5)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_he1640[i])['T_aniso']],next(markers),label=typeline('he1640'),lw=2.5)
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.ylabel(r'$\tau_\epsilon\,\mathrm{[days]}$',fontsize=plot_fontsize)
plt.xlabel(r'$\log(n_\mathrm{H}/1 \mathrm{cm}^{-3})$',fontsize=plot_fontsize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.figtext(0.2,0.25,r'$\log(n_\mathrm{col})=\mathrm{'+str(plotNcol)+'}$',fontsize=plot_fontsize)
plt.axis([7,14,0,180])
plt.axvline(x=10.75,color='black',lw=1.5)
plt.legend(loc='upper right',numpoints=1,fontsize=plot_fontsize)
plt.tight_layout()
plt.savefig('lineresults_s0/Tau_aniso_nH.eps')
plt.close()

markersymbols = ['b','r--','k:','m-.','co','gs']
markers = cycle(markersymbols)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_lya[i])['T_aniso_resp']],next(markers),label=typeline('lya'),lw=2.5)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_civ[i])['T_aniso_resp']],next(markers),label=typeline('civ'),lw=2.5)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_halpha[i])['T_aniso_resp']],next(markers),label=typeline('halpha'),lw=2.5)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_hbeta[i])['T_aniso_resp']],next(markers),label=typeline('hbeta'),lw=2.5)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_he4686[i])['T_aniso_resp']],next(markers),label=typeline('he4686'),lw=2.5)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_he1640[i])['T_aniso_resp']],next(markers),label=typeline('he1640'),lw=2.5)
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.ylabel(r'$\tau_\eta\mathrm{,\,anisotropic\,clouds\,(days)}$',fontsize=plot_fontsize)
plt.xlabel(r'$\log(n_\mathrm{H}/1 \mathrm{cm}^{-3})$',fontsize=plot_fontsize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.figtext(0.35,0.25,r'$\log(n_\mathrm{col})=\mathrm{'+str(plotNcol)+'}$',fontsize=plot_fontsize)
plt.axis([7,14,-100,350])
plt.axvline(x=10.75,color='black',lw=1.5)
plt.legend(loc='upper right',numpoints=1,fontsize=plot_fontsize)
plt.tight_layout()
plt.savefig('lineresults_s0/Tau_aniso_resp_nH.eps')
plt.close()

markersymbols = ['b','r--','k:','m-.','co','gs']
markers = cycle(markersymbols)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_lya[i])['R']],next(markers),label=typeline('lya'),lw=2.5)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_civ[i])['R']],next(markers),label=typeline('civ'),lw=2.5)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_halpha[i])['R']],next(markers),label=typeline('halpha'),lw=2.5)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_hbeta[i])['R']],next(markers),label=typeline('hbeta'),lw=2.5)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_he4686[i])['R']],next(markers),label=typeline('he4686'),lw=2.5)
plt.plot((data_lya[i])['nHarr'],[tau/lightday for tau in (data_he1640[i])['R']],next(markers),label=typeline('he1640'),lw=2.5)
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.xlabel(r'$\log(n_\mathrm{H}/1 \mathrm{cm}^{-3})$',fontsize=plot_fontsize)
plt.ylabel(r'$r_\epsilon\mathrm{\,(lightdays)}$',fontsize=plot_fontsize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.figtext(0.2,0.25,r'$\log(n_\mathrm{col})=\mathrm{'+str(plotNcol)+'}$',fontsize=plot_fontsize)
plt.axvline(x=10.75,color='black',lw=1.5)
plt.axis([7,14,0,180])
plt.tight_layout()
plt.legend(loc='upper right',numpoints=1,fontsize=plot_fontsize)
plt.savefig('lineresults_s0/r_isotropic_nH.eps')
plt.close()
