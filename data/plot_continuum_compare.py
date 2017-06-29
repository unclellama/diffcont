from ewgrids import *
import math as m
import matplotlib.pyplot as plt
from matplotlib import ticker
from ewgrids import *
import json

contfiles = ['continuum_s0_nH9_Ncol22.5.js','continuum_s0_nH10_Ncol22.5.js','continuum_s0_nH11_Ncol22.5.js',
             'continuum_s0_nH12_Ncol22.5.js','continuum_s0_nH13_Ncol22.5.js']
legend_list = [r'$\log(n_\mathrm{H})=9$',r'$\log(n_\mathrm{H})=10$',r'$\log(n_\mathrm{H})=11$',
               r'$\log(n_\mathrm{H})=12$',r'$\log(n_\mathrm{H})=13$']
s = '0'

plot_fontsize = 20
plot_ticksize = 19
plot_lwidth = 2.5
xtick_padding = 15
tickwidth=2.5
markersize=14

diffcont_list = []
contfrac_list = []
lag_list = []
dil_lag_list = []
centroid_list = []
dil_centroid_list = []


for contfile in contfiles:
    with open(contfile) as json_data:
        print('contfile:',contfile)
        contdata = json.load(json_data)
        wl = contdata['contwl']
        diffcont = contdata['diffcont']
        cent = contdata['cent']
        respcent = contdata['respcent']
        ionizing_cont = contdata['ionizing_cont']
        cont_ratio = contdata['cont_ratio']
        weighted_resp = contdata['weighted_resp']
        weighted_inwd = contdata['weighted_inwd']
        cent_tau_aniso = contdata['cent_tau_aniso']
        cent_tau_resp_aniso = contdata['cent_tau_resp_aniso']
        lags = contdata['mean_lags']
        e_lags = contdata['std_lags']
        ccfcent = contdata['mean_centroids']
        e_ccfcent = contdata['std_centroids']
        ii = range(len(diffcont))
        total_cont = [diffcont[i]+ionizing_cont[i] for i in ii]
        contfrac_tot = [diffcont[i]/total_cont[i] for i in ii]
        dil_cent = [contfrac_tot[i]*cent[i] for i in ii]
        dil_lags = [contfrac_tot[i]*lags[i] for i in ii]
        dil_ccfcent = [contfrac_tot[i]*ccfcent[i] for i in ii]
        dil_respcent = [contfrac_tot[i]*respcent[i] for i in ii]
        dil_cent_tau_aniso = [contfrac_tot[i]*cent_tau_aniso[i] for i in ii]
        dil_cent_tau_resp_aniso = [contfrac_tot[i]*cent_tau_resp_aniso[i] for i in ii]
        for k in range(len(wl)):
            print('wl:',wl[k],' lag:',lags[k])

        diffcont_list.append(diffcont)
        contfrac_list.append(contfrac_tot)
        lag_list.append(lags)
        dil_lag_list.append(dil_lags)
        centroid_list.append(ccfcent)
        dil_centroid_list.append(dil_ccfcent)

p_wl_multi(wl,contfrac_list,title='cont_compare_s0',log='no',
         ylabel=r'$L_{\nu,\mathrm{diff}}/L_{\nu,\mathrm{tot}}$',label=r'',
         legends=legend_list,axes=[1000,10000,0,1.1])
p_wl_multi(wl,lag_list,title='cont_compare_s0_lags',log='no',
         ylabel=r'CCF Lag (days)',label=r'',
         legends=legend_list,axes=[1000,10000,0,140])
p_wl_multi(wl,centroid_list,title='cont_compare_s0_centroids',log='no',
         ylabel=r'CCF Centroid (days)',label=r'',
         legends=legend_list,axes=[1000,10000,0,140])
p_wl_multi(wl,dil_lag_list,title='cont_compare_s0_dil_lags',log='no',
         ylabel=r'CCF Lag$\times F_{\mathrm{diff}}$ (days)',label=r'',
         legends=legend_list,axes=[1000,10000,0,18])
p_wl_multi(wl,dil_centroid_list,title='cont_compare_s0_dil_centroids',log='no',
         ylabel=r'CCF Centroid $\times F_{\mathrm{diff}}$ (days)',label=r'',
         legends=legend_list,axes=[1000,10000,0,18])

