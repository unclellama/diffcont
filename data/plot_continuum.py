
contfile = 'continuum_s0_nH10.75_Ncol22.5_Tchar10'+'.js'
# contfile = 'continuum_s0_nH10.75_Ncol22.5_Tchar5'+'.js'
s = r'$0$'
#plotlabel = r'Model 1 $(s=$'+s+',$\log(n_{\mathrm{H}})=10.75$)'
plotlabel = r'Model 1 $(s=$'+s+',$\log(n_{\mathrm{H}})=10.75$,Tchar=5)'
#plotlabel = r'Model 2 $(s=$'+s+',$\log(U)=-1.23$)'

from ewgrids import *
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

with open(contfile) as json_data:
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

p_wl_multi(wl,[diffcont,ionizing_cont,total_cont],title='diff_cont_lumin_s'+s,log='yes',
         ylabel=r'Continuum $\nu L_\nu$ [erg s$^{-1}$]',label=plotlabel,
         legends=[r'Diffuse continuum','Ionizing continuum','Total continuum'],axes=[1000,10000,1e42,1e44])
p_wl(wl,cont_ratio,title='cont_frac_s'+s,label=plotlabel,axes=[1000,10000,0,1.1],
         ylabel=r'$\nu L_{\nu,\mathrm{diff.}}/\nu L_{\nu,\mathrm{nuc.}}$',log='no')
p_wl(wl,contfrac_tot,title='cont_frac_tot_s'+s,label=plotlabel,axes=[1000,10000,0,1],
         ylabel=r'Diffuse continuum fraction, $F_{\mathrm{diff}}$',log='no')
p_wl(wl,weighted_resp,title='weighted_resp_s'+s,label=plotlabel,
         ylabel=r'Flux-weighted responsivity, $\eta_\epsilon$',log='no')
p_wl(wl,weighted_inwd,title='inwd_frac_s'+s,label=plotlabel,
         ylabel=r'$\epsilon_{\mathrm{in}}/\epsilon_{\mathrm{tot}}(\lambda)$',log='no')
p_wl_multi(wl,[cent,respcent,cent_tau_aniso,cent_tau_resp_aniso],
               title='diff_cont_r_s'+s,log='no',label=plotlabel,
               ylabel=r'Centroid [lightdays]',axes=[1000,10000,0,80],
               legends=[r'$r_{\epsilon}$',r'$r_{\eta}$',r'$\tau_{\epsilon}$ (anisotropic)',
                        r'$\tau_{\eta}$  (anisotropic)'])
p_wl_multi(wl,[lags,ccfcent],title='diff_cont_lags_s'+s,log='no',label=plotlabel,
               ylabel=r'Lag (days)',axes=[1000,10000,0,80],
               legends=['CCF Peak','CCF Centroid'])
p_wl_multi(wl,[dil_lags,dil_ccfcent],
               title='diluted_cont_lags_s'+s,log='no',label=plotlabel,
               ylabel=r'Diluted lag (days)',axes=[1000,10000,0,15],
               legends=[r'CCF Peak $\times\,F_{\mathrm{diff}}$',r'CCF Centroid $\times\,F_{\mathrm{diff}}$'])
#p_wl_multi(wl,[dil_cent,dil_respcent,dil_cent_tau_aniso,dil_cent_tau_resp_aniso],
#                title='diluted_cont_r_s'+s,log='no',label=plotlabel,
#                ylabel=r'BLR Centroid [lightdays]',
#                legends=[r'$r_{\epsilon}$',r'$r_{\eta}$',r'$r_{\epsilon}$  (anisotropic)',
#                        r'$r_{\eta}$  (anisotropic)'])
