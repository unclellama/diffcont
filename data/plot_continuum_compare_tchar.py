from ewgrids import *
import math as m
import matplotlib.pyplot as plt
from matplotlib import ticker
from ewgrids import *
import json

contfiles = ['continuum_s0_nH10.75_Ncol22.5_Tchar5.js','continuum_s0_nH10.75_Ncol22.5_Tchar10.js',
             'continuum_s0_nH10.75_Ncol22.5_Tchar20.js',
             'continuum_s0_nH10.75_Ncol22.5.js','continuum_s0_nH10.75_Ncol22.5_Tchar80.js']
legend_list = [r'$T_\mathrm{char}=5$',r'$T_\mathrm{char}=10$',r'$T_\mathrm{char}=20$',r'$T_\mathrm{char}=40$',
               r'$T_\mathrm{char}=80$',r'$\log(n_\mathrm{H})=13$',r'$\log(n_\mathrm{H})=14$']
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
wl_list = []

respfunc_cent_list = []

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

        wl_list.append(wl)
        diffcont_list.append(diffcont)
        contfrac_list.append(contfrac_tot)
        lag_list.append(lags)
        dil_lag_list.append(dil_lags)
        centroid_list.append(ccfcent)
        dil_centroid_list.append(dil_ccfcent)
        respfunc_cent_list.append(cent_tau_resp_aniso)

#p_wl_multi(wl,contfrac_list,title='cont_compare_s0',log='no',
#         ylabel=r'$L_{\nu,\mathrm{diff}}/L_{\nu,\mathrm{tot}}$',label=r'',
#         legends=legend_list,axes=[1000,10000,0,1.1])
#p_wl_multi(wl,lag_list,title='cont_compare_s0_lags',log='no',
#         ylabel=r'CCF Lag (days)',label=r'',
#         legends=legend_list,axes=[1000,10000,0,140])
#p_wl_multi(wl,centroid_list,title='cont_compare_s0_centroids',log='no',
#         ylabel=r'CCF Centroid (days)',label=r'',
#         legends=legend_list,axes=[1000,10000,0,140])
#p_wl_multi(wl,dil_lag_list,title='cont_compare_s0_dil_lags',log='no',
#         ylabel=r'CCF Lag$\times F_{\mathrm{diff}}$ (days)',label=r'',
#         legends=legend_list,axes=[1000,10000,0,18])
#p_wl_multi(wl,dil_centroid_list,title='cont_compare_s0_dil_centroids',log='no',
#         ylabel=r'CCF Centroid $\times F_{\mathrm{diff}}$ (days)',label=r'',
#         legends=legend_list,axes=[1000,10000,0,18])

plot_fontsize = 20
plot_ticksize = 19
plot_lw = 2.5
xtick_padding = 15
tickwidth=2.5
markersize=14

#plt.plot(wl,ionizing_cont,'b:',label='Ionizing continuum',lw=plot_lw)
#plt.plot(wl,diffcont_list[0],lw=plot_lw,label=r'$\log(n_{\mathrm{H}}=8)$')
#plt.plot(wl,diffcont_list[1],lw=plot_lw,label=r'$\log(n_{\mathrm{H}}=9)$')
#plt.plot(wl,diffcont_list[2],lw=plot_lw,label=r'$\log(n_{\mathrm{H}}=10)$')
#plt.plot(wl,diffcont_list[3],lw=plot_lw,label=r'$\log(n_{\mathrm{H}}=11)$')
#plt.plot(wl,diffcont_list[4],lw=plot_lw,label=r'$\log(n_{\mathrm{H}}=12)$')
#plt.plot(wl,diffcont_list[5],lw=plot_lw,label=r'$\log(n_{\mathrm{H}}=13)$')
#plt.xlabel(r'Wavelength [\AA]',fontsize=plot_fontsize)
#plt.ylabel(r'$L_\nu$ (erg s$^{-1}$)',fontsize=plot_fontsize)
#plt.xticks(np.arange(2000, 12000, 2000))
#plt.axis([1000,10000,2*1e41,2*1e45])
#plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
#plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
#plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
#plt.yscale('log')
#plt.legend(ncol=2)
#plt.show()
#plt.close()

from pylab import rcParams
rcParams['figure.figsize'] = 15, 20
rcParams['legend.handlelength'] = 0
max_tau = 100
nrows = 5

# plot of lumins, respcent

plt.plot([1,2,3])

plt.subplot(nrows,2,1)
lab_left = r'$L_\nu\,\mathrm{[erg\,s}^{-1}\mathrm{]}$'
plt.ylabel(lab_left,fontsize=plot_fontsize)
plt.plot(wl,ionizing_cont,'b:',label='Ionizing continuum',lw=plot_lw)
plt.plot(wl,diffcont_list[0],lw=plot_lw,label=r'$\log(n_{\mathrm{H}}=8)$')
plt.axis([1000,10000,1e41,2*1e45])
plt.xticks([], [])
plt.yscale('log')
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.subplot(nrows,2,3)
plt.ylabel(lab_left,fontsize=plot_fontsize)
plt.plot(wl,ionizing_cont,'b:',label='Ionizing continuum',lw=plot_lw)
plt.plot(wl,diffcont_list[1],lw=plot_lw,label=r'$\log(n_{\mathrm{H}}=9)$')
plt.axis([1000,10000,1e41,2*1e45])
plt.xticks([], [])
plt.yscale('log')
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.subplot(nrows,2,5)
plt.ylabel(lab_left,fontsize=plot_fontsize)
plt.plot(wl,ionizing_cont,'b:',label='Ionizing continuum',lw=plot_lw)
plt.plot(wl,diffcont_list[2],lw=plot_lw,label=r'$\log(n_{\mathrm{H}}=10)$')
plt.axis([1000,10000,1e41,2*1e45])
plt.xticks([], [])
plt.yscale('log')
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.subplot(nrows,2,7)
plt.ylabel(lab_left,fontsize=plot_fontsize)
plt.plot(wl,ionizing_cont,'b:',label='Ionizing continuum',lw=plot_lw)
print(len(diffcont_list[3]))
plt.plot(wl_list[3],diffcont_list[3],lw=plot_lw,label=r'$\log(n_{\mathrm{H}}=11)$')
plt.axis([1000,10000,1e41,2*1e45])
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.yscale('log')
plt.ylabel(lab_left,fontsize=plot_fontsize)
plt.subplot(nrows,2,9)
plt.plot(wl,ionizing_cont,'b:',label='Ionizing continuum',lw=plot_lw)
plt.plot(wl,diffcont_list[4],lw=plot_lw,label=r'$\log(n_{\mathrm{H}}=14)$')
plt.xlabel(r'Wavelength [\AA]',fontsize=plot_fontsize)
plt.ylabel(lab_left,fontsize=plot_fontsize)
plt.xticks(np.arange(2000, 12000, 2000))
plt.axis([1000,10000,1e41,2*1e45])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.yscale('log')
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)

plt.subplot(nrows,2,2)
lab_right = r'$\tau_\eta\,\mathrm{[days]}$'
plt.ylabel(lab_right,fontsize=plot_fontsize)
plt.plot(wl,respfunc_cent_list[0],lw=plot_lw,label=legend_list[0])
plt.axis([1000,10000,0,max_tau])
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.legend(markerscale=0)
plt.subplot(nrows,2,4)
plt.ylabel(lab_right,fontsize=plot_fontsize)
plt.plot(wl,respfunc_cent_list[1],lw=plot_lw,label=legend_list[1])
plt.axis([1000,10000,0,max_tau])
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.legend(markerscale=0)
plt.subplot(nrows,2,6)
plt.ylabel(lab_right,fontsize=plot_fontsize)
plt.plot(wl,respfunc_cent_list[2],lw=plot_lw,label=legend_list[2])
plt.axis([1000,10000,0,max_tau])
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.legend(markerscale=0)
plt.subplot(nrows,2,8)
plt.ylabel(lab_right,fontsize=plot_fontsize)
plt.plot(wl_list[3],respfunc_cent_list[3],lw=plot_lw,label=legend_list[3])
plt.axis([1000,10000,0,max_tau])
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.legend(markerscale=0)
plt.subplot(nrows,2,10)
plt.ylabel(lab_right,fontsize=plot_fontsize)
plt.plot(wl,respfunc_cent_list[4],lw=plot_lw,label=legend_list[4])
plt.axis([1000,10000,0,max_tau])
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.xticks(np.arange(2000, 12000, 2000))
plt.legend(markerscale=0)
plt.xlabel(r'Wavelength [\AA]',fontsize=plot_fontsize)
plt.tight_layout()
plt.savefig('cont_comparison_lumin_resp.eps')
plt.close()

# plot of the lags - centroid (left panels), and cent*F_diff (right)

plt.plot([1,2,3])
max_tau = 60
ylab = r'$\tau_{\mathrm{CCF,cent.}}\,\mathrm{[days]}$'
plt.subplot(nrows,2,1)
plt.ylabel(ylab,fontsize=plot_fontsize)
plt.plot(wl,centroid_list[0],lw=plot_lw,label=r'$\log(n_{\mathrm{H}}=8)$')
plt.axis([1000,10000,0,max_tau])
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.subplot(nrows,2,3)
plt.ylabel(ylab,fontsize=plot_fontsize)
plt.plot(wl,centroid_list[1],lw=plot_lw,label=r'$\log(n_{\mathrm{H}}=9)$')
plt.axis([1000,10000,0,max_tau])
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.subplot(nrows,2,5)
plt.ylabel(ylab,fontsize=plot_fontsize)
plt.plot(wl,centroid_list[2],lw=plot_lw,label=r'$\log(n_{\mathrm{H}}=10)$')
plt.axis([1000,10000,0,max_tau])
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.subplot(nrows,2,7)
plt.ylabel(ylab,fontsize=plot_fontsize)
plt.plot(wl_list[3],centroid_list[3],lw=plot_lw,label=r'$\log(n_{\mathrm{H}}=11)$')
plt.axis([1000,10000,0,max_tau])
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.ylabel(ylab,fontsize=plot_fontsize)
plt.subplot(nrows,2,9)
plt.plot(wl,centroid_list[4],lw=plot_lw,label=r'$\log(n_{\mathrm{H}}=14)$')
plt.xlabel(r'Wavelength [\AA]',fontsize=plot_fontsize)
plt.ylabel(ylab,fontsize=plot_fontsize)
plt.xticks(np.arange(2000, 12000, 2000))
plt.axis([1000,10000,0,max_tau])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)

max_tau = 15
plt.subplot(nrows,2,2)
ylab_right = r'$\tau_{\mathrm{CCF,cent.}}\times F_{\mathrm{diff.}}$'
plt.ylabel(ylab_right,fontsize=plot_fontsize)
plt.plot(wl,dil_centroid_list[0],lw=plot_lw,label=legend_list[0])
plt.axis([1000,10000,0,max_tau])
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.legend(markerscale=0)
plt.subplot(nrows,2,4)
plt.ylabel(ylab_right,fontsize=plot_fontsize)
plt.plot(wl,dil_centroid_list[1],lw=plot_lw,label=legend_list[1])
plt.axis([1000,10000,0,max_tau])
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.legend(markerscale=0)
plt.subplot(nrows,2,6)
plt.ylabel(ylab_right,fontsize=plot_fontsize)
plt.plot(wl,dil_centroid_list[2],lw=plot_lw,label=legend_list[2])
plt.axis([1000,10000,0,max_tau])
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.legend(markerscale=0)
plt.subplot(nrows,2,8)
plt.ylabel(ylab_right,fontsize=plot_fontsize)
plt.plot(wl_list[3],dil_centroid_list[3],lw=plot_lw,label=legend_list[3])
plt.axis([1000,10000,0,max_tau])
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.legend(markerscale=0)
plt.subplot(nrows,2,10)
plt.ylabel(ylab_right,fontsize=plot_fontsize)
plt.plot(wl,dil_centroid_list[4],lw=plot_lw,label=legend_list[4])
plt.axis([1000,10000,0,max_tau])
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.xticks(np.arange(2000, 12000, 2000))
plt.legend(markerscale=0)
plt.xlabel(r'Wavelength [\AA]',fontsize=plot_fontsize)
plt.tight_layout()
plt.savefig('cont_comparison_lags.eps')
plt.close()
