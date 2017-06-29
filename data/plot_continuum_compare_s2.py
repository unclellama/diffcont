from ewgrids import *
import math as m
import matplotlib.pyplot as plt
from matplotlib import ticker
from ewgrids import *
import json

contfiles = ['continuum_s2_startnH9_startNcol22.5.js','continuum_s2_startnH10_startNcol22.5.js',
             'continuum_s2_startnH11_startNcol22.5.js','continuum_s2_startnH11.5_startNcol22.5.js']
legend_list = [r'$\log(U)=0.52$',r'$\log(U)=-0.48$',r'$\log(U)=-1.48$',r'$\log(U)=-1.98$']
s = '2'

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

        diffcont_list.append(diffcont)
        contfrac_list.append(contfrac_tot)
        lag_list.append(lags)
        dil_lag_list.append(dil_lags)
        centroid_list.append(ccfcent)
        dil_centroid_list.append(dil_ccfcent)
        respfunc_cent_list.append(cent_tau_resp_aniso)

plot_fontsize = 20
plot_ticksize = 19
plot_lw = 2.5
xtick_padding = 15
tickwidth=2.5
markersize=14

from pylab import rcParams
rcParams['figure.figsize'] = 15, 11.4
rcParams['legend.handlelength'] = 0
max_tau = 100
nrows = 4

# plot of lumins, respcent

plt.plot([1,2,3])

plt.subplot(nrows,2,1)
lab_left = r'$L_\nu\,\mathrm{[erg\,s}^{-1}\mathrm{]}$'
plt.ylabel(lab_left,fontsize=plot_fontsize)
plt.plot(wl,ionizing_cont,'b:',label='Ionizing continuum',lw=plot_lw)
plt.plot(wl,diffcont_list[0],lw=plot_lw)
plt.axis([1000,10000,1e41,2*1e45])
plt.xticks([], [])
plt.yscale('log')
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.subplot(nrows,2,3)
plt.ylabel(lab_left,fontsize=plot_fontsize)
plt.plot(wl,ionizing_cont,'b:',label='Ionizing continuum',lw=plot_lw)
plt.plot(wl,diffcont_list[1],lw=plot_lw)
plt.axis([1000,10000,1e41,2*1e45])
plt.xticks([], [])
plt.yscale('log')
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.subplot(nrows,2,5)
plt.ylabel(lab_left,fontsize=plot_fontsize)
plt.plot(wl,ionizing_cont,'b:',label='Ionizing continuum',lw=plot_lw)
plt.plot(wl,diffcont_list[2],lw=plot_lw)
plt.axis([1000,10000,1e41,2*1e45])
plt.xticks([], [])
plt.yscale('log')
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.subplot(nrows,2,7)
plt.ylabel(lab_left,fontsize=plot_fontsize)
plt.plot(wl,ionizing_cont,'b:',label='Ionizing continuum',lw=plot_lw)
plt.plot(wl,diffcont_list[3],lw=plot_lw)
plt.axis([1000,10000,1e41,2*1e45])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.yscale('log')
plt.ylabel(lab_left,fontsize=plot_fontsize)
plt.xlabel(r'Wavelength [\AA]',fontsize=plot_fontsize)
plt.xticks(np.arange(2000, 12000, 2000))

plt.subplot(nrows,2,2)
lab_right = r'$\tau_\eta\,\mathrm{[days]}$'
plt.ylabel(lab_right,fontsize=plot_fontsize)
plt.plot(wl,respfunc_cent_list[0],lw=plot_lw,label=legend_list[0])
plt.axis([1000,10000,0,max_tau])
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.legend(markerscale=0)
plt.subplot(nrows,2,4)
plt.ylabel(lab_right,fontsize=plot_fontsize)
plt.plot(wl,respfunc_cent_list[1],lw=plot_lw,label=legend_list[1])
plt.axis([1000,10000,0,max_tau])
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.legend(markerscale=0)
plt.subplot(nrows,2,6)
plt.ylabel(lab_right,fontsize=plot_fontsize)
plt.plot(wl,respfunc_cent_list[2],lw=plot_lw,label=legend_list[2])
plt.axis([1000,10000,0,max_tau])
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.legend(markerscale=0)
plt.subplot(nrows,2,8)
plt.ylabel(lab_right,fontsize=plot_fontsize)
plt.plot(wl,respfunc_cent_list[3],lw=plot_lw,label=legend_list[3])
plt.axis([1000,10000,0,max_tau])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.legend(markerscale=0)
plt.xticks(np.arange(2000, 12000, 2000))
plt.xlabel(r'Wavelength [\AA]',fontsize=plot_fontsize)
plt.tight_layout()
plt.savefig('cont_comparison_lumin_resp_s2.eps')
plt.close()

# plot of the lags - centroid (left panels), and cent*F_diff (right)

plt.plot([1,2,3])
max_tau = 60
ylab = r'$\tau_{\mathrm{CCF,cent.}}\,\mathrm{[days]}$'
plt.subplot(nrows,2,1)
plt.ylabel(ylab,fontsize=plot_fontsize)
plt.plot(wl,centroid_list[0],lw=plot_lw)
plt.axis([1000,10000,0,max_tau])
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.subplot(nrows,2,3)
plt.ylabel(ylab,fontsize=plot_fontsize)
plt.plot(wl,centroid_list[1],lw=plot_lw)
plt.axis([1000,10000,0,max_tau])
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.subplot(nrows,2,5)
plt.ylabel(ylab,fontsize=plot_fontsize)
plt.plot(wl,centroid_list[2],lw=plot_lw)
plt.axis([1000,10000,0,max_tau])
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.subplot(nrows,2,7)
plt.ylabel(ylab,fontsize=plot_fontsize)
plt.plot(wl,centroid_list[3],lw=plot_lw)
plt.axis([1000,10000,0,max_tau])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.xlabel(r'Wavelength [\AA]',fontsize=plot_fontsize)
plt.ylabel(ylab,fontsize=plot_fontsize)
plt.xticks(np.arange(2000, 12000, 2000))

max_tau = 15
plt.subplot(nrows,2,2)
ylab_right = r'$\tau_{\mathrm{CCF,cent.}}\times F_{\mathrm{diff.}}$'
plt.ylabel(ylab_right,fontsize=plot_fontsize)
plt.plot(wl,dil_centroid_list[0],lw=plot_lw,label=legend_list[0])
plt.axis([1000,10000,0,max_tau])
plt.xticks([], [])
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.legend(markerscale=0)
plt.subplot(nrows,2,4)
plt.ylabel(ylab_right,fontsize=plot_fontsize)
plt.plot(wl,dil_centroid_list[1],lw=plot_lw,label=legend_list[1])
plt.axis([1000,10000,0,max_tau])
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
plt.legend(markerscale=0)
plt.subplot(nrows,2,8)
plt.ylabel(ylab_right,fontsize=plot_fontsize)
plt.plot(wl,dil_centroid_list[3],lw=plot_lw,label=legend_list[3])
plt.axis([1000,10000,0,max_tau])
plt.legend(markerscale=0)
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.xticks(np.arange(2000, 12000, 2000))
plt.legend(markerscale=0)
plt.xlabel(r'Wavelength [\AA]',fontsize=plot_fontsize)
plt.tight_layout()
plt.savefig('cont_comparison_lags_s2.eps')
plt.close()
