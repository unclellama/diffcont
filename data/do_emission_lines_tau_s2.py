# =======================================================================================
# get BEL luminosities + weighted radii for a grid of nH, Ncol, for s=0 models
# - basically calls do_line() for each line, nH, Ncol combination.
# =======================================================================================

r_in = 1. # in lightdays

from ewgrids import *

L_lya,eL_lya = observed_lumin('lya')
L_civ,eL_civ = observed_lumin('civ')
L_hbeta,eL_hbeta = observed_lumin('hbeta')
L_he4686,eL_he4686 = observed_lumin('he4686')
start_nHlist = [9.5,9.75,10,10.25,10.5,10.75,11,11.25]
start_logPhi = 19.5
#start_nHlist = [9.5,10]

logU = [start_logPhi-nH-m.log10(c_cms) for nH in start_nHlist]

L_lya_U = []
L_civ_U = []
L_hbeta_U = []
L_he4686_U = []
RL_lya_U = []
RL_civ_U = []
RL_hbeta_U = []
RL_he4686_U = []
RR_lya_U = []
RR_civ_U = []
RR_hbeta_U = []
RR_he4686_U = []
RL_an_lya_U = []
RL_an_civ_U = []
RL_an_hbeta_U = []
RL_an_he4686_U = []
RR_an_lya_U = []
RR_an_civ_U = []
RR_an_hbeta_U = []
RR_an_he4686_U = []
for start_nH in start_nHlist:
        L,cent,respcent,cent_tau_aniso,cent_tau_resp_aniso = do_line_s2(
            starting_lognH=start_nH,starting_logPhi=start_logPhi,line='lya')
        L_lya_U.append(L)
        RL_lya_U.append(cent/lightday)
        RR_lya_U.append(respcent/lightday)
        RL_an_lya_U.append(cent_tau_aniso/lightday)
        RR_an_lya_U.append(cent_tau_resp_aniso/lightday)
        L,cent,respcent,cent_tau_aniso,cent_tau_resp_aniso = do_line_s2(
            starting_lognH=start_nH,starting_logPhi=start_logPhi,line='civ')
        L_civ_U.append(L)
        RL_civ_U.append(cent/lightday)
        RR_civ_U.append(respcent/lightday)
        RL_an_civ_U.append(cent_tau_aniso/lightday)
        RR_an_civ_U.append(cent_tau_resp_aniso/lightday)
        L,cent,respcent,cent_tau_aniso,cent_tau_resp_aniso = do_line_s2(
            starting_lognH=start_nH,starting_logPhi=start_logPhi,line='hbeta')
        L_hbeta_U.append(L)
        RL_hbeta_U.append(cent/lightday)
        RR_hbeta_U.append(respcent/lightday)
        RL_an_hbeta_U.append(cent_tau_aniso/lightday)
        RR_an_hbeta_U.append(cent_tau_resp_aniso/lightday)
        L,cent,respcent,cent_tau_aniso,cent_tau_resp_aniso = do_line_s2(
            starting_lognH=start_nH,starting_logPhi=start_logPhi,line='he4686')
        L_he4686_U.append(L)
        RL_he4686_U.append(cent/lightday)
        RR_he4686_U.append(respcent/lightday)
        RL_an_he4686_U.append(cent_tau_aniso/lightday)
        RR_an_he4686_U.append(cent_tau_resp_aniso/lightday)

from matplotlib import pyplot as plt
plot_fontsize = 20
legend_fontsize = 20
plot_ticksize = 21
plot_lwidth = 2.5
xtick_padding = 15
tickwidth=2.5
markersize=14
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.plot(logU,L_lya_U,label=typeline('lya'))
plt.plot(logU,L_civ_U,label=typeline('civ'))
plt.plot(logU,L_hbeta_U,label=typeline('hbeta'))
plt.plot(logU,L_he4686_U,label=typeline('he4686'))
plt.axhline(y=L_lya,color='blue',ls='--')
plt.axhline(y=L_civ,color='green',ls='--')
plt.axhline(y=L_hbeta,color='red',ls='--')
plt.axhline(y=L_he4686,color='cyan',ls='--')
plt.yscale('log')
plt.xlabel(r'U',fontsize=plot_fontsize)
plt.ylabel(r'$L_{\mathrm{line}}$ [erg s$^{-1}$]',fontsize=plot_fontsize)
plt.legend(loc='upper right',fontsize=plot_fontsize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.tight_layout()
plt.savefig('Lline_U.eps')
plt.close()

plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.plot(logU,RL_lya_U,label=typeline('lya'))
plt.plot(logU,RL_civ_U,label=typeline('civ'))
plt.plot(logU,RL_hbeta_U,label=typeline('hbeta'))
plt.plot(logU,RL_he4686_U,label=typeline('he4686'))
plt.xlabel(r'U',fontsize=plot_fontsize)
plt.ylabel(r'Em-weighted Centroid, lightdays',fontsize=plot_fontsize)
plt.legend(loc='upper right',fontsize=plot_fontsize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.tight_layout()
plt.savefig('RLline_U.eps')
plt.close()

plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.plot(logU,RR_lya_U,label=typeline('lya'))
plt.plot(logU,RR_civ_U,label=typeline('civ'))
plt.plot(logU,RR_hbeta_U,label=typeline('hbeta'))
plt.plot(logU,RR_he4686_U,label=typeline('he4686'))
plt.xlabel(r'U',fontsize=plot_fontsize)
plt.ylabel(r'Resp-weighted Centroid, lightdays',fontsize=plot_fontsize)
plt.legend(loc='upper right',fontsize=plot_fontsize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.tight_layout()
plt.savefig('RRline_U.eps')
plt.close()

plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.plot(logU,RL_an_lya_U,label=typeline('lya'))
plt.plot(logU,RL_an_civ_U,label=typeline('civ'))
plt.plot(logU,RL_an_hbeta_U,label=typeline('hbeta'))
plt.plot(logU,RL_an_he4686_U,label=typeline('he4686'))
plt.xlabel(r'U',fontsize=plot_fontsize)
plt.ylabel(r'Aniso. Em-weighted Centroid, lightdays',fontsize=plot_fontsize)
plt.legend(loc='upper right',fontsize=plot_fontsize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.tight_layout()
plt.savefig('tauLline_U.eps')
plt.close()

plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.plot(logU,RR_an_lya_U,label=typeline('lya'))
plt.plot(logU,RR_an_civ_U,label=typeline('civ'))
plt.plot(logU,RR_an_hbeta_U,label=typeline('hbeta'))
plt.plot(logU,RR_an_he4686_U,label=typeline('he4686'))
plt.xlabel(r'U',fontsize=plot_fontsize)
plt.ylabel(r'Aniso. Resp-weighted Centroid, lightdays',fontsize=plot_fontsize)
plt.legend(loc='upper right',fontsize=plot_fontsize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.tight_layout()
plt.savefig('tauRline_U.eps')
plt.close()
