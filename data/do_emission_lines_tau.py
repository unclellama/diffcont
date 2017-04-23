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
Ncollist = ['22','23','24']
#nHlist = [8,9]
#nHplot = [8,9]
nHlist = [7,7.25,7.5,7.75,8,8.25,8.5,8.75,
          9,9.25,9.5,9.75,10,10.25,10.5,10.75,11,11.25,11.5,11.75,
          12,12.25,12.5,12.75,13,13.25,13.5,13.75,14]
nHplot = [7,'','','',8,'','','','','','','','','','','',
          11,'','','','','','','',13,'','','',14]
if len(nHlist) != len(nHplot):
    raise Exception('you screwed up!')
L_lya_Ncol = []
L_civ_Ncol = []
L_hbeta_Ncol = []
L_he4686_Ncol = []
RL_lya_Ncol = []
RL_civ_Ncol = []
RL_hbeta_Ncol = []
RL_he4686_Ncol = []
RR_lya_Ncol = []
RR_civ_Ncol = []
RR_hbeta_Ncol = []
RR_he4686_Ncol = []
for Ncol in Ncollist:
    L_lya_nH = []
    L_civ_nH = []
    L_hbeta_nH = []
    L_he4686_nH = []
    RL_lya_nH = []
    RL_civ_nH = []
    RL_hbeta_nH = []
    RL_he4686_nH = []
    RR_lya_nH = []
    RR_civ_nH = []
    RR_hbeta_nH = []
    RR_he4686_nH = []
    for nH in nHlist:
        print('\nNcol:',Ncol,'nH:',nH)
        L,cent,respcent,cent_tau_aniso,cent_tau_resp_aniso = do_line(
            Ncol=Ncol,line='lya',nH=nH,r_in=r_in)
        L_lya_nH.append(L)
        RL_lya_nH.append(cent_tau_aniso/lightday)
        RR_lya_nH.append(cent_tau_resp_aniso/lightday)
        L,cent,respcent,cent_tau_aniso,cent_tau_resp_aniso = do_line(
            Ncol=Ncol,line='civ',nH=nH,r_in=r_in)
        L_civ_nH.append(L)
        RL_civ_nH.append(cent_tau_aniso/lightday)
        RR_civ_nH.append(cent_tau_resp_aniso/lightday)
        L,cent,respcent,cent_tau_aniso,cent_tau_resp_aniso = do_line(
            Ncol=Ncol,line='hbeta',nH=nH,r_in=r_in)
        L_hbeta_nH.append(L)
        RL_hbeta_nH.append(cent_tau_aniso/lightday)
        RR_hbeta_nH.append(cent_tau_resp_aniso/lightday)
        L,cent,respcent,cent_tau_aniso,cent_tau_resp_aniso = do_line(
            Ncol=Ncol,line='he4686',nH=nH,r_in=r_in)
        L_he4686_nH.append(L)
        RL_he4686_nH.append(cent_tau_aniso/lightday)
        RR_he4686_nH.append(cent_tau_resp_aniso/lightday)
    L_lya_Ncol.append(L_lya_nH)
    L_civ_Ncol.append(L_civ_nH)
    L_hbeta_Ncol.append(L_hbeta_nH)
    L_he4686_Ncol.append(L_he4686_nH)
    RL_lya_Ncol.append(RL_lya_nH)
    RL_civ_Ncol.append(RL_civ_nH)
    RL_hbeta_Ncol.append(RL_hbeta_nH)
    RL_he4686_Ncol.append(RL_he4686_nH)
    RR_lya_Ncol.append(RR_lya_nH)
    RR_civ_Ncol.append(RR_civ_nH)
    RR_hbeta_Ncol.append(RR_hbeta_nH)
    RR_he4686_Ncol.append(RR_he4686_nH)
print(nHlist,RL_civ_Ncol[0])
plot_lya_civ_allnH(Ncollist,nHplot,L_lya_Ncol,L_civ_Ncol,L_lya,eL_lya,L_civ,
                   eL_civ,r_in=1)
plot_lya_hbeta_allnH(Ncollist,nHplot,L_lya_Ncol,L_hbeta_Ncol,L_lya,eL_lya,L_hbeta,
                   eL_hbeta,r_in=1)
plot_lya_he4686_allnH(Ncollist,nHplot,L_lya_Ncol,L_he4686_Ncol,L_lya,eL_lya,L_he4686,
                   eL_he4686,r_in=1)

from matplotlib import pyplot as plt
plot_fontsize = 25
legend_fontsize = 20
plot_ticksize = 24
plot_lwidth = 2.5
xtick_padding = 15
tickwidth=2.5
markersize=14
plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.plot(nHlist,L_lya_Ncol[1],label=r'Ly$\alpha$')
plt.plot(nHlist,L_civ_Ncol[1],label=r'C$IV$')
plt.plot(nHlist,L_hbeta_Ncol[1],label=r'H$\beta$')
plt.yscale('log')
plt.xlabel(r'$\log(n_H/1$cm$^{-3})$',fontsize=plot_fontsize)
plt.ylabel(r'$L_{\mathrm{line}}$ [erg s$^{-1}$]',fontsize=plot_fontsize)
plt.axis([6.5,14.5,5e37,5e43])
plt.legend(loc='lower right',fontsize=plot_fontsize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.figtext(0.3,0.3,'Anisotropic clouds',backgroundcolor='white',
            fontsize=plot_fontsize-5)
plt.tight_layout()
plt.savefig('L_nH_aniso.eps')
plt.close()

plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.plot(nHlist,RL_lya_Ncol[1],label=r'Ly$\alpha$')
plt.plot(nHlist,RL_civ_Ncol[1],label=r'C$IV$')
plt.plot(nHlist,RL_hbeta_Ncol[1],label=r'H$\beta$')
plt.xlabel(r'$\log(n_H/1$cm$^{-3})$',fontsize=plot_fontsize)
plt.ylabel(r'$\epsilon$-weighted delay [lightdays]',fontsize=plot_fontsize)
plt.axis([6.5,14.5,0,140])
plt.legend(loc='upper right',fontsize=plot_fontsize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.figtext(0.3,0.3,'Anisotropic clouds',backgroundcolor='white',
            fontsize=plot_fontsize-5)
plt.tight_layout()
plt.savefig('RL_nH_tau.eps')
plt.close()

plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
plt.plot(nHlist,RR_lya_Ncol[1],label=r'Ly$\alpha$')
plt.plot(nHlist,RR_civ_Ncol[1],label=r'C$IV$')
plt.plot(nHlist,RR_hbeta_Ncol[1],label=r'H$\beta$')
plt.xlabel(r'$\log(\mathrm{nH}/1$cm$^{-3})$',fontsize=plot_fontsize)
plt.ylabel(r'$\eta$-weighted delay [lightdays]',fontsize=plot_fontsize)
plt.axis([6.5,14.5,0,300])
plt.legend(loc='upper right',fontsize=plot_fontsize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
plt.figtext(0.3,0.3,'Anisotropic clouds',backgroundcolor='white',
            fontsize=plot_fontsize-5)
plt.tight_layout()
plt.savefig('RR_nH_tau.eps')
plt.close()
