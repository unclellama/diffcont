# codes to work with the NGC 5548 photoionization grids supplied by Otho Adam Ulrich
import math as m
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
from itertools import cycle

plot_fontsize = 20
plot_ticksize = 19
plot_lwidth = 2.5
xtick_padding = 15
tickwidth=2.5
markersize=14

#import seaborn as sns

lightday = 2.59e15 # cm - lazy/evil duplication...

def p_grid(lognH,logPhi,lognH_cut,logPhi_cut,starting_lognH,starting_logPhi):
    """ plot an EW grid with a point on it and a line drawn through.
    not sure if this is used anywhere! probably kill me."""
    plt.plot(lognH,logPhi,color='white')
    plt.plot(lognH_cut,logPhi_cut,color='red')
    plt.plot([starting_lognH]*2,[starting_logPhi]*2,'*b')
    plt.xlabel(r'$\log[n_H/1$cm$^{-3}]$')
    plt.ylabel(r'$\log[Phi/1$ cm$^{-2}$ s$^{-1}]$')
    plt.show()
    plt.close()

def p_EWcontours_overlay(range_Hden,range_Phi,energy_grid,lognH_cut,logPhi_cut,
                         Ncol_cut_descending):
    """ this is for the s=2 models: plots the column density along the slice,
    overlaid on EW contours for whatever the starting-point Ncol is. """
    plt.figure()
    CS = plt.contourf(range_Hden, range_Phi, energy_grid,
                      locator=ticker.LogLocator())
    plt.colorbar(CS)
    plt.plot(lognH_cut,logPhi_cut,color='black')
    #plt.xlabel(r'$\log[$n_H / 1 cm$^{-3}]$')
    #plt.ylabel(r'$\log[\Phi$(H) / 1 cm$^{-2}$ s$^{-1}]$',rotation=90)
    Ncol_cut_descending = ["{0:.2f}".format(m.log10(N)) if (N > 21.00 & N<24.25) else '' for
                           N in Ncol_cut_descending]
    for i in range(len(lognH_cut)):
        plt.text(lognH_cut[i],logPhi_cut[i],Ncol_cut_descending[i])
    plt.savefig('EWcontours_overlay.eps', format='eps')
    plt.close()

def p_EWcontours_radius_overlay(range_Hden,range_radius,energy_grid,lognH_cut,logr_cut,
                         Ncol_cut_descending,r_in=1*lightday,r_out=140*lightday):
    """ same as p_EWcontours_overlay, but with log radius up the the y axis."""
    plt.figure()
    CS = plt.contourf(range_Hden,range_radius,energy_grid,
                      locator=ticker.LogLocator())
    plt.plot(lognH_cut,logr_cut,color='black')
    plt.axis([7,14,m.log10(300*lightday),m.log10(0.3*lightday)])
    plt.colorbar(CS)
    plt.xlabel('$\log[$n_H / 1 cm$^{-3}]$')
    plt.ylabel(r'$\log[R / 1 cm]$',rotation=90)
    Ncol_cut_descending = ["{0:.2f}".format(m.log10(N)) if ((m.log10(N) > 21.00) & (m.log10(N) < 24.25))
                            else '' for N in Ncol_cut_descending]
    for i in range(len(lognH_cut)):
        plt.text(lognH_cut[i],logr_cut[i],Ncol_cut_descending[i])
    plt.axhline(y=m.log10(r_in),color='black')
    plt.axhline(y=m.log10(14.8*2.59e15),color='pink')
    plt.axhline(y=m.log10(r_out),color='black')
    plt.savefig('EWcontours_overlay_radius.eps', format='eps')
    plt.close()
    
def p_ld(r_ld,y,title='something vs r(lightdays)',vlines='none',
         log='yes',axes='none',hlines='none',label='a label',legends=['1','2','3','4'],
         xtitle='$R$ (lightdays)',ytitle='y',xlog='no',label_xy=[0.3,0.3]):
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
        plt.legend(fontsize=plot_fontsize-5,loc='best')
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

def p_wl(wl,y,ylabel='y',title='something vs wavelength',xmax=7000,log='yes',xlog='no',
         legends='none',label='',range='none',axes='none'):
    plt.plot(wl,y,lw=2.5)
    if log == 'yes':
        plt.yscale('log')
    if xlog == 'yes':
        plt.xscale('log')
    #plt.title(title)
    plt.xlabel(r'Wavelength [\AA]',fontsize=plot_fontsize)
    plt.ylabel(ylabel,rotation=90,fontsize=plot_fontsize)
    if xlog == 'no':
        plt.xticks(np.arange(2000, 12000, 2000))
    plt.figtext(0.5,0.2,label,backgroundcolor='white',fontsize=plot_fontsize-2)
    if range != 'none':
        plt.axis(range)
    plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
    plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
    plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
    if axes != 'none':
        plt.axis(axes)
    plt.tight_layout()
    plt.savefig(title+'.eps')
    plt.close()

def p_wl_multi(wl,yarrays,ylabel='y',title='something vs wavelength',xmax=7000,log='yes',
         legends='none',axes='none',label='label'):
    markersymbols = ['b','r--','k:','m-.','g','cv']
    markers = cycle(markersymbols)
    for i in range(len(yarrays)):
        plt.plot(wl,yarrays[i],next(markers),label=legends[i],lw=2.5)
    if log == 'yes':
        plt.yscale('log')
    plt.xlabel(r'Wavelength [\AA]',fontsize=plot_fontsize)
    plt.ylabel(ylabel,rotation=90,fontsize=plot_fontsize)
    plt.figtext(0.5,0.2,label,backgroundcolor='white',fontsize=plot_fontsize-2)
    plt.legend(fontsize=plot_fontsize-5)
    plt.xticks(np.arange(2000, 12000, 2000))
    plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
    plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
    plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
    if axes != 'none':
        plt.axis(axes)
    plt.tight_layout()
    plt.savefig(title+'.eps')
    plt.close()

def p_nH_multi(nH,y,ylabel='y',title='something vs log nH',range='none',log='yes',
         legends='none',label=''):
    for i in range(len(y)):
        print(i)
        plt.plot(nH,y[i])
    if log == 'yes':
        plt.yscale('log')
    plt.title(title)
    plt.xlabel(r'$\log(\mathrm{nH}/1$cm$^{-3})$')
    plt.ylabel(ylabel)
    plt.figtext(0.8,0.2,label,backgroundcolor='white')
    #plt.legend()
    plt.savefig(title+'.eps')
    plt.close()

def plot_EWcontours(range_Hden,range_Phi,energy_grid,filename):
    plt.figure()
    CS = plt.contourf(range_Hden, range_Phi, energy_grid,
                      locator=ticker.LogLocator())
    plt.colorbar(CS)
    plt.title(filename+' '+header)
    plt.xlabel('$\log[$n(H) / 1 cm$^{-3}]$')
    plt.ylabel(r'$\log[\Phi$(H) / 1 cm$^{-2}$ s$^{-1}]$',rotation=90)
    pltfile = filename+'.eps'
    plt.savefig(pltfile, format='eps')
    plt.close()

def plot_EWcontours_radius(range_Hden,range_radius,energy_grid,filename):
    plt.figure()
    CS = plt.contourf(range_Hden,range_radius, energy_grid,
                      locator=ticker.LogLocator())
    plt.colorbar(CS)
    plt.title(filename+' '+header)
    plt.xlabel('$\log[$n(H) / 1 cm$^{-3}]$')
    plt.ylabel(r'$\log[R / 1 cm]$',rotation=90)
    plt.axhline(y=m.log10(2.59e15),color='black')
    plt.axhline(y=m.log10(13.11*2.59e15),color='pink',linewidth=5)
    plt.axhline(y=m.log10(2.59e16),color='black')
    plt.axhline(y=m.log10(2.59e17),color='black')
    plt.axhline(y=m.log10(dust_radius()),color='red')
    plt.axis([7,14,max(range_radius),min(range_radius)])
    pltfile = filename+'_radial.eps'
    plt.savefig(pltfile, format='eps')
    plt.close()

def plot_phi_emissivity(range_Phi,emissivities,Hden,line='line',label='test',plotfolder=''):
    plt.plot(range_Phi,emissivities)
    plt.ylabel('emissivity [cm-2 s-1]',rotation=90)
    plt.xlabel(r'$\log[\Phi$(H) / 1 cm$^{-2}$ s$^{-1}]$')
    plt.yscale('log')
    plt.axis((17,24,min(emissivities),max(emissivities)))
    plt.savefig(plotfolder+'/'+line+'_logPhi_emissivity.eps')
    plt.close()

def plot_logr_emissivity(logr,emissivities,Hden,line='line',label='test',plotfolder='',
                         dust_radius=1e100):
    plt.plot(logr,emissivities)
    plt.ylabel('emissivity [cm-2 s-1]',rotation=90)
    plt.xlabel(r'$\log[R $/ 1 cm$]$')
    plt.yscale('log')
    plt.axis((14,18.5,min(emissivities),max(emissivities)))
    plt.axvline(x=m.log10(2.59e15),color='black')
    plt.axvline(x=m.log10(13.11*2.59e15),color='yellow')
    plt.axvline(x=m.log10(2.59e16),color='black')
    plt.axvline(x=m.log10(2.59e17),color='black')
    plt.axvline(x=m.log10(dust_radius),color='red')
    plt.savefig(plotfolder+'/'+line+'_radius_emissivity.eps')
    plt.close()

def plot_phi_relative_emissivity(logPhi,emissivities,Hden,line='line',
                                 label='test',plotfolder=''):
    plt.plot(logPhi,emissivities)
    plt.ylabel('relative emissivity [cm-2 s-1]',rotation=90)
    plt.xlabel(r'$\log[Phi]$')
    plt.yscale('log')
    plt.savefig(plotfolder+'/'+line+'_Phi_rel_emissivity.eps')
    plt.close()

def plot_logr_relative_emissivity(logr,emissivities,Hden,line='line',label='test',
                                  plotfolder='',dust_radius=1e100):
    plt.plot(logr,emissivities)
    plt.ylabel('relative emissivity [cm-2 s-1]',rotation=90)
    plt.xlabel(r'$\log[R $/ 1 cm$]$')
    plt.yscale('log')
    plt.axis((14,18.5,min(emissivities),max(emissivities)))
    plt.axvline(x=m.log10(2.59e15),color='black')
    plt.axvline(x=m.log10(13.11*2.59e15),color='yellow')
    plt.axvline(x=m.log10(2.59e16),color='black')
    plt.axvline(x=m.log10(2.59e17),color='black')
    plt.axvline(x=m.log10(dust_radius),color='red')
    plt.savefig(plotfolder+'/'+line+'_radius_rel_emissivity.eps')
    plt.close()

def plot_logr_nuFnu(logr,nuFnu,Hden,line='line',label='test',plotfolder='',
                    dust_radius=1e100):
    plt.plot(logr,nuFnu)
    plt.ylabel(r'Incident $\nu F_{\nu}$',rotation=90)
    plt.xlabel(r'$\log[R $/ 1 cm$]$')
    plt.yscale('log')
    plt.axis((14,18.5,min(nuFnu),max(nuFnu)))
    plt.axvline(x=m.log10(2.59e15),color='black')
    plt.axvline(x=m.log10(13.11*2.59e15),color='yellow')
    plt.axvline(x=m.log10(2.59e16),color='black')
    plt.axvline(x=m.log10(2.59e17),color='black')
    plt.axvline(x=m.log10(dust_radius),color='red')
    plt.savefig(plotfolder+'/'+line+'_radius_nuFnu.eps')
    plt.close()

def plot_logr_lumin(logr,Lr,Hden,line='line',label='test',plotfolder='',
                    dust_radius=1e100,ion_lumin=10.**44.26):
    plt.plot(logr,Lr)
    plt.yscale('log')
    plt.ylabel('Enclosed luminosity',rotation=90)
    plt.xlabel(r'$\log[R / 1$ cm$]$')
    plt.axvline(x=m.log10(2.59e15),color='black')
    plt.axvline(x=m.log10(2.59e16),color='black')
    plt.axvline(x=m.log10(2.59e17),color='black')
    plt.axvline(x=m.log10(dust_radius),color='red')
    plt.axhline(y=ion_lumin,color='blue')
    plt.savefig(plotfolder+'/'+line+'_cumulative_luminosity.eps')
    plt.close()

def plot_logr_C(logr,Cr,Hden,line='line',label='test',plotfolder='',
                    dust_radius=1e100):
    linr_lightdays = [10.**(r)/2.59e15 for r in logr]
    plt.plot(linr_lightdays,[C/(4.*3.14) for C in Cr])
    plt.ylabel(r'Cumulative solid angle / 4$\pi$',rotation=90)
    plt.xlabel('R (light-days)$')
    plt.axhline(y=1,color='yellow')
    plt.axvline(x=1,color='black')
    plt.axvline(x=10,color='black')
    plt.axvline(x=100,color='black')
    plt.axvline(x=dust_radius/2.59e15,color='red')
    plt.savefig(plotfolder+'/'+line+'_cumulative_C.eps')
    plt.close()

def plot_logr_lumin_allnH(logr,nH_Lr,nHlist,n,line='line',label='test',plotfolder='',
                    dust_radius=1e100,ion_lumin=10.**44.26):
    for i in range(len(nHlist)):
        plt.plot(logr,nH_Lr[i],label='nH = '+str(nHlist[i]))
    plt.xlabel(r'$\log[R / 1$ cm$]$')
    plt.ylabel('Enclosed luminosity',rotation=90)
    plt.yscale('log')
    plt.axvline(x=m.log10(2.59e15),color='black')
    plt.axvline(x=m.log10(2.59e16),color='black')
    plt.axvline(x=m.log10(2.59e17),color='black')
    plt.axvline(x=m.log10(dust_radius),color='red')
    plt.axhline(y=ion_lumin,color='blue')
    plt.legend(loc='lower right')
    plt.savefig(line+'_cumulative_luminosity_allnH_n'+n+'.eps')
    plt.close()

def plot_logr_emissivities_allnH(logr,nH_em,nHlist,n,line='line',label='test',plotfolder='',
                    dust_radius=1e100):
    for i in range(len(nHlist)):
        plt.plot(logr,nH_Lr[i],label='nH = '+str(nHlist[i]))
    plt.xlabel(r'$\log[R / 1$ cm$]$')
    plt.ylabel('Line flux',rotation=90)
    plt.yscale('log')
    plt.axvline(x=m.log10(2.59e15),color='black')
    plt.axvline(x=m.log10(2.59e16),color='black')
    plt.axvline(x=m.log10(2.59e17),color='black')
    plt.axvline(x=m.log10(dust_radius),color='red')
    plt.legend(loc='lower right')
    plt.text(15,0.1*max(nH_Lr[-1]), label, fontsize=20)
    plt.savefig(line+'_cumulative_luminosity_allnH_n'+n+'.eps')
    plt.close()

def plot_rin_lumin(r_in,Lrin,n='nn',line='line',label='test',plotfolder='',
                   ion_lumin=10.**44.26,Hden='Hden'):
    plt.plot(r_in,Lrin)
    plt.xlabel(r'$\log[R_{in} / 1$ cm$]$')
    plt.ylabel('Line luminosity',rotation=90)
    plt.yscale('log')
    plt.axvline(x=m.log10(2.59e15),color='black')
    plt.axvline(x=m.log10(2.59e16),color='black')
    plt.axhline(y=ion_lumin,color='blue')
    plt.savefig(plotfolder+'/'+line+'_Lmax_vs_Rin_n'+n+'_nH'+str(Hden)+'.eps')
    plt.close()

def plot_rout_lumin(r_out,Lrout,n='nn',line='line',label='test',plotfolder='',
                   ion_lumin=10.**44.26,Hden='Hden'):
    plt.plot(r_out,Lrout)
    plt.xlabel(r'$\log[R_{out} / 1$ cm$]$')
    plt.ylabel('Line luminosity',rotation=90)
    plt.yscale('log')
    plt.axvline(x=m.log10(2.59e15),color='black')
    plt.axvline(x=m.log10(2.59e16),color='black')
    plt.axhline(y=ion_lumin,color='blue')
    plt.savefig(plotfolder+'/'+line+'_Lmax_vs_Rout_n'+n+'_nH'+str(Hden)+'.eps')
    plt.close()

def plot_rin_lumin_allnH(r_in,Lrin,nHlist,n='nn',line='line',label='test',
                   ion_lumin=10.**44.26,obslumin=0.1):
    for i in range(len(nHlist)):
        plt.plot([10.**(r)/2.59e15 for r in r_in[i]],Lrin[i],
                 label='log(nH)='+str(nHlist[i]))
    plt.xlabel(r'$R_{in}$ (lightdays)')
    plt.ylabel('Line luminosity',rotation=90)
    plt.yscale('log')
    plt.xscale('log')
    #plt.axhline(y=ion_lumin,color='blue')
    plt.axhline(y=obslumin,color='black')
    plt.legend(loc='lower right')
    plt.title('Total line luminosity within r_out, as function of imposed r_in.')
    plt.figtext(0.8,0.8,label)
    plt.savefig(line+'_Lmax_vs_Rin_n'+n+'_allnH.eps')
    plt.close()

def plot_rout_lumin_allnH(r_out,Lrout,nHlist,n='nn',line='line',label='test',
                   ion_lumin=10.**44.26,obslumin=0.1):
    for i in range(len(nHlist)):
        plt.plot([10.**(r)/2.59e15 for r in r_out[i]],Lrout[i],
                 label='log(nH)='+str(nHlist[i]))
    plt.xlabel(r'$R_{out}$ (lightdays)')
    plt.ylabel('Line luminosity',rotation=90)
    plt.yscale('log')
    plt.xscale('log')
    plt.axhline(y=obslumin,color='black')
    plt.legend(loc='lower right')
    plt.title('Total line luminosity within r_out, as function of imposed r_out.')
    plt.figtext(0.8,0.8,label)
    plt.savefig(line+'_Lmax_vs_Rout_n'+n+'_allnH.eps')
    plt.close()

def plot_lya_civ_allnH(Ncollist,nHlist,L_lya_Ncol,L_civ_Ncol,L_lya,eL_lya,L_civ,
                       eL_civ,r_in=1):
    for i in range(len(Ncollist)):
        plt.plot(L_lya_Ncol[i],L_civ_Ncol[i],label='log(Ncol)='+str(Ncollist[i]))
        for j in range(len(nHlist)):
            plt.text((L_lya_Ncol[0])[j],(L_civ_Ncol[0])[j],nHlist[j],
                     fontsize=plot_fontsize)
    plt.xlabel(r'L(Ly-$\alpha$) (erg s$^{-1}$)',fontsize=plot_fontsize)
    plt.ylabel(r'L(C$IV$) (erg s$^{-1}$)',rotation=90,fontsize=plot_fontsize)
    plt.yscale('log')
    plt.xscale('log')
    plt.axis([1e40,2*1e43,1e37,3e43])
    plt.plot([L_lya,L_lya],[L_civ,L_civ],
                 color='black',marker='*',markersize=20)
    plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
    plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
    plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
    plt.tight_layout()
    plt.savefig('lineL_civ_rin'+str(r_in)+'.eps')
    plt.close()

def plot_lya_hbeta_allnH(Ncollist,nHlist,L_lya_Ncol,L_civ_Ncol,L_lya,eL_lya,L_civ,
                       eL_civ,r_in=1):
    for i in range(len(Ncollist)):
        plt.plot(L_lya_Ncol[i],L_civ_Ncol[i],label='log(Ncol)='+str(Ncollist[i]))
        for j in range(len(nHlist)):
            plt.text((L_lya_Ncol[0])[j],(L_civ_Ncol[0])[j],nHlist[j],
                     fontsize=plot_fontsize)
    plt.xlabel(r'L(Ly-$\alpha$) (erg s$^{-1}$)',fontsize=plot_fontsize)
    plt.ylabel(r'L(H$\beta$) (erg s$^{-1}$)',rotation=90,fontsize=plot_fontsize)
    plt.yscale('log')
    plt.xscale('log')
    plt.axis([1e40,2*1e43,1e37,3e43])
    plt.plot([L_lya,L_lya],[L_civ,L_civ],
                 color='black',marker='*',markersize=20)
    plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
    plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
    plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
    plt.tight_layout()
    plt.savefig('lineL_hbeta_rin'+str(r_in)+'.eps')
    plt.close()

def plot_lya_he4686_allnH(Ncollist,nHlist,L_lya_Ncol,L_civ_Ncol,L_lya,eL_lya,L_civ,
                       eL_civ,r_in=1):
    for i in range(len(Ncollist)):
        plt.plot(L_lya_Ncol[i],L_civ_Ncol[i],label='log(Ncol)='+str(Ncollist[i]))
        for j in range(len(nHlist)):
            plt.text((L_lya_Ncol[0])[j],(L_civ_Ncol[0])[j],nHlist[j],
                     fontsize=plot_fontsize)
    plt.xlabel(r'L(Ly-$\alpha$) (erg s$^{-1}$)',fontsize=plot_fontsize)
    plt.ylabel(r'L(He 4676) (erg s$^{-1}$)',rotation=90,fontsize=plot_fontsize)
    plt.yscale('log')
    plt.xscale('log')
    plt.axis([1e40,2*1e43,1e37,3e43])
    plt.plot([L_lya,L_lya],[L_civ,L_civ],
                 color='black',marker='*',markersize=20)
    plt.gca().tick_params(pad=xtick_padding,width=tickwidth,length=2.5*tickwidth)
    plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
    plt.tick_params(axis='both', which='major', labelsize=plot_ticksize)
    plt.tight_layout()
    plt.savefig('lineL_he4686_rin'+str(r_in)+'.eps')
    plt.close()
