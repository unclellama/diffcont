import numpy as np
from glob import glob
import shutil
import os
import math as m
from scipy.integrate import simps

from ewgrids_plot import *

lightday = 2.59e15

def intline_tau_old(rr,em,s=0,f_factor=['none'],resp=['none'],Ac_scaling=1.,
                    theta_bins=500,
                do_plots='no',label='a label'):
    """ integrates radial functs over time delay tau relative to continuum LOS."""
    from matplotlib import pyplot as plt
    tau_rebinning = 10
    r_binsize = rr[1]-rr[0]
    tau_binsize = tau_rebinning*r_binsize
    tau = np.arange(0,2*max(rr)+r_binsize,tau_binsize) # maximum time delay is 2*r_out
    theta_binsize = m.pi/float(theta_bins)
    theta = np.arange(0.,theta_bins*theta_binsize,theta_binsize)# angle relative to LOS
    theta = [th+0.5*theta_binsize for th in theta[:-1]]
    # integrate luminosity over r (and get dL in each spherical shell)
    dL = [4*m.pi*em[i]*(Ac_scaling)*(rr[i]**(2.*s/3.-3./2.))*(rr[i]**2.)
          for i in range(len(rr))]
    L = np.trapz(dL,x=rr)
    cumdL = cumtrapz(dL,x=rr,initial=0.)
    # integrate lumin-weighted radius
    Lr = np.trapz([dL[i]*rr[i] for i in range(len(rr))],x=rr)
    cent = Lr/L
    cumdLr = cumtrapz([dL[i]*rr[i] for i in range(len(rr))],x=rr,initial=0.)
    # integrate responsivity-weighted radius
    if resp[0] == 'none':
        resp = np.ones_like(rr)
    dretaL = [resp[i]*dL[i]*rr[i] for i in range(len(rr))]
    retaL = np.trapz(dretaL,x=rr)
    detaL = [resp[i]*dL[i] for i in range(len(rr))]
    etaL = np.trapz(detaL,x=rr)
    cumdetaL = cumtrapz(detaL,x=rr,initial=0.)
    respcent = retaL/etaL
    # loop over radial shells, binning into bins of tau for each shell
    dL_tau = np.zeros_like(tau) # remap dL(r) to dL(tau)
    dL_tau_aniso = np.zeros_like(tau) # dL(tau) including anisotropic clouds
    dtauL_tau = np.zeros_like(tau) # tau*dL (integrate to get the centroid)
    dtauL_tau_aniso = np.zeros_like(tau) # tau*dL for aniso clouds
    detaL_tau = np.zeros_like(tau) # tau*dL (integrate to get the centroid)
    detaL_tau_aniso = np.zeros_like(tau) # tau*dL for aniso clouds
    dtauetaL_tau = np.zeros_like(tau) # tau*dL (integrate to get the centroid)
    dtauetaL_tau_aniso = np.zeros_like(tau) # tau*dL for aniso clouds
    for i in range(len(rr)):   # i loops over spherical shells of radius r
        # setup shells to bin
        tau_shell = [rr[i]*(1-m.cos(th)) for th in theta]
        # bin shells
        for j in range(len(tau_shell)):   # j loops over angles to LOS, 0..pi
            k = int(np.round(tau_shell[j]/tau_binsize))
            # factor 0.5 due to: int(sin(theta))=2 for theta=0..pi
            dL_tau[k] = dL_tau[k]+0.5*dL[i]*m.sin(theta[j])*theta_binsize/tau_rebinning
            dtauL_tau[k] = dtauL_tau[k]+0.5*dL[i]*tau_shell[j]*(
                m.sin(theta[j])*theta_binsize/tau_rebinning)
            detaL_tau[k] = detaL_tau[k]+0.5*dL[i]*resp[i]*(
                m.sin(theta[j])*theta_binsize/tau_rebinning)
            dtauetaL_tau[k] = dtauetaL_tau[k]+0.5*dL[i]*resp[i]*tau_shell[j]*(
                m.sin(theta[j])*theta_binsize/tau_rebinning)
            if f_factor[0] != 'none':
                dL_tau_aniso[k] = dL_tau_aniso[k] + theta_binsize/tau_rebinning*(
                    0.5*dL[i]*m.sin(theta[j])*(1-(2*f_factor[i]-1)*m.cos(theta[j])))
                detaL_tau_aniso[k] = detaL_tau_aniso[k] + theta_binsize/tau_rebinning*(
                    0.5*resp[i]*dL[i]*m.sin(theta[j])*(
                        1-(2*f_factor[i]-1)*m.cos(theta[j])))
                dtauL_tau_aniso[k] = dtauL_tau_aniso[k] + theta_binsize/tau_rebinning*(
                    tau_shell[j]*
                    0.5*dL[i]*m.sin(theta[j])*(1-(2*f_factor[i]-1)*m.cos(theta[j])))
                dtauetaL_tau_aniso[k] = dtauetaL_tau_aniso[k] +theta_binsize*(
                    (1/tau_rebinning)*(tau_shell[j]*resp[i]*
                    0.5*dL[i]*m.sin(theta[j])*(1-(2*f_factor[i]-1)*m.cos(theta[j]))))
    cumdL_tau = cumtrapz(dL_tau,x=tau,initial=0.)
    cumdL_tau_aniso = cumtrapz(dL_tau_aniso,x=tau,initial=0.)
    cumdL_resp_tau = cumtrapz(dL_tau,x=tau,initial=0.)
    cumdL_resp_tau_aniso = cumtrapz(dL_tau_aniso,x=tau,initial=0.)
    L_tau = np.trapz(dL_tau,x=tau)
    L_tau_resp = np.trapz(detaL_tau,x=tau)
    L_tau_aniso = np.trapz(dL_tau_aniso,x=tau)
    L_tau_resp_aniso = np.trapz(detaL_tau_aniso,x=tau)
    tauL_tau = np.trapz(dtauL_tau,x=tau)
    tauL_tau_resp = np.trapz(dtauetaL_tau,x=tau)
    tauL_tau_aniso = np.trapz(dtauL_tau_aniso,x=tau)
    tauL_tau_resp_aniso = np.trapz(dtauetaL_tau_aniso,x=tau)
    cent_tau = tauL_tau/L_tau
    cent_tau_resp = tauL_tau_resp/L_tau_resp
    cent_tau_aniso = tauL_tau_aniso/L_tau_aniso
    cent_tau_resp_aniso = tauL_tau_resp_aniso/L_tau_resp_aniso
    if do_plots == 'yes':
        p_ld([t/lightday for t in tau],[dL_tau,dL_tau_aniso,detaL_tau,detaL_tau_aniso],
             title='dL(tau)',
             vlines='none',label=label,legends=['L-iso','L-aniso','R-iso','R-aniso'],
             xtitle=r'$\tau$ (lightdays)',ytitle=r'$dL(\tau)$')
        p_ld([t/lightday for t in tau],[cumdL_tau,cumdL_tau_aniso],title='L(tau)',
             vlines='none',label=label,legends=['isotropic','anisotropic'],
             xtitle=r'$\tau$ (lightdays)',ytitle=r'enclosed $L(\tau)$')
    #print('integral over tau: ',"{0:.2f}".format(m.log10(L_tau)),
    #      'over r:',"{0:.2f}".format(m.log10(L)),'aniso:',
    #      "{0:.2f}".format(m.log10(L_tau_aniso)))
    print('int over tau, L: ',"{0:.2f}".format(m.log10(L_tau)),
          'cent:',"{0:.2f}".format(cent_tau/lightday),'cent_ani:',
          "{0:.2f}".format(cent_tau_aniso/lightday),
          'rcent:',"{0:.2f}".format(cent_tau_resp/lightday),'rcent_ani:',
          "{0:.2f}".format(cent_tau_resp_aniso/lightday))
    return [L_tau,cumdL_tau,cent_tau_aniso,cent_tau_resp_aniso]

def intline_s2(rr,em,resp='none',Ac_scaling=1.):
    """ the integrator to be called by do_line """
    # integrate lumin
    dL = [4*m.pi*em[i]*(Ac_scaling)*(rr[i]**(4./3.-3./2))*(rr[i]**2.)
          for i in range(len(rr))]
    L = np.trapz(dL,x=rr)
    cumdL = cumtrapz(dL,x=rr,initial=0.)
    # integrate angular coverage
    dC = [4*m.pi*Ac_scaling*(r**(4./3.-3./2.)) for r in rr]
    C = np.trapz(dC,x=rr)
    cumdC = cumtrapz(dC,x=rr,initial=0.)
    # integrate lumin-weighted radius
    Lr = np.trapz([dL[i]*rr[i] for i in range(len(rr))],x=rr)
    cent = Lr/L
    cumdLr = cumtrapz([dL[i]*rr[i] for i in range(len(rr))],x=rr,initial=0.)
    # integrate responsivity-weighted radius
    if resp == 'none':
        resp = np.ones_like(rr)
    dretaL = [resp[i]*dL[i]*rr[i] for i in range(len(rr))]
    retaL = np.trapz(dretaL,x=rr)
    detaL = [resp[i]*dL[i] for i in range(len(rr))]
    etaL = np.trapz(detaL,x=rr)
    cumdetaL = cumtrapz(detaL,x=rr,initial=0.)
    respcent = retaL/etaL
    return [L,cumdL,C,cumdC,cent,respcent,cumdetaL]

def responsivities(logr,logf,step=1,do_plots='yes',label='',smooth='none',
                   inwd_ratio='none'):
    """ gets responsivites as defined K+G2004, using that d log(Phi) = -2 d log(r)."""
    step = int(step)
    print(step)
    if step > 1:
        logf_padded = np.append(logf,[0]*step)
        dlogf = [logf_padded[i+step]-logf[i] for i in range(len(logf))]
        dfdr = [dlogf[i]/(logr[step]-logr[0]) for i in range(len(dlogf))]        
    else:
        dlogf = np.append(np.diff(logf),1.)
        dfdr = [dlogf[i]/(logr[1]-logr[0]) for i in range(len(dlogf))]

    resp = [-0.5*diff for diff in dfdr]
    if do_plots == 'yes':
        p_ld(logr,[logf,dfdr,resp],title='derivs',hlines=[0],
             label=label,log='no',legends=['fluxes','dfdr','resp'],
             vlines=[m.log10(lightday),m.log10(10*lightday)],range=[15,18,-10,10])
        p_ld(logr,logf,title='fluxes',hlines=[0],vlines=[m.log10(lightday)],label=label,
             log='no')
        p_ld(logr,dfdr,title='dfdr',label=label,log='no',hlines=[0.],
             vlines=[m.log10(lightday),m.log10(10*lightday)])
        p_ld([(10.**r)/lightday for r in logr],resp,title='responsivity',
             range=[0,140,-2,max(resp)],
            log='no',hlines=[0.,1.],label=label,ytitle='$\eta$')
        if inwd_ratio != 'none':
            norm = max(logf)
            print(norm)
            normlogf = [f/norm for f in logf]
            p_ld([(10.**r)/lightday for r in logr],[inwd_ratio,normlogf],
                 title='inwd_flux',range=[0,140,0,1.],
                 log='no',label=label,ytitle='',
                 legends=['Inwd flux ratio','Norm. tot flux'])
    return resp

def do_line_constnH(n='22',interp='raw',line='lya',Hden=10,cont_wl=1215,method='trapz',
                    r_in=1.,do_plots='yes',do_cumulative='yes'):
    """ get the integrated luminosity for a single line.
    Hden: log(nH/cm-2). cont_wl: the wavelength at which the line energies are scaled to.
    r_in: inner radius of BLR (lightdays).
    note: this provides line lumins for a covering angle of 4pi. scale to data later."""
    r_in = m.log10(r_in*lightday)
    r_out = m.log10(dust_radius())
    # find the input file for this line
    fortfolder = gen_foldername(n,interp)
    range_Phi,slice,rel_emissivities = match_and_slice(fortfolder,
                                                       line,Hden=Hden,cont_wl=1215.)
    # multiply relative emissivities with cont flux
    incident_nu_F_nu = get_incident_nu_Fnu(range_Phi)
    emissivities = [incident_nu_F_nu[i]*rel_emissivities[i]/1215
                    for i in range(len(rel_emissivities))]
    # convert phi to radius
    logr = Phi_to_radius(range_Phi)
    # make diagnostic plots (emissivity vs logPhi, vs logr) before integrating
    if do_plots == 'yes':
        plotfolder = 'n_'+n+'_'+interp+'_plots'
        plot_phi_emissivity(range_Phi,emissivities,Hden,line=line,label=line,
                        plotfolder=plotfolder)
        plot_phi_relative_emissivity(range_Phi,rel_emissivities,Hden,line=line,label=line,
                                 plotfolder=plotfolder)
        plot_logr_emissivity(logr,emissivities,Hden,line=line,label=line,plotfolder=plotfolder,
                         dust_radius=dust_radius())
        plot_logr_relative_emissivity(logr,rel_emissivities,Hden,line=line,label=line,
                                  plotfolder=plotfolder,dust_radius=dust_radius())
        plot_logr_nuFnu(logr,incident_nu_F_nu,Hden,line=line,label=line,plotfolder=plotfolder,
                    dust_radius=dust_radius())
    # reverse so that radius goes [small...large].
    logr = logr[::-1]
    emissivities = emissivities[::-1]
    # get scaling of A_c*n_c by assuming full coverage
    Ac_scaling = get_Ac_scaling(logr,s=0,r_in=r_in,r_out=r_out)
    # integrate over the line out to dust radius
    Lmax,C,linemoment = integrate_line(logr,emissivities,s=0,r_in=r_in,r_out=r_out,
                            method=method,verbose='yes',Ac_scaling=Ac_scaling)
    print('log lumin out to dust radius:',m.log10(Lmax))
    print('covering fraction out to dust radius:',C/(4.*m.pi))
    print('raidus centroid:',linemoment/lightday)
    if do_cumulative == 'yes':
        # make cumulative r vs L(R<r) plot
        Lr = [0]
        Cr = [0]
        print('Log inner radius:',logr[0])
        print('Log r, outer edge of grid:',logr[-1])
        print('Log r, dust radius:',m.log10(dust_radius()))
        i_rmin = ((np.where(logr > r_in))[0])[0]
        i_rmax = ((np.where(logr > m.log10(dust_radius())))[0])[0]-1
        logr_range = logr[i_rmin:i_rmax]
        np.append(logr_range,m.log10(dust_radius()))
        for rmax in logr_range:
            L,C,linemoment2 = integrate_line(logr,emissivities,s=0,r_in=r_in,r_out=rmax,
                            method=method,Ac_scaling=Ac_scaling)
            Lr.append(L)
            Cr.append(C)
        print(Lr)
        if do_plots == 'yes':
            plot_logr_lumin(logr[i_rmin-1:i_rmax],Lr,line=line,label='test',
                            plotfolder=plotfolder,dust_radius=dust_radius(),Hden=Hden)
            plot_logr_C(logr[i_rmin-1:i_rmax],Cr,line=line,label='test',
                        plotfolder=plotfolder,dust_radius=dust_radius(),Hden=Hden)
        return [logr_range,Lr,emissivities[:-1],Lmax,linemoment]
    else:
        return [0,0,0,Lmax]

def integrate_line(logr,emissivity,s=0,r_in=15,r_out=m.log10(dust_radius()),
                   Ac_scaling=1.):
    """ e.g. goad+1993 eqs (9) and (10)
    integrate along the emissivity function, using 4pi = int(k*A_c(r)*n_c(r)/r^2 dr)
    to get the normalization k for A_c*n_c - i.e., assume full coverage (scale L later!).
    r_in, r_out, logr need to be in units of log(r/1cm).
    """
    # clip the radius and emissivity arrays to r_in,r_out
    logrclip,emclip = clip_radii(logr,emissivity,r_in,r_out)
    # linearize radii (clipped to r_in, r_out)
    lin_r = np.asarray([10.0**r for r in logrclip])
    # integrate over liminosity using goad+1993 relations
    dL = [4*m.pi*emclip[i]*(Ac_scaling)*(lin_r[i]**(2.*s/3.-3./2.+2.))
                 for i in range(len(lin_r))]
    Lmax = np.trapz(dL,x=lin_r)
    # calculate angular coverage
    dC = [4*m.pi*Ac_scaling*(r**(2.*s/3.-3./2.)) for r in lin_r]
    C_sterad = np.trapz(dC,x=lin_r)
    # calculate lumin-weighted radius
    centroid = np.trapz([dL[i]*lin_r[i] for i in range(len(lin_r))],
                        x=lin_r)/np.trapz(dL,x=lin_r)
    return [Lmax,C_sterad,centroid,lin_r,dL]

def test_integrate_line(r_in=15.5,r_out=18.5,gridsize=0.01):
    """ tests the line integrator using emissivity=r^2."""
    logr = np.arange(14.,20.,gridsize)
    emissivity = np.asarray([r**2 for r in np.asarray([10.0**r for r in logr])])
    Lmax,C_sterad,centroid,lin_r,dL = integrate_line(logr,emissivity,Ac_scaling=1.,
                                                     verbose='yes',r_in=r_in,r_out=r_out)
    print('L integral from r_in,r_out: ',"{0:.6e}".format(Lmax))
    print('analytical solution: ',"{0:.6e}".format((8.*m.pi/7.)*
                                (((10.**r_out)**(7./2.))-((10.**r_in)**(7./2.)))))
    print('centroid from r_in,r_out: ',"{0:.6e}".format(centroid))
    print('analytical solution: ',"{0:.6e}".format((14./18.)*((((10.**r_out)**(9./2.))-((10.**r_in)**(9./2.)))/(((10.**r_out)**(7./2.))-((10.**r_in)**(7./2.))))))


def do_Lr_line(n='22',interp='raw',line='lya',cont_wl=1215,method='trapz'):
    nH_Lr = []
    nH_em = []
    emwt_radius = []
    nHlist = [8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14]
    for nH in nHlist:
        logr,Lr,em,Lmax = do_line_constnH(n=n,interp=interp,line=line,Hden=nH,
                                          cont_wl=cont_wl,method=method)
        nH_Lr.append(Lr)
        nH_em.append(em)
        linr = [10.**r for r in logr]
        print(len(linr))
        print(len(Lr))
        Lr_prev = [0]+Lr
        emwt_radius.append((np.trapz([linr[i]*(Lr[i]-Lr_prev[i])
                                      for i in range(len(linr))])/Lmax)/2.59e15)
    print('emission-weighted radii:')
    for i in range(len(nHlist)):
        print(nHlist[i],emwt_radius[i])
    plt.plot(nHlist,emwt_radius)
    plt.xlabel('log nH')
    plt.ylabel('Luminosity-weighted radius')
    plt.savefig(line+'_lumin_wt_radius.eps')
    #plot_logr_lumin_allnH(logr,nH_Lr,nHlist,n,line=line,label=line,plotfolder='',
    #                dust_radius=dust_radius())
    #plot_logr_emissivities_allnH(logr,nH_em,nHlist,n,line=line,label=line,plotfolder='',
    #                dust_radius=dust_radius())

def slice_density(raw_Hden,raw_Phi,raw_EW,Hden=10):
    range_Hden,range_Phi,grid = arrange_grid(raw_Hden,raw_Phi,raw_EW)
    x = ((np.where(range_Hden == Hden))[0])[0]
    Hden_slice = grid[:,x]
    return [range_Phi,Hden_slice]



def match_and_slice(fortfolder,line,Hden=10,cont_wl=1215.):
    match = 0
    filelist = grab_files(fortfolder)
    for fortfile in filelist:
        header,raw_Hden,raw_Phi,raw_EW = read_fort(fortfile)
        if header == header_id(line):
            match = match+1
            print('matched header:'+fortfile)
            range_Phi,slice = slice_density(raw_Hden,raw_Phi,raw_EW,Hden = Hden)
            rel_emissivities = [cont_wl*eqwid for eqwid in slice]
    if match > 1:
        print('\n\nfound >1 matches for this line! using last match...\n\n')
    if match < 1:
        raise Exception('failed to match this line!')
    return [range_Phi,slice,rel_emissivities]

def compare_continua(n='22',interp='raw',wl='1000.00A',Hden=10,cont_wl=1215):
    
    # make a filelist
    fortfolder = gen_foldername(n,interp)
    filelist = grab_files(fortfolder)
    for fortfile in filelist:
        header,raw_Hden,raw_Phi,raw_EW = read_fort(fortfile)
        if (wl in header) & ('nFnu' in header):
            range_Phi_nFnu,Hden_slice_nFnu=slice_density(raw_Hden,raw_Phi,raw_EW,
                                                         Hden = Hden)
        if (wl in header) & ('nInu' in header):
            range_Phi_nInu,Hden_slice_nInu=slice_density(raw_Hden,raw_Phi,raw_EW,
                                                         Hden = Hden)
        if (wl in header) & ('InwT' in header):
            range_Phi_InwT,Hden_slice_InwT=slice_density(raw_Hden,raw_Phi,raw_EW,
                                                         Hden = Hden)
        if (wl in header) & ('InwC' in header):
            range_Phi_InwC,Hden_slice_InwC=slice_density(raw_Hden,raw_Phi,raw_EW,
                                                         Hden = Hden)
        if (wl in header) & ('InwD' in header):
            range_Phi_InwD,Hden_slice_InwD=slice_density(raw_Hden,raw_Phi,raw_EW,
                                                         Hden = Hden)
    plt.plot(range_Phi_nFnu,[cont_wl*i for i in Hden_slice_nFnu],label='nFnu')
    plt.plot(range_Phi_nInu,[cont_wl*i for i in Hden_slice_nInu],label='nInu')
    plt.plot(range_Phi_InwT,[cont_wl*i for i in Hden_slice_InwT],label='InwT')
    plt.plot(range_Phi_InwC,[cont_wl*i for i in Hden_slice_InwC],label='InwC')
    plt.plot(range_Phi_InwD,[cont_wl*i for i in Hden_slice_InwD],label='InwD')
    plt.plot(range_Phi_nFnu,[cont_wl*i for i in [Hden_slice_InwC[i]+Hden_slice_InwD[i] for
                         i in range(len(Hden_slice_InwD))]],label='InwC+InwD',linestyle='--')
    plt.ylabel('$W_r$ (scaled by 1215A incident flux)',rotation=90)
    plt.xlabel(r'$\log[\Phi$(H) / 1 cm$^{-2}$ s$^{-1}]$')
    plt.yscale('log')
    plt.axis((17,24,1e-2,1e4))
    plt.title(r'nH=1e'+str(Hden)+'cm$^{-3}$, cont wl='+wl+', ColDens: '+n+'cm$^{-2}$')
    plt.legend(loc='lower left')
    pltfile = 'slice_Hden'+str(Hden)+'_wl'+wl+'_'+n+'.eps'
    plt.savefig(pltfile, format='eps')
    plt.close()

    def run_ewcont(n='22',interp='raw',upto=-1):
    """ makes relative-flux contour plots for each fort.xxx file in a cloudy output."""
    # open a logfile to put the headers in
    logfile = 'n_'+n+'_logs/'+interp+'_logfile'
    fortfolder = gen_foldername(n,interp)
    plotfolder = fortfolder+'_plots/'
    if not os.path.exists(plotfolder):
        os.makedirs(plotfolder)
    f = open(logfile, 'w')
    # go through files, making rel.flux plots
    filelist = grab_files(fortfolder,upto=upto)
    linelist = ['inci1215','inci4860','lya','mgii','hbeta','he1640','he4686','civ']
    for fortfile in filelist:
        # make an EW contour plot for each file
        header = EW_contours(fortfile)
        # copy some of the files to get legible names
        for line in linelist:
            if header == header_id('line'):
                shutil.copy2(fortfile+'.eps', plotfolder+line+'.eps')
                shutil.copy2(fortfile+'_radial.eps', plotfolder+line+'_radial.eps')
        f.write(fortfile+' '+header) # write this header to log
    f.close()
