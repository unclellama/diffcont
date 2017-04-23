from ewgrids import *

def intline_tau(rr,em,f_factor='none',resp='none',Ac_scaling=1.,theta_bins=20):
    """ integrates radial functs over time delay tau for a given line-of-sight."""
    from matplotlib import pyplot as plt
    tau_rebinning = 100
    r_binsize = rr[1]-rr[0]
    tau_binsize = tau_rebinning*r_binsize
    tau = np.arange(0,2*max(rr)+1,tau_binsize) # maximum time delay is 2*r_out
    theta_binsize = m.pi/float(theta_bins)
    theta = np.arange(0.,theta_bins*theta_binsize,theta_binsize)# angle relative to LOS
    theta = [th+0.5*theta_binsize for th in theta[:-1]]
    # set up the r-integral (to get the luminosity in each radial shell)
    dL = [4*m.pi*em[i]*(Ac_scaling)*(rr[i]**(-3./2))*(rr[i]**2.)
          for i in range(len(rr))]
    L = np.trapz(dL,x=rr)
    cumdL = cumtrapz(dL,x=rr,initial=0.)
    # loop over radial shells, binning into bins of tau for each shell
    dL_tau = np.zeros_like(tau) # remap dL(r) to dL(tau)
    dL_tau_aniso = np.zeros_like(tau)
    for i in range(len(rr)):   # i loops over spherical shells of radius r
        # setup shells to bin
        tau_shell = [rr[i]*(1-m.cos(th)) for th in theta]
        # bin shells
        for j in range(len(tau_shell)):   # j loops over angles to LOS, 0..pi
            k = int(np.round(tau_shell[j]/tau_binsize))
            # factor 0.5 due to: int(sin(theta))=2 for theta=0..pi
            dL_tau[k] = dL_tau[k]+0.5*dL[i]*m.sin(theta[j])*theta_binsize/tau_rebinning
            if f_factor != 'none':
                dL_tau_aniso[k] = dL_tau_aniso[k] + theta_binsize/tau_rebinning*(
                    0.5*dL[i]*m.sin(theta[j])*(1-(2*f_factor[i]-1)*m.cos(theta[j])))
    cumdL_tau = cumtrapz(dL_tau,x=tau,initial=0.)
    cumdL_tau_aniso = cumtrapz(dL_tau_aniso,x=tau,initial=0.)
    L_tau = np.trapz(dL_tau,x=tau)
    L_tau_aniso = np.trapz(dL_tau_aniso,x=tau)
    plt.plot(tau,dL_tau)
    plt.xlabel(r'$\tau$',fontsize=20)
    plt.ylabel(r'$dL(\tau)$',fontsize=20)
    plt.show()
    plt.close()
    plt.plot(tau,dL_tau_aniso)
    plt.xlabel(r'$\tau$',fontsize=20)
    plt.ylabel(r'$dL(\tau) (\mathrm{anisotropic})$',fontsize=20)
    plt.show()
    plt.close()
    print('integral over tau: ',"{0:.2f}".format(L_tau),
          'over r:',"{0:.2f}".format(L),'over tau (aniso):',
          "{0:.2f}".format(L_tau_aniso))
    return [L,cumdL]

tb = 500

# test the integration over time delay (tau)

rr = np.arange(50.,100.,0.01)
em = [1. for r in rr]
f_factor = [0.5 for r in rr]

L,cumdL = intline_tau(rr,em,f_factor=f_factor,resp='none',Ac_scaling=1.,theta_bins=tb)
