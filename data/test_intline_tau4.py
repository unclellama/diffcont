from ewgrids import *


# testing

rin = 1.*lightday
rout = 20.*lightday
rres = 0.01*lightday
rr = np.arange(rin,rout+rres,rres)
em = np.zeros_like(rr)
em[-1] = (1.e38/rres)/(rout**2)
#em[-100] = (1./rres)/(rout**2)
ff = [1.]*len(rr)
resp = [1.]*len(rr)

A = Ac_scaling_analytic(m.log10(rin),m.log10(rout),s=0)

L,cumdL,C,cumdC,cent,respcent,cumdetaL = intline(rr,em,resp=resp,Ac_scaling=A,s=0.,
                                                 testing='no')

print('r integral: L:',L,'C:',cent)

L_tau,cumdL,cent_aniso,cent_resp_aniso = intline_tau(rr,em,s=0,f_factor=ff,resp=resp,
    theta_bins=500,Ac_scaling=A,do_plots='yes',label='',testing='no')

print('log L from r:',m.log10(L),'cent:',cent/lightday)
print('log L from tau (in package):',m.log10(L_tau),'cent:',cent/lightday)
