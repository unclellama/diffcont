# test lumin and radius integrator for C iv
from ewgrids import *
import matplotlib.pyplot as plt

r_in = m.log10(1*lightday)
r_out = m.log10(140.*lightday)
nH = 10.
f = 'test_civ_n23/fort.146_n23' # raw data
f_mike = 'test_civ_n23/l_s0_c4_nh10_fort.16' # mike's luminosity curves

# ==================================================================
# load and slice data - this part works!
# ==================================================================
# load data
header,lognH_tab,logPhi_tab,EW_tab = read_fort(f)
# arrange EWs on grid
lognH,logPhi,EWgrid = arrange_grid(lognH_tab,logPhi_tab,EW_tab)
# get a constant-nH slice
x = ((np.where(lognH == nH))[0])[0]
EW_slice = EWgrid[:,x] # checked this code against data files - slices ok.

# ==================================================================
# convert equivalent widths to emissivities
# ==================================================================
# multiply with continuum wavelength to get a relative energy
cont_wl = 1215.
rel_emissivities = [cont_wl*eqwid for eqwid in EW_slice]
# to get incident continuum:
# use the normalization log(nu*F_nu)=10.577 at radius where logPhi=20
incident_nu_F_nu = [10.**(P-10.577) for P in logPhi]
# multiply EW from data with actual incident continuum
emissivities = [incident_nu_F_nu[i]*EW_slice[i] for i in range(len(logPhi))]

# ==================================================================
# convert ionizing cont Phi to radius - this part works! (fluxes match mike)
# ==================================================================
# convert the logPhi to a logr (in cm)
r20_lightdays = 14.81
logr = [-0.5*(P-20)+m.log10(lightday)+m.log10(r20_lightdays) for P in logPhi]
# reverse axes so logr goes from low to high. make np.arrays.
logr = np.asarray(logr[::-1])
emissivities = np.asarray(emissivities[::-1])

# load mike's data for comparison
col1 = []
col2 = []
with open('/Users/danielpe/Dropbox/dark/leicester/data/test_civ_n23/c4_nh10_xs.dat') as f:
    for line in f:
        col1.append(float((line.split())[0]))
        col2.append(float((line.split())[1]))
logr_mike = col1
emissivities_mike = col2

# clip the radius and emissivity arrays to r_in,r_out
logrclip,emclip = clip_radii(logr,emissivities,r_in,r_out)
# linearize radii (clipped to r_in, r_out)
lin_r = np.asarray([10.0**r for r in logrclip])
r_interp = np.arange(10.**(r_in),10.**(r_out),0.1*lightday)
r_interp_log = [m.log10(r) for r in r_interp]
em_log = [m.log10(em+1e-19) for em in emclip] # best to interpolate fluxes in logspace
em_interp_log = np.interp(r_interp,lin_r,em_log)
em_interp = [10.**em for em in em_interp_log]

lin_r_mike = np.asarray([10.0**r for r in logr_mike])
r_interp_mike = np.arange(lin_r_mike[0],lin_r_mike[-1],0.2*lightday)
em_mike_interp = np.interp(r_interp_mike,lin_r_mike,emissivities_mike)
em_mike_interp = [10.**(em+1e-19) for em in em_mike_interp]

# ==================================================================
# get d f(r) / dr, which is used to calculate responsivities
# ==================================================================

dfdr = []
for i in range(len(r_interp)-1):
    df = em_interp_log[i+1]-em_interp_log[i]
    dr = r_interp_log[i+1]-r_interp_log[i]
    dfdr.append(df/dr)
p_ld([r/lightday for r in r_interp[:-1]],dfdr,title='dfdr') # plotting

# ==================================================================
# get responsivities (use that for s=0, d log(Phi) = -2 d log(r) )
# ==================================================================

responsivity = [-0.5*diff for diff in dfdr]
p_ld([r/lightday for r in r_interp[:-1]],responsivity,title='responsivity',
     log='no',range=[-10,2],hlines=[0.,1.])

# ==================================================================
# do the integrals once to get L_tot. *the L(r) shape is ok, scaling may be off!*
# ==================================================================

# get A_c scaling for s=0

scale = Ac_scaling_analytic(r_in,r_out)

# do integral

def intline(rr,em,Ac_scaling=1.):

    dL = [4*m.pi*em[i]*(Ac_scaling)*(rr[i]**(-3./2))*(rr[i]**2.) for i in range(len(rr))]
    Lmax = simps(dL,x=rr)
    # calculate angular coverage
    dC = [4*m.pi*Ac_scaling*(r**(-3./2.)) for r in rr]
    C_sterad = np.trapz(dC,x=rr)
    # calculate lumin-weighted radius
    centroid = np.trapz([dL[i]*rr[i] for i in range(len(rr))],x=rr) / np.trapz(dL,x=rr)
    return [Lmax,C_sterad,centroid,dL]

test_em = [r**1.5 for r in np.arange(1,10,0.1)]
test,C_test,centroid_test,dL = intline(np.arange(1,10,0.1),test_em,Ac_scaling=1)
print('TEST emissivity as r^1.5 : logLtot:',"{0:.2f}".format(m.log10(test)))
test_analytical = 4*m.pi*(((10**3.)/3)-((0.**3.)/3))
print('TEST analytical: logLtot:',"{0:.2f}".format(m.log10(test)))

Ltot,C_sterad,centroid,dL = intline(r_interp,em_interp,Ac_scaling=scale)
print('logLtot:',"{0:.2f}".format(m.log10(Ltot)),' centroid/ld:',
      "{0:.2f}".format(centroid/lightday))

# ==================================================================
# loop over the integrals to get L(r), dL/dr
# ==================================================================

Lr =  [0]
rLr = [0]
Cr =  [0]
dLr = [0]
Lprev = 0

for i in range(len(r_interp)):
    if i>0:
        rcut = r_interp[0:i]
        emcut = em_interp[0:i]
        L,C,junk,junk2 = intline(rcut,emcut,Ac_scaling=scale)
        Lr.append(L)
        dLr.append(L-Lprev)
        Cr.append(C)
        Lprev = L

p_ld([r/lightday for r in r_interp],dLr,title='dL(r)') # plotting

# ==================================================================
# compare with mike's flux curve
# ==================================================================

plt.plot([(10.**r)/lightday for r in logr_mike],emissivities_mike,color='blue')
plt.plot([r/lightday for r in r_interp],[m.log10(em+1e-10) for em in em_interp],
         color='red')
plt.axis((0,25,4,9.5))
plt.title('fluxes (blue:mike, red:me)')
plt.savefig('test_civ_fluxcompare.eps')
plt.close()

# ==================================================================
# compare with mike's luminosity curve
# ==================================================================

col2 = []
col7 = []
with open(f_mike) as f:
    for line in f:
        col2.append(float((line.split())[1]))
        col7.append(float((line.split())[6]))
import matplotlib.pyplot as plt
plt.plot(col2,col7,color='blue')
r_axis = [r/lightday for r in r_interp]
plt.plot(r_axis,[m.log10(L+1e-5) for L in Lr],color='red')
plt.xlabel('radius (lightdays)')
plt.ylabel('enclosed luminosity (erg s-1)',rotation=90)
plt.axis([0,150,35,43])
plt.savefig('test_civ_lumincompare.eps')
plt.close()
