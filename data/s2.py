# attempts at doing the pressure-law BLR s=2 case (constant ionization parameter)

from ewgrids import *
from matplotlib import pyplot as plt

# line on which to test this stuff!

line = 'he4686'
interp = 'raw'

# choose a point to start at in nH-Phi. Ncol goes as r^(-4/3), so choose small Ncol
# for large r, i.e., for small Phi

starting_lognH = 11
starting_logPhi = 19.5
starting_logNcol_float = 22

r_in = m.log10(1.*lightday)
r_out = m.log10(140.*lightday)

interpolation=0.1 # linear interpolation stepsize, in lightdays

# first steps

c_cms = 2.998e+10 # speed of light in cm s-1
starting_logU = starting_logPhi-starting_lognH-m.log10(c_cms)
print('log U:',starting_logU)
starting_logNcol = lookup_Ncol(starting_logNcol_float)

# load the grid containing the starting datapoint

gridfile = match_line(line=line,Ncol=starting_logNcol,interp=interp)
header,lognH,logPhi,raw_EW = read_fort(gridfile)
step = 0.25
lognH,logPhi,EW_grid = arrange_grid(lognH,logPhi,raw_EW,stepsize_x=step,stepsize_y=step)

# get a diagonal slice of the Ncol = start grid, going through starting point
lognH_cut,logPhi_cut,EW_cut,i_cut = diagonal_slice_grid(lognH,logPhi,starting_lognH,
                        starting_logPhi,EW_grid,stepsize_x=step,stepsize_y=step)
# here, i_cut is the index of the starting datapoint with respect to the diagonal cut!

# convert equivalent widths to emissivities

incident_nu_F_nu = get_incident_nu_Fnu(logPhi_cut)
fluxes = [incident_nu_F_nu[i]*EW_cut[i] for i in range(len(logPhi_cut))]

# convert ionizing cont Phi to radius
logr_descending = Phi_to_radius(logPhi_cut)
logr = np.asarray(logr_descending[::-1])

# get column density as function of radius for s=2

kcol = (10.**starting_logNcol_float)/((10.**logr_descending[i_cut])**(-4./3.))
Ncol_cut = [kcol*((10.**r)**(-4./3.)) for r in logr]
Ncol_cut_descending = Ncol_cut[::-1]
#for i in range(len(lognH_cut)):
#    print('lognH:',lognH_cut[i],'logPhi:',logPhi_cut[i],'rel EW:',EW_cut[i],
#          'r_ld:',10.**logr_descending[i]/lightday,'logNcol:',
#          m.log10(Ncol_cut_descending[i]))
#p_grid(lognH,logPhi,lognH_cut,logPhi_cut) # plotting
#p_EWcontours_overlay(lognH,logPhi,EW_grid,lognH_cut,logPhi_cut,Ncol_cut_descending)
range_logr = Phi_to_radius(logPhi)
p_EWcontours_radius_overlay(lognH,range_logr,EW_grid,lognH_cut,logr_descending,
                         Ncol_cut_descending)
# now load grids for each logNcol

Ncol_list = np.asarray([22,22.25,22.5,22.75,23,23.25,23.5,23.75,24])
Ncol_str_list = [lookup_Ncol(N) for N in Ncol_list]
fluxes_Nblock = []
for i in range(len(Ncol_list)):
    gridfile_N = match_line(line=line,Ncol=Ncol_str_list[i],interp=interp)
    header_N,lognH_N,logPhi_N,raw_EW_N = read_fort(gridfile_N)
    lognH_N,logPhi_N,EW_grid_N = arrange_grid(lognH_N,logPhi_N,raw_EW_N,stepsize_x=step,
                                        stepsize_y=step)
    lognH_cut_N,logPhi_cut_N,EW_cut_N,i_cut_N = diagonal_slice_grid(lognH_N,logPhi_N,
                                        starting_lognH,starting_logPhi,EW_grid_N)
    fluxes_N = [incident_nu_F_nu[i]*EW_cut_N[i] for i in range(len(logPhi_cut))]
    fluxes_Nblock.append(fluxes_N)

# get flux of 'closest' Ncol for each r

nearest_Ncol_cut = np.asarray([takeClosest(Ncol_list,m.log10(N)) for N in Ncol_cut])
best_fluxes = []
interpolated_fluxes = []
interpolate_Ncol = 'yes'

for i in range(len(logr)):
    farr = [(fluxes_Nblock[j])[i] for j in range(len(Ncol_list))]
    interp_flux = np.interp([m.log10(Ncol_cut[i])],Ncol_list,farr)
    interpolated_fluxes.append((interp_flux)[0])
    #print(10.**logr[i]/lightday,Ncol_cut[i],interp_flux[0],
    #          interp_flux[0]-(fluxes_Nblock[iblock])[i])
    iblock = (np.where(Ncol_list == takeClosest(Ncol_list,m.log10(Ncol_cut[i])))[0])[0]
    best_fluxes.append((fluxes_Nblock[iblock])[i])
interp_diff = [(interpolated_fluxes[i]-best_fluxes[i])/interpolated_fluxes[i]
               for i in range(len(best_fluxes))]

for ff in fluxes_Nblock:
    plt.plot([10.**r/lightday for r in logr],ff)
plt.plot([10.**r/lightday for r in logr],best_fluxes,'ko')
plt.xlabel('radius, lighdays')
plt.ylabel('emissivity')
plt.yscale('log')
plt.axis([0,150,1e8,5e9])
plt.show()
plt.close()

p_ld([10.**r/lightday for r in logr],fluxes_Nblock,title='fluxes_Ncol_interp',
         log='yes',axes='none',hlines='none',label='interp over Ncol',
         legends=Ncol_str_list)

# get responsivities, based on flux vs radius.
fine_logr = np.arange(min(logr),max(logr),0.0001)
fine_logf = np.interp(fine_logr,logr,[m.log10(f+1e-19) for f in best_fluxes])
fine_inwd = np.interp(fine_logr,logr,[m.log10(f+1e-19) for f in fluxes]) # dummy
resp = responsivities_smooth(fine_logr,fine_logf,smoothing=1000,
                          do_plots='yes',label=typeline(line))
# Clip the radius and emissivity arrays to r_in,r_out
logrclip,fclip = clip_radii(logr,best_fluxes,r_in,r_out)
# linearize radii and interpolate to get even spacing
linr,linf,logr,logf = interpolate_fluxes(logrclip,fclip,interpolation)
resp = np.interp(linr,[10.**r for r in fine_logr],resp)
# get A_c scaling for s=0
scale = Ac_scaling_analytic(r_in,r_out,s=2)
# do integral
L,cumdL,C,cumdC,cent,respcent,cumdetaL = intline(linr,linf,resp=resp,Ac_scaling=scale,
                                                 s=2.)
dL = [L/interpolation for L in np.insert(np.diff(cumdL),0,0.)]
detaL = [L/interpolation for L in np.insert(np.diff(cumdetaL),0,0.)]
dC = [cov/interpolation for cov in np.insert(np.diff(cumdL),0,0.)]
print('logLtot:',"{0:.2f}".format(m.log10(L)),' centroid/ld:',
        "{0:.2f}".format(cent/lightday),' resp_cent/ld:',"{0:.2f}".format(
        respcent/lightday))
p_ld([r/lightday for r in linr],[m.log10(L+1e-19) for L in cumdL],
         title='r_L(r)',label=typeline(line),log='no',
         ytitle=r'Enclosed $\log[L / $erg s$^{-1}$]',axes=[0,140,38,46])
p_ld([r/lightday for r in linr],dL,title='r_dL(r)',label=typeline(line))
p_ld([r/lightday for r in linr],detaL,title='r_detaL(r)',label=typeline(line),log='no')
