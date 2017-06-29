# codes to work with the NGC 5548 photoionization grids supplied by Otho Adam Ulrich

import numpy as np
import os
import sys
import math as m 
from glob import glob                # glob file search patterns
import shutil                        # move, rename files
from scipy.integrate import simps    # simpson's rule integrator
from scipy.integrate import cumtrapz # cumulative sum
import json                          # to save/load dictionaries as json data
import pexpect                       # allows feeding commands into shell scripts
from astropy.io import ascii         # working with ascii tables
from ewgrids_plot import *           # plotting codes for this project
import matplotlib as mpl             # tells mpl to use local latex install for encoding - allows amsmath \textsc
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

lightday = 2.59e15 # cm
c_cms = 2.998e+10 # cm s-1

# =======================================================================================
# utility
# =======================================================================================

def normalize(v):
    norm=np.linalg.norm(v)
    if norm==0: 
       return v
    return v/norm

def find_nearest(array,value):
    """ returns index and array value of entry nearest to value """
    i = (np.abs(array-value)).argmin()
    return [i,array[i]]

from bisect import bisect_left

def takeClosest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.
    If two numbers are equally close, return the smallest number.
    code stolen from lauritz v thaulow
    http://stackoverflow.com/questions/12141150/from-list-of-integers-get
    -number-closest-to-a-given-value
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return after
    else:
       return before

def takeBounding(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.
    If two numbers are equally close, return the smallest number.
    code stolen from lauritz v thaulow
    http://stackoverflow.com/questions/12141150/from-list-of-integers-get
    -number-closest-to-a-given-value
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    return [before,after]

# =======================================================================================
# file/header reading and sorting
# =======================================================================================

def foldername(n,interp):
    return 'n_'+n+'_'+interp

def grab_files(fortfolder,upto=-1):
    """ find all fort.xxx files, in numerical order of the extension"""
    filelist = glob(fortfolder+'/fort.?')
    filelist.extend(glob(fortfolder+'/fort.??'))
    filelist.extend(glob(fortfolder+'/fort.???'))
    filelist.extend(glob(fortfolder+'/fort.????'))
    if upto == -1:
        return filelist
    else:
        return filelist[:upto]

def read_fort(filename):
    """ reads the data and header from a single fort.xx file.
    assumes that the header dontains the line 'relative to' - shitty coding!
    the files should be tables of nH | Phi(H) | EW, as supplied for this project.
    """
    if os.path.isfile(filename):
        Hden = []
        Phi = []
        EW = []
        with open(filename) as f:
            for line in f:
                if "relative to" in line:
                    header = line
                elif 'Hden   Phi(H)  Eq_Width (A)' in line:
                    pass
                else:
                    Hden.append(float((line.split())[0]))
                    Phi.append(float((line.split())[1]))
                    EWcheck = float((line.split())[2])
                    if EWcheck > 2.0000e-07: # as the data has minimum value 1e-7...
                        EW.append(float((line.split())[2]))
                    else:
                        EW.append(0.) # ...zero low vals to avoid runaway fluxes at high nu
            return [header,np.array(Hden),np.array(Phi),np.array(EW)]
    else:
        raise Exception("no file found")

def header_id(argument):
    switcher = {
        'inci1215' : 'Inci 1215.00A  relative to Inci 1215.00A scaled to 1215A.\n',
        'inci4860' : 'Inci 4860.00A  relative to Inci 1215.00A scaled to 1215A.\n',
        'lya': 'TOTL 1215.68A  relative to Inci 1215.00A scaled to 1215A.\n',
        'lya_inwd': 'Inwd 1215.68A  relative to Inci 1215.00A scaled to 1215A.\n',
        'mgii': 'TOTL 2798.00A  relative to Inci 1215.00A scaled to 1215A.\n',
        'mgii_inwd': 'Inwd 2798.00A  relative to Inci 1215.00A scaled to 1215A.\n',
        'halpha': 'H  1 6562.85A  relative to Inci 1215.00A scaled to 1215A.\n',
        'halpha_inwd': 'Inwd 6562.85A  relative to Inci 1215.00A scaled to 1215A.\n',
        'hbeta': 'H  1 4861.36A  relative to Inci 1215.00A scaled to 1215A.\n',
        'hbeta_inwd': 'Inwd 4861.36A  relative to Inci 1215.00A scaled to 1215A.\n',
        'he1640': 'He 2 1640.00A  relative to Inci 1215.00A scaled to 1215A.\n',
        'he1640_inwd': 'Inwd 1640.00A  relative to Inci 1215.00A scaled to 1215A.\n',
        'he4686': 'He 2 4686.01A  relative to Inci 1215.00A scaled to 1215A.\n',
        'he4686_inwd': 'Inwd 4686.01A  relative to Inci 1215.00A scaled to 1215A.\n',
        'civ': 'TOTL 1549.00A  relative to Inci 1215.00A scaled to 1215A.\n',
        'civ_inwd': 'Inwd 1549.00A  relative to Inci 1215.00A scaled to 1215A.\n'
    }
    return switcher.get(argument, "nothing")

def typeline(argument):
    switcher = {
        'inci1215' : 'Incident flux 1215',
        'inci4860' : 'Incident flux 4860',
        'lya': r'$\mathrm{Ly}$-$\alpha$',
        'mgii': r'$\mathrm{Mg}\,\textsc{ii}$',
        'halpha': r'$\mathrm{H}\alpha$',
        'hbeta': r'$\mathrm{H}\beta$',
        'he1640': r'$\mathrm{He}\,\textsc{ii}\,1640$',
        'he4686': r'$\mathrm{He}\,\textsc{ii}\,4686$',
        'civ': r'$\mathrm{C}\,\textsc{iv}$',
        'continuum': 'Diff. cont.'
    }
    return switcher.get(argument, "nothing")

def match_line(line='lya',Ncol='23',interp='raw'):
    folder = foldername(Ncol,interp)
    match = 0
    filelist = grab_files(folder)
    for file in filelist:
        header,raw_Hden,raw_Phi,raw_EW = read_fort(file)
        if header == header_id(line):
            return file
    raise Exception('failed to match header: '+header_id(line)+' in folder: ',folder)

def match_cont(contwl,Ncol='23',interp='raw',key='nFnu'):
    folder = foldername(Ncol,interp)
    filelist = grab_files(folder)
    for file in filelist:
        header,raw_Hden,raw_Phi,raw_EW = read_fort(file)
        if key in header:
            wl = float(get_wl(header))
            if ((wl < contwl+1) & (wl > contwl-1)):
                return file
    raise Exception('couldnt find continuum file.')

def get_wl(header):
    """ extract the line wavelength from header"""
    line_wl = float(header[5:11])
    line_unit = header[12]
    if line_unit == 'm':
        line_wl = line_wl*10000.
    return line_wl

def diff_cont_files(Ncol='23',interp='raw',minwl=1000.,maxwl=10000.,key='nFnu'):
    folder = foldername(Ncol,interp)
    filelist = grab_files(folder)
    contfiles = []
    contwl = []
    for file in filelist:
        header,raw_Hden,raw_Phi,raw_EW = read_fort(file)
        if key in header:
            wl = float(get_wl(header))
            if ((wl < maxwl) & (wl > minwl)):
                contwl.append(wl)
                contfiles.append(file)
    if len(contwl) == 0:
        raise Exception('couldnt find any continuum files.')
    else:
        return [contwl,contfiles]

def lookup_Ncol(argument):
    switcher = {
        22 : '22',
        22.25 : '22p25',
        22.5 : '22p5',
        22.75 : '23p75',
        23 : '23',
        23.25: '23p25',
        23.5 : '23p5',
        23.75 : '23p75',
        24 : '24',
    }
    return switcher.get(argument, "nothing")

# =======================================================================================
# measurements / cosmology
# =======================================================================================

def DL(DL_Mpc=72.5):
    # default DL_Mpc value from misty bentz MBH database
    Mpc = 3.086e+24 # 1 Mpc in cm
    return DL_Mpc*Mpc

def calc_line_lumin(F_int=5*1e-13,wl=1215,z=0.01717):
    """ convert integrated line flux (ergs s-1 cm-2) to line lumin ergs s-1"""
    Mpc = 3.086e+24 # 1 Mpc in cm
    DL_cm = DL()
    return 4*m.pi*(DL_cm**2)*(1+z)*F_int

def calc_lambda_L(F_lambda=42.64*1e-15,wl=1367,beta=-1.5):
    """ calculates the radius (in light-days) for which log(Phi(H))=20.
    [F_lambda] = ergs s-1 cm-2 Ang-1
    default F_lambda_1367 is from deRosa+2015"""
    DL_cm = DL()
    F_lambda_1215 = F_lambda*((1215/wl)**(-1.5))
    return 4*m.pi*(DL_cm**2)*F_lambda_1215*1215

def observed_lumin(line,z=0.01717):
    # observations
    F_he1640 = 7.93*1e-13
    eF_he1640 = 1.13*1e-13
    F_lya = 41.22*1e-13  # deRosa+15
    eF_lya = 8.60*1e-13
    F_civ = 53.28*1e-13 # deRosa+15
    eF_civ = 3.91*1e-13
    F_hbeta = 738.49*1e-15 # pei+17
    eF_hbeta = 28.29*1e-15
    F_he4686 = 78.71*1e-15 # pei+17
    eF_he4686 = 35.14*1e-15
    # narrow-line subtraction
    F_lya = F_lya -89.5*1e-14 # K+G2000
    F_civ = F_civ-6.97*1e-13 # K+G2000
    F_hbeta = F_hbeta-8.43*1e-14 # peterson+1999
    # dereddening
    F_lya = F_lya*1.3140856 # seaton1979, howarth1983
                            # (E(B-V)=0.03 lookup table supplied by MG)
    F_civ = F_civ*1.2498338
    F_hbeta = F_hbeta*1.1055139
    F_he4686 = F_he4686*1.1102941

    # line fluxes to luminosities
    if line == 'lya':
        return [calc_line_lumin(F_int=F_lya,wl=1215,z=0.01717),
                calc_line_lumin(F_int=eF_lya,wl=1215,z=0.01717)]
    if line == 'civ':
        return [calc_line_lumin(F_int=F_civ,wl=1550,z=0.01717),
                calc_line_lumin(F_int=eF_civ,wl=1551,z=0.01717)]
    if line == 'hbeta':
        return [calc_line_lumin(F_int=F_hbeta,wl=4861,z=0.01717),
            calc_line_lumin(F_int=eF_hbeta,wl=4861,z=0.01717)]
    if line == 'he1640':
        return [calc_line_lumin(F_int=F_he1640,wl=1640,z=0.01717),
                calc_line_lumin(F_int=eF_he1640,wl=1640,z=0.01717)]
    if line == 'he4686':
        return [calc_line_lumin(F_int=F_he4686,wl=4861,z=0.01717),
                calc_line_lumin(F_int=eF_he4686,wl=4861,z=0.01717)]
    else:
        raise Exception('line not defined!')

# =======================================================================================
# manipulating EW grids
# =======================================================================================

def arrange_grid(raw_x,raw_y,raw_z,stepsize_x=0.25,stepsize_y=0.25):
    """ given list data of type [x,y,z], arranges the z values on an x vs y grid.
    the stepsize needs to be constant, and needs to be specified - the latter i should
    probably automize later.
    """
    range_x = np.arange(np.min(raw_x),np.max(raw_x)+stepsize_x,stepsize_x)
    range_y = np.arange(np.min(raw_y),np.max(raw_y)+stepsize_y,stepsize_y)
    grid = np.empty([len(range_x),len(range_y)])
    for x in range(len(range_x)):
        for y in range(len(range_y)):
            grid[x,y] = raw_z[np.where((raw_x == range_x[x]) & (raw_y == range_y[y]))]
    grid = np.transpose(grid)
    return range_x,range_y,grid

def EW_contours(filename,cont_wl=1215):
    ''' plots the equivalent width EW(Hden,Phi) as an intensity map'''
    # get data+header, arrange EWs on grid
    header,raw_Hden,raw_Phi,raw_EW = read_fort(filename)
    line_wl = get_wl(header)
    range_Hden,range_Phi,EW_grid = arrange_grid(raw_Hden,raw_Phi,raw_EW)
    # multiply with 1215 Ang
    energy_grid = np.multiply(cont_wl,EW_grid)
    # plot Hden, Phi, EW
    plot_EWcontours(range_Hden,range_Phi,energy_grid,filename)
    # plot Hden, log(R), EW
    plot_EWcontours_radius(range_Hden,Phi_to_radius(range_Phi),energy_grid,filename)
    return header

def diagonal_slice_grid(lognH,logPhi,starting_lognH,starting_logPhi,EW_grid,
                        stepsize_x=0.25,stepsize_y=0.25):
    """ makes a diagonal (1:1 logphi-lognh) slice through an EW grid, going through
    a given point. for the s=2 pressure law model."""
    # find coords of the starting point
    i_lognH = ((np.where(lognH == starting_lognH))[0])[0]
    i_logPhi = ((np.where(logPhi == starting_logPhi))[0])[0]
    # find the min, max grid values for diagonal line through this point
    min_lognH = starting_lognH-min(i_logPhi,i_lognH)*stepsize_x
    min_logPhi = starting_logPhi-min(i_logPhi,i_lognH)*stepsize_y
    max_lognH = starting_lognH+(len(lognH)-1-max(i_logPhi,i_lognH))*stepsize_x
    max_logPhi = starting_logPhi+(len(logPhi)-1-max(i_logPhi,i_lognH))*stepsize_x
    i_min_lognH = ((np.where(lognH == min_lognH))[0])[0]
    i_min_logPhi = ((np.where(logPhi == min_logPhi))[0])[0]
    i_max_lognH = ((np.where(lognH == max_lognH))[0])[0]
    i_max_logPhi = ((np.where(logPhi == max_logPhi))[0])[0]
    # get ii and value arrays for logPhi, lognH along the cut
    lognH_cut = np.arange(min_lognH,max_lognH+stepsize_x,stepsize_x)
    logPhi_cut = np.arange(min_logPhi,max_logPhi+stepsize_y,stepsize_y)
    ii_lognH_cut = range(i_min_lognH,i_max_lognH+1)
    ii_logPhi_cut = range(i_min_logPhi,i_max_logPhi+1)
    test_lognH_cut = np.asarray([lognH[i] for i in ii_lognH_cut])
    test_logPhi_cut = np.asarray([logPhi[i] for i in ii_logPhi_cut])
    if not (test_lognH_cut[0] == lognH_cut[0]) & (test_lognH_cut[-1] == lognH_cut[-1]):
        raise Exception('you screwed up with the grid coordinates!')
    EW_cut = [EW_grid[ii_logPhi_cut[i],ii_lognH_cut[i]] for i in range(len(lognH_cut))]
    return [lognH_cut,logPhi_cut,EW_cut,i_logPhi]

# =======================================================================================
# radii; inner and outer BLR radii; logPhi to logr
# =======================================================================================

def get_incident_nu_Fnu(range_logPhi,nu_Fnu_Phi20=10.577):
    """ given the incident nu*F_nu(1215) corresponding to logPhi=20, spits out an
    array of nu*F_nu. the constant depends on the SED shape, and was provided by
    mike goad. """
    return [10.**(logPhi-nu_Fnu_Phi20) for logPhi in range_logPhi]

def Phi_to_radius(logPhi,r20_lightdays = 14.81):
    """convert Phi(H), ionizing flux per unit area, to radius from source.
    r20_lightdays is the radius at which logPhi=20."""
    logr = [-0.5*(P-20)+m.log10(lightday)+m.log10(r20_lightdays)
            for P in logPhi]
    return np.asarray(logr)

def dust_radius(L=10**44.26):
    """ returns dust sublimation radius in cm. default L taken from K+G 2001."""
    #pc_cm = 3.086e+18 # 1 pc in cm
    #return 0.4*((L*1e-45)**0.5)*pc_cm
    return 140.*lightday

def clip_radii(logr,emissivity,r_in,r_out):
    """ needs min(logr)<r_in and max(logr)>r_out.
    clips logr and emissivity arrays to [r_in...r_out], [F(rin)...F(r_out)].
    does linear interpolation in log-log space. give r_in, out in log(cm) ! """
    if not logr[-1] == r_out:
        try:
            ((np.where(logr > r_out))[0])[0]
            ((np.where(logr < r_in))[0])[0]
        except:
            raise Exception(
                'cannot make radial cuts - check that logr covers r_in, r_out')
    em_in = np.interp(r_in, logr, emissivity)
    em_out = np.interp(r_out, logr, emissivity)
    i_in = min((np.where(logr > r_in))[0])
    if not logr[-1] == r_out:
        i_out = max((np.where(logr < r_out))[0])
        lrclip = logr[i_in:i_out+1] # weird python slicing syntax is weird.
        emclip = emissivity[i_in:i_out+1]
    else:
        lrclip = logr[i_in:]
        emclip = emissivity[i_in:]
    lrclip = np.insert(lrclip, 0, r_in)
    emclip = np.insert(emclip, 0, em_in)
    if not logr[-1] == r_out:
        lrclip = np.append(lrclip,r_out)
        emclip = np.append(emclip,em_out)
    return [lrclip,emclip]

# =======================================================================================
# scaling of A_c*n_c integral to achieve total coverage at r_out
# =======================================================================================

def Ac_scaling_analytic(r_in,r_out,s=0):
    """ provides the scaling factor k to get C=4pi when integrating
    dC = k*r^(2s/3-3/2) over the volume r_in...r_out"""
    rlin_in = 10.**r_in
    rlin_out = 10.**r_out
    if s == 0:
        inverse = ((2./m.sqrt(rlin_in))-(2./m.sqrt(rlin_out)))
        return 1./inverse
    if s == 2:
        inverse = (6./5.)*((rlin_out**(5./6.))-(rlin_in**(5./6.)))
        return 1./inverse
    else:
        print('this s not implemented!')

def Ac_scaling_numerical(logr,s=0,r_in=15,r_out=m.log10(dust_radius()),verbose='no'):
    emissivity = np.asarray(logr) # dummy values - don't need emissivities here
    # clip input to r_in, r_out
    logrclip,emclip = clip_radii(logr,emissivity,r_in,r_out)
    # linearize radii
    lin_r = np.asarray([10.0**r for r in logrclip])
    # normalize by angular coverage constraint
    scaling = 1.e8 # scale numerically
    dC = [4*m.pi*scaling*(r**(2.*s/3.-3./2.)) for r in lin_r]
    C = np.trapz(dC,x=lin_r)
    scorr = scaling/(C/(4*m.pi))
    if verbose == 'yes':
        print('scaling:',scorr)
        dCcorr = [4*m.pi*(scorr)*(r**(-3./2.)) for r in lin_r]
        Ccorr = np.trapz(dCcorr,x=lin_r)
        print('angular coverage/4pi at r_out:',Ccorr/(4.*m.pi))
    return scorr

# =======================================================================================
# luminosity (etc) integrators
# =======================================================================================

def interpolate_fluxes(lr,f,interpolation):
    r_interp = np.arange(10.**(lr[0]),10.**(lr[-1]),interpolation*lightday)
    r_interp_log = [m.log10(r) for r in r_interp]
    f_log = [m.log10(flux+1e-19) for flux in f] # best to interpolate fluxes in logspace
    f_interp_log = np.interp(r_interp_log,lr,f_log)
    f_interp = [10.**flux for flux in f_interp_log]
    return [r_interp,f_interp,r_interp_log,f_interp_log]

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def responsivities_smooth(logr,logf,do_plots='yes',label='',step=10,smoothing=10):
    """ gets responsivites as defined K+G2004, using that d log(Phi) = -2 d log(r)."""
    # smooth emissivity function
    logf_smooth = smooth(logf,smoothing)
    # step up in radii
    logf_padded = np.append(logf_smooth,[0]*step)
    dlogf = [logf_padded[i+step]-logf_padded[i] for i in range(len(logf))]
    dfdr_up = [dlogf[i]/(logr[step]-logr[0]) for i in range(len(dlogf))]
    #step down in radii
    logf_padded = np.append([0]*step,logf_smooth)
    dlogf = [logf_padded[i]-logf_padded[i-step] for i in [j+step for j in range(len(logf))]]
    dfdr_down = [dlogf[i]/(logr[step]-logr[0]) for i in range(len(dlogf))]
    dfdr_down = np.append([0]*step,dfdr_down)
    # average the up and down gradients
    dfdr = [(dfdr_up[i]+dfdr_down[i])/2. for i in range(len(dfdr_up))]
    # responsivity = 1/2dFdr
    resp = [-0.5*diff for diff in dfdr]
    resp_smooth = smooth(resp,smoothing)
    if do_plots == 'yes':
        p_ld([(10**r)/lightday for r in logr],[logf,logf_smooth,dfdr,resp],
             title='derivs',hlines=[0],
             label=label,log='no',legends=['log(F)','l(F)smooth','dfdr','resp'],
             vlines=[m.log10(lightday),m.log10(140*lightday)],
             axes=[0,140,-10,10],label_xy=[0.5,0.8])
        p_ld([(10**r)/lightday for r in logr],logf,title='fluxes',hlines=[0],
             vlines=[m.log10(lightday),m.log10(140*lightday)],label=label,
             log='no',ytitle=r'$\log[\epsilon/1$erg cm$^{-2}$s$^{-1}]$')
        p_ld([(10.**r)/lightday for r in logr],resp,title='responsivity',
             axes=[0,140,-10,2],
            log='no',hlines=[0.,1.],label=label,ytitle='$\eta$')
    return resp

def intline(rr,em,resp=['none'],Ac_scaling=1.,s=0.,testing='no'):
    """ the integrator to be called by do_line """
    # integrate lumin
    dL = [4*m.pi*em[i]*(Ac_scaling)*(rr[i]**(2.*s/3.-3./2))*(rr[i]**2.)
          for i in range(len(rr))]
    if testing == 'yes':
        dL = [4*m.pi*em[i]*(Ac_scaling)*(rr[i]**2.) for i in range(len(rr))]        
    L = np.trapz(dL,x=rr)
    cumdL = cumtrapz(dL,x=rr,initial=0.)
    # integrate angular coverage
    dC = [4*m.pi*Ac_scaling*(r**(2.*s/3.-3./2.)) for r in rr]
    C = np.trapz(dC,x=rr)
    cumdC = cumtrapz(dC,x=rr,initial=0.)
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
    return [L,cumdL,C,cumdC,cent,respcent,cumdetaL]

def intline_tau(rr,em,s=0,f_factor=['none'],resp=['none'],Ac_scaling=1.,
                theta_bins=500,do_plots='no',label='a label',testing='no'):
    """ integrates radial functs over time delay tau relative to continuum LOS."""
    from matplotlib import pyplot as plt
    r_binsize = rr[1]-rr[0]
    costheta_binsize = 2./theta_bins
    costheta = np.arange(-1.,1+costheta_binsize,costheta_binsize) # angle rel. to LOS
    costheta = costheta[::-1]
    tau = [rr[-1]*(1-th) for th in costheta]
    tau_binsize = tau[1]-tau[0]
    # loop over radial shells, binning into bins of tau for each shell
    # dLtheta is the dL for a 'ring' of the shell at a given cos(theta).
    dLtheta = [2*m.pi*em[i]*(Ac_scaling)*(rr[i]**(2.*s/3.-3./2.))*(rr[i]**2.)
          for i in range(len(rr))]
    if testing == 'yes':
        dLtheta = [2*m.pi*em[i]*(Ac_scaling)*(rr[i]**2.)
          for i in range(len(rr))]
    dL = np.zeros_like(tau) # remap dL(r) to dL(tau)
    dL_resp = np.zeros_like(tau)
    dL_an = np.zeros_like(tau)
    dL_resp_an = np.zeros_like(tau)
    f_rebin = costheta_binsize*(r_binsize/tau_binsize)
    for i in range(len(rr)):   # i loops over spherical shells of radius r
        # setup spherical shells
        tau_shell = [rr[i]*(1-th) for th in costheta]
        # bin over spherical shells
        for j in range(len(tau_shell)):   # j loops over angles to LOS, 0..pi
            f_trap = 0.5 if (j == 0 or j == len(costheta)-1) else 1.
            k = int(np.round(tau_shell[j]/tau_binsize))
            f_an = (1-(2*f_factor[i]-1)*costheta[j])
            dL[k] = dL[k]+dLtheta[i]*f_rebin*f_trap
            dL_resp[k] = dL_resp[k]+dLtheta[i]*resp[i]*f_rebin*f_trap
            dL_an[k] = dL_an[k]+dLtheta[i]*f_an*f_rebin*f_trap
            dL_resp_an[k] = dL_resp_an[k]+(
                dLtheta[i]*resp[i]*f_an*f_rebin*f_trap)
    cumdL = cumtrapz(dL,x=tau,initial=0.)
    L = np.trapz(dL,x=tau)
    cumdL_resp = cumtrapz(dL_resp,x=tau,initial=0.)
    L_resp = np.trapz(dL_resp,x=tau)
    cumdL_an = cumtrapz(dL_an,x=tau,initial=0.)
    L_an = np.trapz(dL_an,x=tau)
    cumdL_resp_an = cumtrapz(dL_resp_an,x=tau,initial=0.)
    L_resp_an = np.trapz(dL_resp_an,x=tau)
    ii = range(len(tau))
    tauL = np.trapz([tau[i]*dL[i] for i in ii],x=tau)
    tauL_resp = np.trapz([tau[i]*dL_resp[i] for i in ii],x=tau)
    tauL_an = np.trapz([tau[i]*dL_an[i] for i in ii],x=tau)
    tauL_resp_an = np.trapz([tau[i]*dL_resp_an[i] for i in ii],x=tau)
    cent = tauL/L
    cent_resp = tauL_resp/L_resp
    cent_an = tauL_an/L_an
    cent_resp_an = tauL_resp_an/L_resp_an
    rfL = np.append(np.diff(cumdL_an),0.)
    rfR = np.append(np.diff(cumdL_resp_an),0.)
    if do_plots == 'yes':
        p_ld([t/lightday for t in tau],[dL,dL_an,dL_resp,dL_resp_an],
             title='dL(tau)',axes=[0.3,tau[-1]/lightday,-0.6*max(dL),1.1*max(dL)],
             vlines=[2*max(rr)/lightday],label=label,label_xy=[0.2,0.25],
             legends=[r'$\epsilon$-weighted, iso.',r'$\epsilon$-weighted, aniso.',
                      r'$\eta$-weighted, iso.',r'$\eta$-weighted, aniso.'],
             xtitle=r'$\tau$ [lightdays]',ytitle=r'd$L(\tau)$ (arbitrary units)',log='no')
        p_ld([t/lightday for t in tau],[dL,dL_an,dL_resp,dL_resp_an],
             title='dL(tau)xlog',xlog='yes',
             axes=[0.3,tau[-1]/lightday,
                    -max(dL),1.1*max(dL)],
             vlines=[2*max(rr)/lightday],label=label,
             legends=['E-iso','E-aniso','R-iso','R-aniso'],
             xtitle=r'$\tau$ [lightdays]',ytitle=r'$dL(\tau)$',log='no')
        p_ld([t/lightday for t in tau],[cumdL,cumdL_an],
             title='L(tau)',label=label,axes=[0,300,1e40,5*max(cumdL)],label_xy=[0.3,0.5],
             legends=['Isotropic emission','Anisotropic emission'],
             xtitle=r'$\tau$ (lightdays)',ytitle=r'Enclosed $\log(L/1$ erg s$^{-1})$ ')
    print('int over tau, L: ',"{0:.2f}".format(m.log10(L)),
          'cent:',"{0:.2f}".format(cent/lightday),'cent_ani:',
          "{0:.2f}".format(cent_an/lightday),
          'rcent:',"{0:.2f}".format(cent_resp/lightday),'rcent_ani:',
          "{0:.2f}".format(cent_resp_an/lightday))
    return [L,cumdL,cent_an,cent_resp_an,rfL,rfR]

# =======================================================================================
# s=0
# =======================================================================================

def do_line(line='lya',nH=10.75,r_in=1.,r_out=140.,cont_wl=-1,
            Ncol='22p5',interp='raw',is_cont='no',interpolation=0.1,file='indef',
            inwdfile='indef',extrapolate_fluxes='yes'):
    """integrates over a line (or continuum) from r_in to r_out [lightdays]."""

    r_in = m.log10(r_in*lightday)
    r_out = m.log10(140.*lightday)
    # load data
    if is_cont == 'yes':
        filepath=file
        inwdpath=inwdfile
        line = 'continuum'
    else:
        filepath = match_line(line,Ncol,interp)
        inwdpath = match_line(line+'_inwd',Ncol,interp)
    print('doing line:',line,'nH:',nH)
    header,lognH_tab,logPhi_tab,EW_tab = read_fort(filepath)
    junk,lognH_tab_inwd,logPhi_tab_inwd,EW_tab_inwd = read_fort(inwdpath)
    # arrange EWs on grid
    lognH,logPhi,EWgrid = arrange_grid(lognH_tab,logPhi_tab,EW_tab)
    junk,junk2,EWgrid_inwd = arrange_grid(lognH_tab_inwd,logPhi_tab_inwd,EW_tab_inwd)
    # get a constant-nH slice
    x = ((np.where(lognH == nH))[0])[0]
    EW_slice = EWgrid[:,x]
    EW_slice_inwd = EWgrid_inwd[:,x]
    # convert equivalent widths to emissivities
    incident_nu_F_nu = get_incident_nu_Fnu(logPhi)
    fluxes = [incident_nu_F_nu[i]*EW_slice[i] for i in range(len(logPhi))]
    fluxes_inwd = [incident_nu_F_nu[i]*EW_slice_inwd[i] for i in range(len(logPhi))]
    # convert ionizing cont Phi to radius
    logr = Phi_to_radius(logPhi)
    logr = np.asarray(logr[::-1])
    fluxes = np.asarray(fluxes[::-1])
    fluxes_inwd = np.asarray(fluxes_inwd[::-1])
    # interpolate from edge of flux plot
    if extrapolate_fluxes == 'yes':
        logfluxes = [m.log10(f+1e-19) for f in fluxes]
        logfluxes_inwd = [m.log10(f+1e-19) for f in fluxes_inwd]
        i = ((np.nonzero(fluxes))[0])[0]
        i_inwd = ((np.nonzero(fluxes_inwd))[0])[0]
        slope = logfluxes[i+1]-logfluxes[i]
        slope_inwd = logfluxes_inwd[i_inwd+1]-logfluxes_inwd[i_inwd]
        for j in range(i):
            logfluxes[j] = logfluxes[i]-(i-j)*slope
        for j in range(i_inwd):
            logfluxes_inwd[j] = logfluxes_inwd[i_inwd]-(i_inwd-j)*slope_inwd
        rawfluxes = list(fluxes)
        rawfluxes_inwd = list(fluxes_inwd)
        fluxes = [10.**f for f in logfluxes]
        fluxes_inwd = [10.**f for f in logfluxes_inwd]
        plt.plot(logr,rawfluxes,'b-',label='raw')
        plt.plot(logr,fluxes,'b--',label='extrapolated')
        plt.plot(logr,rawfluxes_inwd,'r-',label='raw inwd')
        plt.plot(logr,fluxes_inwd,'r--',label='extrap inwd')
        plt.axvline(x=m.log10(lightday),color='black')
        plt.yscale('log')
        plt.legend(loc='best')
        plt.savefig('fluxextrap.eps')
        plt.close()
    print(fluxes)
    inwd_ratio = np.asarray([fluxes_inwd[i]/(fluxes[i]+1e-31) for i in range(len(fluxes))])
    # get responsivities, based on flux vs radius.
    print(min(logr))
    fine_logr = np.arange(min(logr),max(logr),0.001)
    fine_logf = np.interp(fine_logr,logr,[m.log10(f+1e-31) for f in fluxes])
    fine_logf_inwd = np.interp(fine_logr,logr,[m.log10(f+1e-31) for f in fluxes_inwd])
    fine_inwd = np.interp(fine_logr,logr,inwd_ratio)
    # discard geometrical ratios for very small fluxes
    fine_inwd = [fine_inwd[i] if fine_logf_inwd[i] > 1e-8 else 0.5 for i in range(len(fine_inwd))]
    # discard geometrical ratios >1 (numerical problems)
    fine_inwd = [fine_inwd[i] if fine_inwd[i] < 1.0 else 1.0 for i in range(len(fine_inwd))]
    # plot fluxes inwd, total
    p_ld(fine_logr,[fine_logf,fine_logf_inwd],title='fluxes_compare',log='yes',
         xtitle=r'$\log(r/1$ cm$)$',label=typeline(line)+r', $\log(n_\mathrm{H})=$'+str(nH),
         ytitle=r'$\epsilon$',vlines=[r_in,r_out],
         legends=[r'$\epsilon_{\mathrm{tot}}$',r'$\epsilon_{\mathrm{inwd}}$'])
    p_ld(fine_logr,fine_inwd,title='F_factor',log='no',hlines=[0.5,1.],
         xtitle=r'$\log(r/1$ cm$)$',label=typeline(line)+r', $\log(n_\mathrm{H})=$'+str(nH),
         ytitle=r'$F_{\mathrm{an}}$',vlines=[r_in,r_out],axes=[14.5,18.5,0.3,1.1])
    resp = responsivities_smooth(fine_logr,fine_logf,smoothing=100,step=1,
                          do_plots='yes',label=typeline(line)+r', $N_{\mathrm{col}}=$'+
                          Ncol+r', $n_{\mathrm{H}}=$'+str(nH))
    # Clip the radius and emissivity arrays to r_in,r_out
    logrclip,fclip = clip_radii(logr,fluxes,r_in,r_out)
    # linearize radii and interpolate to get even spacing
    linr,linf,logr,logf = interpolate_fluxes(logrclip,fclip,interpolation)
    resp = np.interp(linr,[10.**r for r in fine_logr],resp)
    inwd_ratio = np.interp(linr,[10.**r for r in fine_logr],fine_inwd)
    wt_resp = np.trapz([resp[i]*linf[i] for i in range(len(linf))])/np.trapz(linf)
    wt_inwd = np.trapz([inwd_ratio[i]*linf[i] for i in range(len(linf))])/np.trapz(linf)
    p_ld([m.log10(r) for r in linr],resp,title='logr_resp',
         label=typeline(line)+r', $\log(n_\mathrm{H})=$'+str(nH)+', s=0',
         log='no',axes=[15,m.log10(200*lightday),-10,10],hlines=[0,1],
         vlines=[m.log10(1*lightday),m.log10(140*lightday)],xtitle=r'$\log(r/1$ cm$)$',
         ytitle=r'$\eta$')
    p_ld([r/lightday for r in linr],resp,title='lin_resp',label=typeline(line),
         log='no',axes=[0,140,-10,5],hlines=[0,1])
    p_ld(logr,inwd_ratio,title='F_factor_interp',log='no',hlines=[0.5,1.],
         xtitle='log(r/cm)',label=typeline(line)+' log nH:'+str(nH),
         ytitle=r'Cloud $L_{\mathrm{inwd}}/L_{\mathrm{tot}}$',vlines=[r_in,r_out])
    # get A_c scaling
    scale = Ac_scaling_analytic(r_in,r_out,s=0)
    # do integral over r
    L,cumdL,C,cumdC,cent,respcent,cumdetaL = intline(
        linr,linf,resp=resp,Ac_scaling=scale)
    ascii.write([[r/lightday for r in linr],cumdC],line+'_s0_cumdC.dat',names=['r','cumdC'])
    detaL = np.append(np.asarray([0]),np.diff(cumdetaL))
    p_ld([r/lightday for r in linr],[m.log10(L+1e-31) for L in cumdL],
         title='r_L(r)',label=typeline(line),log='no',xtitle=r'$\log(r/1$cm)',
         ytitle=r'Enclosed $\log(L/1$ erg s$^{-1})$',axes=[0,140,38,44.5])
    p_ld([r/lightday for r in linr],cumdetaL,title='r_cumdetaL(r)',
         label=typeline(line),log='no')
    p_ld([r/lightday for r in linr],detaL,title='r_detaL(r)',
         label=typeline(line),log='no')
    print('integral over r, L:',"{0:.2f}".format(m.log10(L)),' centroid/ld:',
        "{0:.2f}".format(cent/lightday),' resp_cent/ld:',"{0:.2f}".format(
            respcent/lightday))
    # do integral over tau
    L_tau,cumdL_tau,cent_tau_aniso,cent_tau_resp_aniso,rfL,rfR = intline_tau(linr,linf,
            resp=resp,Ac_scaling=scale,f_factor=inwd_ratio,do_plots='yes',
            label=typeline(line)+r', s=0, $n_{\mathrm{H}}=$'+str(nH),theta_bins=0.25*len(linr))
    return [L_tau,cent,respcent,cent_tau_aniso,cent_tau_resp_aniso,rfL,rfR,wt_resp,wt_inwd]

#=========================================
# s=2
#=========================================

def U(lognH,logPhi=20):
    return logPhi-lognH-m.log10(c_cms)

def do_line_s2(line='lya',interp='raw',starting_lognH=10.75,starting_logPhi=20,
               starting_logNcol_float = 22.5,r_in=1.,r_out=140.,interpolation=0.1,
               cont_wl=-1,is_cont='no',do_extrapolate='yes',extrapolate_fluxes='yes'):
    step = 0.25 # probably leave this at logphi, lognH grid stepsize
    r_in = m.log10(r_in*lightday)
    r_out = m.log10(r_out*lightday)
    logU = "{0:.2f}".format(starting_logPhi-starting_lognH-m.log10(c_cms))
    print('doing line:',line,'log U:',logU)
    starting_logNcol = lookup_Ncol(starting_logNcol_float)

    # load stating grid (to set slice, etc)
    if is_cont == 'yes':
        gridfile = match_cont(cont_wl,Ncol=starting_logNcol,interp=interp)
    else:
        gridfile = match_line(line=line,Ncol=starting_logNcol,interp=interp)
    header,lognH,logPhi,raw_EW = read_fort(gridfile)
    lognH,logPhi,EW_grid = arrange_grid(lognH,logPhi,raw_EW,stepsize_x=step,
                                        stepsize_y=step)
    lognH_cut,logPhi_cut,EW_cut,i_cut = diagonal_slice_grid(lognH,logPhi,
        starting_lognH,starting_logPhi,EW_grid,stepsize_x=step,stepsize_y=step)
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
    # overlay the s=2 slice on the EW grid plot (for testing)
    range_logr = Phi_to_radius(logPhi)
    #p_EWcontours_radius_overlay(lognH,range_logr,EW_grid,lognH_cut,logr_descending,
    #                     Ncol_cut_descending,r_in=10.**r_in,r_out=10.**r_out)

    # now load a 'block' of grids for flux and inwd spanning the range of Ncol
    Ncol_list = np.asarray([22,22.25,22.5,22.75,23,23.25,23.5,23.75,24])
    Ncol_str_list = [lookup_Ncol(N) for N in Ncol_list]
    fluxes_block = []
    inwd_block = []
    for i in range(len(Ncol_list)):
        if is_cont == 'yes':
            gridfile_N = match_cont(cont_wl,Ncol=Ncol_str_list[i],interp=interp,
                                    key='nFnu')
            inwdfile_N = match_cont(cont_wl,Ncol=Ncol_str_list[i],interp=interp,
                                    key='InwT')
        else:
            gridfile_N = match_line(line=line,Ncol=Ncol_str_list[i],interp=interp)
            inwdfile_N = match_line(line=line+'_inwd',Ncol=Ncol_str_list[i],
                                    interp=interp)
        header_N,lognH_N,logPhi_N,raw_EW_N = read_fort(gridfile_N)
        header_inwd,j2,j3,raw_inwd_EW_N = read_fort(inwdfile_N)
        j1,j2,EW_grid_N = arrange_grid(lognH_N,logPhi_N,raw_EW_N,
                                                  stepsize_x=step,stepsize_y=step)
        lognH_N,logPhi_N,inwd_EW_grid_N = arrange_grid(lognH_N,logPhi_N,raw_inwd_EW_N,
                                                  stepsize_x=step,stepsize_y=step)
        lognH_cut_N,logPhi_cut_N,EW_cut_N,i_cut_N = diagonal_slice_grid(lognH_N,
                                logPhi_N,starting_lognH,starting_logPhi,EW_grid_N)
        j1,j2,inwd_EW_cut_N,j3 = diagonal_slice_grid(lognH_N,
                                logPhi_N,starting_lognH,starting_logPhi,inwd_EW_grid_N)
        fluxes_N = [incident_nu_F_nu[i]*EW_cut_N[i] for i in range(len(logPhi_cut))]
        fluxes_block.append((fluxes_N[::-1]))
        inwd_N = [incident_nu_F_nu[i]*inwd_EW_cut_N[i] for i in range(len(logPhi_cut))]
        inwd_block.append((inwd_N)[::-1])

    # extrapolate grids to low Ncol
    if do_extrapolate == 'yes':
        extrapolate_list = [21,21.25,21.5,21.75]
        extragrids = []
        extragrids_inwd = []
        nsteps = [4,3,2,1]
        zeropoint = np.log10(fluxes_block[0])
        zeropoint_inwd = np.log10(inwd_block[0])
        diffgrid = np.subtract(zeropoint,np.log10(fluxes_block[1]))
        diffgrid_inwd = np.subtract(zeropoint_inwd,np.log10(inwd_block[1]))
        for i in range(len(extrapolate_list)):
            extragrids.append(np.add(zeropoint,np.multiply(nsteps[i],diffgrid)))
            extragrids_inwd.append(np.add(zeropoint_inwd,np.multiply(nsteps[i],
                                                                diffgrid_inwd)))
        extragrids = np.power(10.,extragrids)
        extragrids_inwd = np.power(10.,extragrids_inwd)
        fluxes_block = np.append(extragrids,fluxes_block,axis=0)
        inwd_block = np.append(extragrids_inwd,inwd_block,axis=0)
        Ncol_list = np.append(extrapolate_list,Ncol_list)
    #inwd_block = np.divide(inwd_block,fluxes_block)
    Ncol_str_list = [str(n) for n in Ncol_list]

    # interpolate between fluxes of adjacent Ncol for each r
    nearest_Ncol_cut = np.asarray([takeClosest(Ncol_list,m.log10(N))
                                   for N in Ncol_cut])
    best_fluxes = []
    interp_fluxes = []
    interp_inwd = []
    for i in range(len(logr)):
        farr = [(fluxes_block[j])[i] for j in range(len(Ncol_list))]
        inwdarr = [(inwd_block[j])[i] for j in range(len(Ncol_list))]
        f = np.interp([m.log10(Ncol_cut[i])],Ncol_list,farr)
        interp_fluxes.append((f)[0])
        inwd = np.interp([m.log10(Ncol_cut[i])],Ncol_list,inwdarr)
        interp_inwd.append((inwd)[0])
        iblock = (np.where(Ncol_list == takeClosest(
            Ncol_list,m.log10(Ncol_cut[i])))[0])[0]
        best_fluxes.append((fluxes_block[iblock])[i])
    interp_diff = [(interp_fluxes[i]-best_fluxes[i])/interp_fluxes[i]
               for i in range(len(best_fluxes))] # for testing.
    # interpolate from edge of flux plot
    if extrapolate_fluxes == 'yes':
        logfluxes = [m.log10(f+1e-19) for f in interp_fluxes]
        logfluxes_inwd = [m.log10(f+1e-19) for f in interp_inwd]
        i = ((np.nonzero(interp_fluxes))[0])[0]
        i_inwd = ((np.nonzero(interp_inwd))[0])[0]
        slope = logfluxes[i+1]-logfluxes[i]
        slope_inwd = logfluxes_inwd[i_inwd+1]-logfluxes_inwd[i_inwd]
        for j in range(i):
            logfluxes[j] = logfluxes[i]-(i-j)*slope
        for j in range(i_inwd):
            logfluxes_inwd[j] = logfluxes_inwd[i_inwd]-(i_inwd-j)*slope_inwd
        rawfluxes = list(interp_fluxes)
        rawfluxes_inwd = list(interp_inwd)
        interp_fluxes = [10.**f for f in logfluxes]
        interp_inwd = [10.**f for f in logfluxes_inwd]
        plt.plot(logr,rawfluxes,'b-')
        plt.plot(logr,interp_fluxes,'b--')
        plt.plot(logr,rawfluxes_inwd,'r-')
        plt.plot(logr,interp_inwd,'r--')
        plt.axvline(x=m.log10(lightday),color='black')
        plt.yscale('log')
        plt.savefig('fluxextrap.eps')
        plt.close()
    inwd_ratio = np.asarray([interp_inwd[i]/(interp_fluxes[i]+1e-19) for i in range(len(interp_fluxes))])
    p_ld(logr,inwd_ratio,title='interp_geoF',
         label=typeline(line),log='no',ytitle=r'$F_{\mathrm{an}}$',xtitle=r'$\log(r/1$ cm$)$')
    p_ld([10.**r/lightday for r in logr],interp_fluxes,title='interp_flux',
         label=typeline(line),log='yes')
    # get responsivities, based on flux vs radius.
    fine_logr = np.arange(min(logr),max(logr),0.0001)
    fine_logf = np.interp(fine_logr,logr,[m.log10(f+1e-19) for f in interp_fluxes])
    fine_inwd = np.interp(fine_logr,logr,inwd_ratio) # dummy
    # discard geometrical ratios for very small fluxes
    fine_inwd = [fine_inwd[i] if fine_logf[i] > 1e-8 else 0.5 for i in range(len(fine_inwd))]
    # discard geometrical ratios >1 (numerical problems)
    fine_inwd = [fine_inwd[i] if fine_inwd[i] < 1.0 else 1.0 for i in range(len(fine_inwd))]
    p_ld(fine_logr,fine_inwd,title='F_factor',log='no',hlines=[0.5,1.],
         xtitle=r'$\log(r/1$ cm$)$',label=typeline(line)+r', s=2, $\log(U)=$'+logU,
         ytitle=r'$F_{\mathrm{an}}$',vlines=[r_in,r_out],axes=[14.5,18.5,0.3,1.1])
    resp = responsivities_smooth(fine_logr,fine_logf,smoothing=1000,
                          do_plots='yes',label=typeline(line))
    # Clip the radius and emissivity arrays to r_in,r_out
    logrclip,fclip = clip_radii(logr,interp_fluxes,r_in,r_out)
    # linearize radii and interpolate to get even spacing
    linr,linf,logr,logf = interpolate_fluxes(logrclip,fclip,interpolation)
    resp = np.interp(linr,[10.**r for r in fine_logr],resp)
    inwd_ratio = np.interp(linr,[10.**r for r in fine_logr],fine_inwd)
    wt_resp = np.trapz([resp[i]*linf[i] for i in range(len(linf))])/np.trapz(linf)
    wt_inwd = np.trapz([inwd_ratio[i]*linf[i] for i in range(len(linf))])/np.trapz(linf)
    p_ld([m.log10(r) for r in linr],resp,title='logr_resp',
         label=typeline(line)+r', s=2, $\log(U)=$'+logU,
         log='no',axes=[15,m.log10(200*lightday),-3,3],hlines=[0,1],
         vlines=[m.log10(1*lightday),m.log10(140*lightday)],xtitle=r'$\log(r/1$ cm$)$',
         ytitle=r'$\eta$')
    # get A_c scaling
    scale = Ac_scaling_analytic(r_in,r_out,s=2)
    # do integral
    L,cumdL,C,cumdC,cent,respcent,cumdetaL = intline(linr,linf,resp=resp,
                                        Ac_scaling=scale,s=2.)
    ascii.write([[r/lightday for r in linr],cumdC],line+'_s2_cumdC.dat',names=['r','cumdC'])
    dL = [L/interpolation for L in np.insert(np.diff(cumdL),0,0.)]
    detaL = [L/interpolation for L in np.insert(np.diff(cumdetaL),0,0.)]
    dC = [cov/interpolation for cov in np.insert(np.diff(cumdL),0,0.)]
    print('logLtot:',"{0:.2f}".format(m.log10(L)),' centroid/ld:',
        "{0:.2f}".format(cent/lightday),' resp_cent/ld:',"{0:.2f}".format(
        respcent/lightday))
    L_tau,cumdL_tau,cent_tau_aniso,cent_tau_resp_aniso,rfL,rfR = intline_tau(
        linr,linf,s=2.,label=typeline(line)+r', s=2, $\log(U)=$'+logU,
        resp=resp,Ac_scaling=scale,f_factor=inwd_ratio,
            do_plots='yes')
    p_ld([m.log10(r) for r in linr],[m.log10(L+1e-19) for L in cumdL],
         title='r_L(r)',label=typeline(line),log='no',xtitle=r'$\log(r/1$cm)',
         ytitle=r'Enclosed $\log(L/1$ erg s$^{-1})$',axes=[0,140,38,44.5])
    p_ld([r/lightday for r in linr],dL,title='r_dL(r)',label=typeline(line))
    p_ld([r/lightday for r in linr],detaL,title='r_detaL(r)',
         label=typeline(line),log='no')
    return [L_tau,cent,respcent,cent_tau_aniso,cent_tau_resp_aniso,rfL,rfR,wt_resp,wt_inwd]

# =======================================================================================
# drive response function with continua
# =======================================================================================

def cont_5548(timestep=1):
    """ loads HST continuum lightcurve for NGC 5548. always interpolates to get normaliztion (time avg);
    if interp='no', the original data points are output after flux normalization."""
    filename = 'HST_1157.5.txt'
    # load data
    mjd = []
    flux = []
    e_flux = []
    with open(filename) as f:
        for line in f:
            mjd.append(float((line.split())[0]))
            flux.append(float((line.split())[1]))
            e_flux.append(float((line.split())[2]))
    # interpolate data
    date = [mjd[i]-mjd[0] for i in range(len(mjd))]
    date_interp = np.arange(0,max(date),timestep)
    flux_interp = np.interp(date_interp,date,flux)
    e_flux_interp = np.interp(date_interp,date,e_flux)
    # normalize flux
    mean_flux = flux[0] # np.mean(flux_interp)
    flux_norm = [f/mean_flux for f in flux]
    e_flux_norm = [f/mean_flux for f in e_flux]
    flux_norm_interp = [f/mean_flux for f in flux_interp]
    e_flux_norm_interp = [f/mean_flux for f in e_flux_interp]
    return [date_interp,flux_norm_interp,e_flux_norm_interp,date,flux_norm,e_flux_norm]

def RandomWalk(N=100,amplitude=1):
    """
    Use numpy.cumsum and numpy.random.uniform to generate
    a 2D random walk of length N, each of which has a random DeltaX and
    DeltaY between -1/2 and 1/2.  You'll want to generate an array of 
    shape (N,d), using (for example), random.uniform(min, max, shape).
    """
    return np.cumsum(np.random.normal(0.0,amplitude,N))

def DampedRandomWalk(N=100,tau_char=40,sigma_char=0.04,is_log='yes'):
    """ for random walks, phi = exp(1-/tau), where tau is the characteristic timescale.
    amplitude^2 = tau*sigma^2(1-exp(-2/tau))/2 """
    phi = m.exp(-1./tau_char)
    amplitude = m.sqrt(tau_char*(sigma_char**2.)*(1-m.exp(-2./tau_char))/2.)
    walk = np.zeros(N)
    for i in range(walk.size-1):
        walk[i+1] = phi*walk[i]+np.random.normal(0.0,amplitude)
    p_ld(range(len(walk)),walk,title='drwalk',log='no',xtitle='Days since t=0',ytitle='rel. flux',
         label_xy=[0.2,0.2],
         label=r'DRW, $\tau_{\mathrm{char}}=$'+str(tau_char)+r' $\sigma_{\mathrm{char}}=$'+str(sigma_char))
    p_ld(range(len(walk)),[10.**f for f in walk],title='exp_drwalk',log='no',label_xy=[0.2,0.2],
         xtitle='Days since t=0',ytitle='rel. flux',
        label=r'exp(DRW), $\tau_{\mathrm{char}}=$'+str(tau_char)+r' $\sigma_{\mathrm{char}}=$'+str(sigma_char))
    if is_log == 'yes':
        return [10.**f for f in walk]
    return walk

def make_continuum(time=1000,amplitude=0.01,spike='no'):
    steady = np.zeros(time)
    rand = RandomWalk(time,amplitude)
    if spike == 'yes':
        steady[100] = 1
        steady[200] = 2
        return steady
    else:
        return rand
    
def trimheader(file,nhead=1):
    with open(file, 'r') as fin:
        data = fin.read().splitlines(True)
    with open(file, 'w') as fout:
        fout.writelines(data[nhead:])

def run_ccf(contfile,linefile,datfolder='gen_lightcurves',sampling='1',plot='yes'):
    """ calls White+Peterson1994 CCF code from shell, prints results and plots CCF"""
    pathnow = os.getcwd()
    os.chdir(pathnow+'/'+datfolder)
    child = pexpect.spawnu('ccf')
    child.logfile = sys.stdout
    child.expect('Enter number of files to be read:')
    child.sendline('2')
    child.expect('(?i)cross-correlation:')
    child.sendline(contfile)
    child.expect('(?i)from data file:')
    child.sendline('1')
    child.expect('(?i)correlation:')
    child.sendline(linefile)
    child.expect('(?i)from data file:')
    child.sendline('2')
    child.expect('(?i)divisor(?i)')
    child.sendline(sampling)
    child.expect('(?i)(Y/N)?:')
    child.sendline('Y')
    child.expect('(?i)shift:')
    child.sendline('200')
    child.expect('(?i)for centroid calculation:')
    child.sendline('0.8')
    child.expect('(?i)option:')
    child.sendline('0')
    child.timeout=300
    child.expect(pexpect.EOF)
    #sys.stdout.write (child.after)
    ccfdat = pathnow+'/'+datfolder+'/ccf.dat'
    tau = []
    ccf = []
    if plot == 'yes':
        with open(ccfdat, 'r') as fin:
            data = fin.read().splitlines(True)
            for line in data:
                t,c = [float(x) for x in line.split()]
                tau.append(t)
                ccf.append(c)
        os.chdir(pathnow)
        p_ld(tau,ccf,title='ccf',log='no',xtitle=r'$\tau$ (days)',ytitle='cross-correlation',label='CCF')
    

def drive_line(line='lya',nH=10.75,Ncol='23',time=1000,amplitude=0.1,obscont = 'yes',spike='no',
               damped='yes',tau_char=40,sigma_char=0.04,s=0):
    from scipy.ndimage.filters import convolve1d
    if obscont == 'yes':
        timearray,lc_cont,err_lc_cont,time_orig,lc_cont_orig,err_lc_cont_orig = cont_5548()
        delta_cont = [f-1 for f in lc_cont]
        time = len(timearray)
        contlabel = 'Observed continuum'
    else:
        # make a continuum
        timearray = np.arange(time)
        time_orig = timearray
        if damped == 'yes':
            lc_cont = DampedRandomWalk(N=time,tau_char=tau_char,sigma_char=sigma_char,is_log='yes')
            delta_cont = [f-1 for f in lc_cont]
        else:
            delta_cont = make_continuum(time=time,amplitude=amplitude,spike=spike)
            lc_cont = [1.+dc for dc in delta_cont]
            contlabel = 'RWalk continuum'
        lc_cont_orig = lc_cont
        err_lc_cont = [0.01*f for f in lc_cont]
        err_lc_cont_orig = err_lc_cont
    # retrieve line response
    if s==2:
        L,cent,respcent,cent_tau_aniso,cent_tau_resp_aniso,rfL,rfR,wt_resp,wt_inwd = do_line_s2(
            line=line,starting_lognH=nH,r_in=1.,r_out=140.,cont_wl=1215.,
            interp='raw',interpolation=0.1)
    else:
        L,cent,respcent,cent_tau_aniso,cent_tau_resp_aniso,rfL,rfR,wt_resp,wt_inwd = do_line(
            line=line,nH=nH,r_in=1.,r_out=140.,cont_wl=1215.,
            Ncol=Ncol,interp='raw',interpolation=0.1)
    # convolve response with norm. continuum
    deltaL = np.convolve(rfL,delta_cont)
    deltaL = deltaL[0:time]
    print(len(timearray),len(deltaL))
    lcL = [L+dL for dL in deltaL]
    lcL_norm = [Lt/L for Lt in lcL]
    deltaR = np.convolve(rfR,delta_cont)
    deltaR = deltaR[0:time]
    lcR = [L+dR for dR in deltaR]
    lcR_norm = [Rt/L for Rt in lcR]
    # cut the convolved response to match continuum len
    p_ld(range(len(rfL)),[rfL,rfR],log='no',label=typeline(line),
         legends=['deltaL','deltaR'],axes=[0,75,-0.5*max(rfL),1.1*max(rfL)],
         xtitle='Time [lightdays]',ytitle='Response',title='response_'+line+'_s'+str(s))
    p_ld(timearray,[lc_cont,lcL_norm,lcR_norm],log='no',label='',
         legends=['DRW Continuum',typeline(line)+r', $\epsilon$-weighted resp.',
                typeline(line)+r', $\eta$-weighted resp.'],
         xtitle='Time [lightdays]',ytitle='Relative luminosity',
         title='drive_'+line+'_s'+str(s))
    # write lightcurves to gen-lightcurves; trim headers for 
    datfolder = 'gen_lightcurves'
    ascii.write([time_orig,lc_cont_orig,err_lc_cont_orig],datfolder+'/'+line+'_obs_cont.dat')
    trimheader(datfolder+'/'+line+'_obs_cont.dat')
    ascii.write([timearray,lc_cont,err_lc_cont],datfolder+'/'+line+'_interp_cont.dat')
    trimheader(datfolder+'/'+line+'_interp_cont.dat')
    ascii.write([timearray,lcL_norm,[0.05*L for L in lcL_norm]],datfolder+'/'+line+'_driven_line_L.dat')
    trimheader(datfolder+'/'+line+'_driven_line_L.dat')
    ascii.write([timearray,lcR_norm,[0.05*L for L in lcR_norm]],datfolder+'/'+line+'_driven_line_R.dat')
    trimheader(datfolder+'/'+line+'_driven_line_R.dat')
    # run white+peterson lag code
    run_ccf(line+'_interp_cont.dat',line+'_driven_line_R.dat')

from io import StringIO
class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout

def drive_line_multi(line='lya',nH=10.75,Ncol='22p5',time=1000,tau_char=40,sigma_char=0.04,
                     s=0,N=10,rf=[-1],L_line=-1,is_cont='no',cont_wl=None):
    from scipy.ndimage.filters import convolve1d
    timearray = np.arange(time)
    print('driving line!')
    # retrieve line response
    if rf[0] == -1:
        print('obtaining response function!')
        if s==0:
            L,cent,respcent,cent_tau_aniso,cent_tau_resp_aniso,rfL,rfR,wt_resp,wt_inwd = do_line(is_cont=is_cont,
                line=line,nH=nH,r_in=1.,r_out=140.,Ncol=Ncol,interp='raw',interpolation=0.1,cont_wl=cont_wl)
        elif s==2:
            print('\nnote: startingNcol not implemented yet, defaulting to 22.5\n')
            L,cent,respcent,cent_tau_aniso,cent_tau_resp_aniso,rfL,rfR,wt_resp,wt_inwd = do_line_s2(is_cont=is_cont,
                line=line,starting_lognH=nH,r_in=1.,r_out=140.,starting_logNcol_float=22.5,interp='raw',
                interpolation=0.1,cont_wl=cont_wl)
        else:
            'needs to be s=0 or s=2 !'
    else: # allow piped-in response function
        print('using supplied response function and L')
        rfR = rf
        L = L_line
    # convolve response with norm. continuum
    lags = []
    gperrs = []
    centroids = []
    for i in range(N):
        # generate continuum lc realization
        lc_cont = DampedRandomWalk(N=time,tau_char=tau_char,sigma_char=sigma_char,is_log='yes')
        delta_cont = [f-1 for f in lc_cont]
        err_cont = [0.01*f for f in lc_cont]
        datfolder = 'gen_lightcurves'
        ascii.write([timearray,lc_cont,err_cont],datfolder+'/'+line+'_interp_cont.dat')
        trimheader(datfolder+'/'+line+'_interp_cont.dat')
        # generate line lc realization
        deltaR = np.convolve(rfR,delta_cont)
        deltaR = deltaR[0:time]
        lcR = [L+dR for dR in deltaR]
        lcR_norm = [Rt/L for Rt in lcR]
        ascii.write([timearray,lcR_norm,[0.05*L for L in lcR_norm]],datfolder+'/'+line+'_driven_line_R.dat')
        trimheader(datfolder+'/'+line+'_driven_line_R.dat')
        p_ld(timearray,[lc_cont,lcR_norm],log='no',label='',
            legends=['DRW Continuum',typeline(line)+r', $\eta$-weighted resp.'],
            xtitle='Time [lightdays]',ytitle='Relative luminosity',title='drive_'+line+'_s'+str(s))
        # run white+peterson lag code
        with Capturing() as output:  # capture stdout for manipulation
            run_ccf(line+'_interp_cont.dat',line+'_driven_line_R.dat')
        print(output)
        try:
            lags.append(float(((output[-7]).split())[-1]))
            gperrs.append(float(((output[-2]).split())[-1]))
            centroids.append(float(((output[-1]).split())[-1]))
            print(i,lags[i],gperrs[i],centroids[i])
        except ValueError:
            print('ICCF method cant centroid for iteration ',i)
            lags.append(np.mean(lags))
            gperrs.append(np.mean(gperrs))
            centroids.append(np.mean(centroids))
    print('mean lag:',np.mean(lags),' mean gperr:',np.mean(gperrs),' mean centroid:',np.mean(centroids),
          ' stddev lag:',np.std(lags),' stddev centroid:',np.std(centroids))
    return [np.mean(lags),np.std(lags),np.mean(centroids),np.std(centroids),np.mean(gperrs)]

def drive_line_loop(line='lya',nH=10.75,Ncol='22p5',time=1000,sigma_char=0.04,s=0,N=100,
                    Tmin=10,Tmax=250,deltaT=10):
    mean_lags = []
    std_lags = []
    mean_centroids = []
    std_centroids = []
    mean_errs = []
    T_char_arr = range(Tmin,Tmax,deltaT)
    for T_char in T_char_arr:
        mean_lag,std_lag,mean_centroid,std_centroid,mean_err = drive_line_multi(
            line=line,nH=nH,Ncol=Ncol,time=time,tau_char=T_char,sigma_char=0.04,s=s,N=N)
        mean_lags.append(mean_lag)
        std_lags.append(std_lag)
        mean_centroids.append(mean_centroid)
        std_centroids.append(std_centroid)
        mean_errs.append(mean_err)
        for i in range(len(mean_lags)):
            print(T_char_arr[i],mean_lags[i],std_lags[i],mean_centroids[i],std_centroids[i],mean_errs[i])
    ascii.write([T_char_arr,mean_lags,std_lags,mean_centroids,std_centroids,mean_errs],
                line+'_T_char_lag.dat',names=['Tchar','lag','e_lag','cent','e_cent','e_GP'])
    p_ld(T_char_arr,[mean_lags,mean_centroids],legends=['CCF Peak','CCF Centroid'],
         xtitle=r'Continuum $T_{\mathrm{char}}$ [days]',ytitle='Lag [days]',log='no',
         title=line+'_T_char',label=typeline(line)+r' $n_{\mathrm{H}}=$'+str(nH))
    
# =======================================================================================
# diffuse continua
# =======================================================================================

def incident_continuum(wl,beta=-1.5,ref_wl=1215.,flux='no'):
    logPhi = [18.]
    nu_F_nu_ref = (get_incident_nu_Fnu(logPhi))[0]
    logr = (Phi_to_radius(logPhi))[0]
    nu_F_nu = nu_F_nu_ref*(wl**beta/ref_wl**beta)
    nu_L_nu = 4*m.pi*nu_F_nu*((10.**logr)**2.)
    if flux == 'yes':
        return nu_F_nu
    else:
        return nu_L_nu

def incident_continuum_mehdipour(wl,alt='no'):
    from bisect import bisect_left
    filename = 'continuum_mehdipour.dat'

    cont_radius = 14.8*lightday
    scaling_lumin = 10.**43.689455 # luminosity at scaling_wl
    scaling_wl = 1367. # ang
    energy_rydberg = []
    norm_nuFnu = []
    with open(filename) as f:
        for line in f:
            if "#Cont" in line:
                header = line
            else:
                energy_rydberg.append(float((line.split())[0]))
                norm_nuFnu.append(float((line.split())[1]))
    ryd_eV = 13.605693
    energy_eV = [ryd_eV*nu for nu in energy_rydberg]
    eV = 1.602176565e-19
    energy_J = [energy*eV for energy in energy_eV]
    h = 6.63e-34 # in m2 kg s-1
    c_ms = 299792458
    nu_Hz = [energy/h for energy in energy_J]
    wl_m = [c_ms/nu for nu in nu_Hz]
    wl_ang = [wl*1e10 for wl in wl_m]
    wl_ang_ordered = wl_ang[::-1]
    norm_nuFnu = norm_nuFnu[::-1]
    scaling_wlarr = np.arange(1000.,5000.,1.)
    scaling_interp = np.interp(scaling_wlarr,wl_ang_ordered,norm_nuFnu)
    i = bisect_left(scaling_wlarr,scaling_wl)
    scaling_factor = scaling_lumin/scaling_interp[i]
    interp_nuFnu_output = np.interp(wl,wl_ang_ordered,norm_nuFnu)
    nuLnu = [energy*scaling_factor for energy in interp_nuFnu_output]
    nuLnu_alt = [4*m.pi*(cont_radius**2.)*energy for energy in interp_nuFnu_output]
    j = bisect_left(wl,scaling_wl)
    p_wl(wl_ang_ordered,[F*scaling_factor for F in norm_nuFnu],ylabel=r'Continuum $\lambda L_\lambda$ [erg s$^{-1}$]',
         title='mp_cont_all',log='yes',legends='none',label='',xlog='yes',range=[0.001,1e8,5e25,3e44])
    p_wl(wl,interp_nuFnu_output,ylabel='cont',title='mp_cont_nuFnu',log='yes',
         legends='none',label='')
    p_wl(wl,nuLnu,ylabel='cont',title='mp_cont_interp',log='yes',
         legends='none',label='')
    p_wl(wl,nuLnu_alt,ylabel='cont',title='mp_cont_interp_alt',log='yes',
         legends='none',label='')
    print('nuLnu at scaling wl:',nuLnu_alt[j])
    if alt == 'yes':
        return nuLnu_alt
    else:
        return nuLnu
    
def do_continuum(nH=10.75,start_Ncol=22.5,interp='raw',minwl=1000.,maxwl=7500.,s=0.,N_drive=0,T_char = 40):
    
    contwl,contfiles = diff_cont_files(Ncol=lookup_Ncol(start_Ncol),interp=interp,
                                       minwl=minwl,maxwl=maxwl)
    junk,inwdfiles = diff_cont_files(Ncol=lookup_Ncol(start_Ncol),interp=interp,minwl=minwl,
                                     maxwl=maxwl,key='InwT')
    mp_cont = incident_continuum_mehdipour(contwl)
    diffcont = []
    cent = []
    respcent = []
    cent_tau_aniso = []
    cent_tau_resp_aniso = []
    weighted_resp = []
    weighted_inwd = []
    mean_lags = []
    std_lags = []
    mean_centroids = []
    std_centroids = []
    ii = range(len(contwl))
    folder = 'cont_plots/'+lookup_Ncol(start_Ncol)+'/'+str(nH)
    if not os.path.exists(folder):
        os.makedirs(folder)
    for i in ii:
        print(contwl[i],'/',max(contwl),' percent done: ',100.*i/float(ii[-1]))
        if s == 2:
            L,C,rC,ctau,cresptau,rfL,rfR,wt_resp,wt_inwd =do_line_s2(line='continuum',
                            starting_lognH=nH,starting_logPhi=20,
               starting_logNcol_float = start_Ncol,cont_wl=contwl[i],is_cont='yes')
        else:
            L,C,rC,ctau,cresptau,rfL,rfR,wt_resp,wt_inwd = do_line(file=contfiles[i],nH=nH,
                                r_in=1.,r_out=140.,cont_wl=1215.,is_cont='yes',line='continuum',
                                inwdfile=inwdfiles[i])
        if N_drive > 0:
            time = 500
            mean_lag,std_lag,mean_centroid,std_centroid,mean_err = drive_line_multi(L_line=L,rf=rfR,cont_wl=contwl[i],
            line='continuum',is_cont='yes',nH=nH,Ncol=lookup_Ncol(start_Ncol),time=time,tau_char=T_char,sigma_char=0.04,s=s,N=N_drive)
            mean_lags.append(mean_lag)
            std_lags.append(std_lag)
            mean_centroids.append(mean_centroid)
            std_centroids.append(std_centroid)
        diffcont.append(L)
        cent.append(C/lightday)
        respcent.append(rC/lightday)
        cent_tau_aniso.append(ctau/lightday)
        cent_tau_resp_aniso.append(cresptau/lightday)
        weighted_resp.append(wt_resp)
        weighted_inwd.append(wt_inwd)

    relcontmp = [diffcont[i]/(diffcont[i]+mp_cont[i]) for i in ii]
    ratiompcont = [diffcont[i]/mp_cont[i] for i in ii]
    dil_cent = [relcontmp[i]*cent[i] for i in ii]
    dil_respcent = [relcontmp[i]*respcent[i] for i in ii]
    dil_cent_tau_aniso = [relcontmp[i]*cent_tau_aniso[i] for i in ii]
    dil_cent_tau_resp_aniso = [relcontmp[i]*cent_tau_resp_aniso[i] for i in ii]
    p_wl(contwl,diffcont,title='diff_cont_lumin_nH'+str(nH)+'s'+str(s),log='yes',
         ylabel=r'$L_{\mathrm{C,diff}}$ [ergs s$^{-1}$]')
    p_wl(contwl,mp_cont,title='mp_cont_lumin_nH'+str(nH)+'s'+str(s),log='yes',
         ylabel=r'$L_{\mathrm{C,nuc}}$ [ergs s$^{-1}$]')
    p_wl(contwl,ratiompcont,title='mp_cont_ratio_nH'+str(nH)+'s'+str(s),
         ylabel=r'$\nu L_{\nu}(\mathrm{diff})/\nu L_{\nu}(\mathrm{nuc})$',log='no')
    p_wl_multi(contwl,[cent,respcent,cent_tau_aniso,cent_tau_resp_aniso],
               title='diff_cont_R_nH'+str(nH)+'s'+str(s),log='no',
               ylabel=r'Delay centroid [lightdays]',
               legends=[r'$R_{L}$',r'$R_{resp}$',r'Aniso.$R_{L}$',r'Aniso.$R_{resp}$'])
    p_wl_multi(contwl,
                [dil_cent,dil_respcent,dil_cent_tau_aniso,dil_cent_tau_resp_aniso],
                title='diluted_cont_R_nH'+str(nH)+'s'+str(s),log='no',
                ylabel=r'Delay centroid [lightdays]',
                legends=[r'$R_{L}$',r'$R_{resp}$',r'Aniso.$R_{L}$',r'Aniso.$R_{resp}$'])
    if s == 2:
        contdata = {'s':s,'starting_lognH':nH,'starting_logPhi':20,'starting_logNcol':start_Ncol,
                    'contwl':contwl,'diffcont':diffcont,'cent':cent,'respcent':respcent,
                    'ionizing_cont':mp_cont,'cont_ratio':ratiompcont,'weighted_resp':weighted_resp,
                    'weighted_inwd':weighted_inwd,'cent_tau_aniso':cent_tau_aniso,
                    'cent_tau_resp_aniso':cent_tau_resp_aniso,'mean_lags':mean_lags,
                    'std_lags':std_lags,'mean_centroids':mean_centroids,'std_centroids':std_centroids}
        outfilename = 'continuum_s2_startnH'+str(nH)+'_startNcol'+str(start_Ncol)+'_Tchar'+str(T_char)+'.js'
    else:
        contdata = {'s':0,'lognH':nH,'logPhi':20,'logNcol':start_Ncol,
                    'contwl':contwl,'diffcont':diffcont,'cent':cent,'respcent':respcent,
                    'ionizing_cont':mp_cont,'cont_ratio':ratiompcont,'weighted_resp':weighted_resp,
                    'weighted_inwd':weighted_inwd,'cent_tau_aniso':cent_tau_aniso,
                    'cent_tau_resp_aniso':cent_tau_resp_aniso,'mean_lags':mean_lags,
                    'std_lags':std_lags,'mean_centroids':mean_centroids,'std_centroids':std_centroids}
        outfilename = 'continuum_s0_nH'+str(nH)+'_Ncol'+str(start_Ncol)+'_Tchar'+str(T_char)+'.js'
    with open(outfilename, 'w') as outfile:
        json.dump(contdata,outfile)
    return [contwl,diffcont,cent,respcent,mp_cont,ratiompcont]
