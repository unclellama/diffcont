# codes to work with the NGC 5548 photoionization grids supplied by Otho Adam Ulrich

import numpy as np
from glob import glob
import shutil
import os
import math as m
from scipy.integrate import simps
from ewgrids_plot import *

lightday = 2.59e15 # cm

def test_clip_radii(r_in,r_out):
    lr = np.asarray([1.,2.,3.,4.,5.,6.,7.,8.])
    em = np.asarray([1.,5.,6.,8.,2.,2.,1.,0.5])
    lrclip,emclip = clip_radii(lr,em,r_in,r_out)
    print('input logr:',lr)
    print('input emis:',em)
    print('clipped logr:',lrclip)
    print('clipped emis:',emclip)

def test_Ac_scaling(r_in=15.,r_out=19.):
    lr = np.arange(14.,21.,0.5)
    print(lr)
    scaling = Ac_scaling_numerical(lr,s=0,r_in=r_in,r_out=r_out,verbose='no')
    emissivity = np.asarray(lr)
    logrclip,emclip = clip_radii(lr,emissivity,r_in,r_out)
    lin_r = np.asarray([10.0**r for r in logrclip])
    dC = [4*m.pi*scaling*(r**(-3./2.)) for r in lin_r]
    C = np.trapz(dC,x=lin_r)
    print('scaling: ',scaling,' C/4pi: ',C/(4.*m.pi))
    
