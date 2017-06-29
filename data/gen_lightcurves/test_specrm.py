import numpy as np
import matplotlib.pyplot as plt
from javelin.predict import PredictSignal, PredictRmap, generateLine, generateError, PredictSpear
from javelin.lcio import *
from javelin.zylc import LightCurve, get_data
from javelin.lcmodel import Cont_Model, Rmap_Model, Pmap_Model
from scipy.stats import norm

""" Tests from scratch.
"""

def test(contfile='con.dat',linefile1='yelm.dat',laglimit=200):
    # make a damped-RW continuum model
    # - this provides priors on the continuum DRW params for the next step!
    print('doing lag for '+contfile+' - '+linefile1)
    c = get_data(contfile)
    cmod = Cont_Model(c)
    cmod.do_mcmc(nwalkers=200,nburn=200,nchain=200)
    # make a continuum+lines model
    reverb_data = get_data([contfile,linefile1],names=['continuum','line1'])
    rmod = Rmap_Model(reverb_data)
    rmod.do_mcmc(conthpd=cmod.hpd,fchain='output.dat',laglimit=[[0,laglimit]])
    rmod.get_hpd()
    cont_sigma,cont_tau,lag,width,scale = np.loadtxt('output.dat', skiprows=1, unpack=True)
    lagcent,lagsigma = norm.fit(lag)
    rmod.show_hist(figout='histograms',figext='eps')
    return [lagcent,lagsigma]

def reload_chain(contfile='con.dat',linefile1='yelm.dat',laglimit=50):
    reverb_data = get_data([contfile,linefile1])
    rmod = Rmap_Model(reverb_data)
    rmod.load_chain("output.dat")
    rmod.break_chain([[0, laglimit],])
    cont_sigma,cont_tau,lag,width,scale = np.loadtxt('output.dat', skiprows=1, unpack=True)
    lag = [l for l in lag if l<laglimit]
    lagcent,lagsigma = norm.fit(lag)
    rmod.show_hist()
    print('line1 lag centroid (cut):',lagcent,'sigma:',lagsigma)

def test_2lines(contfile='con.dat',linefile1='yelm.dat',linefile2='zing.dat',laglimit=200):
    # make a damped-RW continuum model
    # - this provides priors on the continuum DRW params for the next step!
    print('doing lag for ',contfile,' - ',linefile1,linefile2)
    c = get_data(contfile)
    cmod = Cont_Model(c)
    cmod.do_mcmc(nwalkers=200,nburn=200,nchain=200)
    # make a continuum+lines model
    reverb_data = get_data([contfile,linefile1,linefile2],names=['continuum','line1','line2'])
    rmod = Rmap_Model(reverb_data)
    rmod.do_mcmc(conthpd=cmod.hpd,fchain='output.dat',laglimit=[[0,laglimit],[0,laglimit]])
    rmod.get_hpd()
    cont_sigma,cont_tau,lag,width,scale,lag2,width2,scale2 = np.loadtxt('output.dat', skiprows=1, unpack=True)
    lagcent,lagsigma = norm.fit(lag)
    rmod.show_hist(figout='histograms_2lines',figext='eps')
    return [lagcent,lagsigma]

if __name__ == "__main__":
    #lagcent,lagsigma = test_2lines(contfile='lyaobs_cont.dat',linefile1='lyadriven_line_L.dat',
    #                               linefile2='civdriven_line_L.dat',laglimit=80)
    lagcent,lagsigma = test(contfile='mgii_obs_cont.dat',linefile1='mgii_driven_line_R.dat',laglimit=100)
    print('Line1 lag centroid:',lagcent,'line1 lag width:',lagsigma)
