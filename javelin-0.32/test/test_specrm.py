import numpy as np
import matplotlib.pyplot as plt
from javelin.predict import PredictSignal, PredictRmap, generateLine, generateError, PredictSpear
from javelin.lcio import *
from javelin.zylc import LightCurve, get_data
from javelin.lcmodel import Cont_Model, Rmap_Model, Pmap_Model
from scipy.stats import norm

""" Tests from scratch.
"""

def test(contfile='con.dat',linefile='yelm.dat'):
    # make a damped-RW continuum model
    # - this provides priors on the continuum DRW params for the next step!
    c = get_data(contfile)
    cmod = Cont_Model(c)
    cmod.do_mcmc(threads=8)
    # make a continuum+line model
    reverb_data = get_data([contfile,linefile],names=['continuum','line'])
    rmod = Rmap_Model(reverb_data)
    rmod.do_mcmc(conthpd=cmod.hpd,fchain='output.dat',lagtobaseline=0.3,laglimit=[[0,150]],threads=8)
    rmod.get_hpd()
    cont_sigma,cont_tau,lag,width,scale = np.loadtxt('output.dat', skiprows=1, unpack=True)
    lagcent,lagsigma = norm.fit(lag)
    #with open('output.dat') as myfile:
    #    results = list(myfile)[-1]
    #    results = [float(x) for x in results.split()]
    #    print(results)
    #    print("%.2f"%float(results[2]),"%.2f"%float(results[3]),"%.2f"%float(results[4]))
    #rmod.show_hist(figout='histograms',figext='eps')
    return [lagcent,lagsigma]
if __name__ == "__main__":
    lagcent,lagsigma = test()
    print(lagcent,lagsigma)
