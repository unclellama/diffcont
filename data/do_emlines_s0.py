# =======================================================================================
# get BEL luminosities + weighted radii for a grid of nH, Ncol, for s=0 models
# - basically calls do_line() for each line, nH, Ncol combination.
# =======================================================================================

r_in = 1. # in lightdays

from ewgrids import *
import json

Ncollist = ['22','23','24']
#nHlist = [8,9]

nHlist = [7,7.25,7.5,7.75,8,8.25,8.5,8.75,
          9,9.25,9.5,9.75,10,10.25,10.5,10.75,11,11.25,11.5,11.75,
          12,12.25,12.5,12.75,13,13.25,13.5,13.75,14]

linelist = ['lya','civ','halpha','hbeta','he4686','he1640','mgii']

for line in linelist:
    linedata = {'line':line,'linename':typeline(line),'Ncol_arr':Ncollist}
    for Ncol in Ncollist:
        lineL = []
        lineR = []
        lineR_resp = []
        lineT_aniso = []
        lineT_aniso_resp = []
        for nH in nHlist:
            print('\nNcol:',Ncol,'nH:',nH)
            L,cent,respcent,cent_tau_aniso,cent_tau_resp_aniso,rfL,rfR,wt_resp,wt_inwd = do_line(
                Ncol=Ncol,line=line,nH=nH,r_in=r_in)
            lineL.append(L)
            lineR.append(cent)
            lineR_resp.append(respcent)
            lineT_aniso.append(cent_tau_aniso)
            lineT_aniso_resp.append(cent_tau_resp_aniso)
        Ncoldata = {'Ncol':Ncol,'nHarr':nHlist,'L':lineL,'R':lineR,'R_resp':lineR_resp,
                    'T_aniso':lineT_aniso,'T_aniso_resp':lineT_aniso_resp}
        linedata['data_Ncol_'+Ncol] = Ncoldata
    with open('lineresults_s0/s0_linedata_'+line+'.js', 'w') as outfile:
        json.dump(linedata,outfile)
