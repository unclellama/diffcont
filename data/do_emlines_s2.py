# =======================================================================================
# get BEL luminosities + weighted radii for a grid of nH, Ncol, for s=0 models
# - basically calls do_line() for each line, nH, Ncol combination.
# =======================================================================================

r_in = 1. # in lightdays

from ewgrids import *
import json

start_nHlist = [9.5,9.75,10,10.25,10.5,10.75,11]
start_logPhi = 20
start_logNcol = 22
#linelist = ['lya','civ','halpha','hbeta','he4686','he1640','mgii']

linelist = ['lya']

logU_list = [start_logPhi-nH-m.log10(c_cms) for nH in start_nHlist]

for line in linelist:
    linedata = {'line':line,'linename':typeline(line),'U_arr':logU_list}
    for nH in start_nHlist:
        lineL = []
        lineR = []
        lineR_resp = []
        lineT_aniso = []
        lineT_aniso_resp = []
        logU = start_logPhi-nH-m.log10(c_cms)
        print('start_nH:',nH,'logU:',logU)
        L,cent,respcent,cent_tau_aniso,cent_tau_resp_aniso,rfL,rfR = do_line_s2(
            starting_lognH=nH,starting_logPhi=start_logPhi,line=line,
            starting_logNcol_float=start_logNcol)
        lineL.append(L)
        lineR.append(cent)
        lineR_resp.append(respcent)
        lineT_aniso.append(cent_tau_aniso)
        lineT_aniso_resp.append(cent_tau_resp_aniso)
        Udata = {'logU':logU,'L':lineL,'R':lineR,'R_resp':lineR_resp,
                    'T_aniso':lineT_aniso,'T_aniso_resp':lineT_aniso_resp}
        linedata['data_logU_'+"{0:.2f}".format(logU)] = Udata
    with open('lineresults_s2/s2_linedata_'+line+'.js', 'w') as outfile:
        json.dump(linedata,outfile)
