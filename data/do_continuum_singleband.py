# drive continuum and read lags for a single continuum wl band.

from ewgrids import *

def do_continuum_singleband(nH=10,start_Ncol=22.5,interp='raw',wl=3400,s=0.,N_drive=100):
    
    contwl,contfiles = diff_cont_files(Ncol=lookup_Ncol(start_Ncol),interp=interp,
                                       minwl=1000,maxwl=10000)
    junk,inwdfiles = diff_cont_files(Ncol=lookup_Ncol(start_Ncol),interp=interp,minwl=1000,
                                     maxwl=10000,key='InwT')

    i,value = find_nearest(np.asarray(contwl),wl)
    print('wavelength:',contwl[i])
    if s == 2:
        L,C,rC,ctau,cresptau,rfL,rfR,wt_resp,wt_inwd =do_line_s2(line='continuum',
                        starting_lognH=nH,starting_logPhi=20,starting_logNcol_float = start_Ncol,
                        cont_wl=contwl[i],is_cont='yes')
    else:
        L,C,rC,ctau,cresptau,rfL,rfR,wt_resp,wt_inwd = do_line(file=contfiles[i],nH=nH,
                    r_in=1.,r_out=140.,cont_wl=1215.,is_cont='yes',line='continuum',
                    inwdfile=inwdfiles[i])
    if N_drive > 0:
        time = 500
        T_char = 40
        mean_lag,std_lag,mean_centroid,std_centroid,mean_err = drive_line_multi(L_line=L,rf=rfR,cont_wl=contwl[i],
        line='continuum',is_cont='yes',nH=nH,Ncol=lookup_Ncol(start_Ncol),time=time,tau_char=T_char,sigma_char=0.04,s=s,N=N_drive)

    return [contwl[i],L,C,rC,ctau,cresptau,mean_lag,std_lag,mean_centroid,std_centroid,mean_err]

do_continuum_singleband(nH=10.,N_drive=50,wl=3400.0)
