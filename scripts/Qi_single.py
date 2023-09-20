import os 
import glob
import shutil
import obspy
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

import scipy
import math

from Esyn_func import *
from plotting import *
from obspy.signal.filter import bandpass


path="/home/kffeng/DVV_kura/COR/UU06-22/RCF_H5/"
sfiles = sorted(glob.glob(os.path.join(path, 'SPU*.h5')))
sta_pair=sfiles[0].split("H5/")[1].split(".h5")[0]
print("# Read in file: ",sfiles,"\n # Station Pair: ",sta_pair)

# component
comp_arr = ["ZN", "ZE","NE", "EN","EZ","NZ"]
num_cmp=len(comp_arr)
fnum=len(sfiles)

lag=100
samp=20
leng=int(lag*samp*2+1)
npts=leng
print("\n# Lag-time (sec): ",lag,"\n# Sampling rate: ",samp,"\n# Data length (npts): ",leng)
# targeted frequency band for waveform monitoring
freq = [0.5, 1, 2, 4]  
# smoothing window lengths corresponding to the frequency bands
winlen=[8,4,2]   
nfreq = len(freq) - 1
print("# Target Frequency Bands: ", freq)
print("# Smoothing window length: ", winlen)


stackf=np.ndarray((fnum,num_cmp,2,leng))
print("Array shape: ",stackf.shape)
vdist=np.zeros((fnum,1))  # S-R distance array
fname=[]                  # file name array
aa=0
# loop through each station-pair
for sfile in sfiles:
    ncmp=0
    fname.append(sta_pair)
    for ccomp in comp_arr:
        print(aa, sfile, ccomp)
        # read stacked waveforms
        if ( read_pyasdf(sfile,ccomp) == None):
            continue
        dist,dt, tvec,sdata = read_pyasdf(sfile,ccomp) # read waveform from pyasdf
        stackf[aa][ncmp]=[tvec,sdata]
        vdist[aa]=dist
        ncmp=ncmp+1
    # -- plot
    plot_waveforms(num_cmp,stackf[aa],fname[aa],comp_arr)

    aa=aa+1
fnum=len(fname)
print("# total files number: ",fnum)


# --- filtering and mean-squared values
indx = npts // 2                      # half-side number of points
MSE=np.ndarray((fnum,num_cmp,nfreq+1,npts)) # filtered two-side averaged stack CF

for aa in range (fnum):
    dafbp=np.ndarray((nfreq,npts))

    for ncmp in range (len(comp_arr)):
        ccomp=comp_arr[ncmp]
        print(fname[aa],ccomp)

        for fb in range(nfreq):
            fmin=freq[fb]
            fmax=freq[fb+1]
            tt = np.arange(0, npts) * dt
            data = stackf[aa][ncmp][1]
            dafbp[fb] = bandpass(data, fmin, fmax, int(1 / dt), corners=4, zerophase=True)

        MSE[aa][ncmp]=[stackf[aa][ncmp][0],dafbp[0],dafbp[1],dafbp[2]]
        # -- plot
        #plot_filtered_waveforms(freq,stackf[aa][ncmp][0],dafbp,fname[aa],ccomp)

# --- calculate average msv
msv=np.ndarray((fnum,num_cmp,nfreq+1,npts))
msv_mean=np.ndarray((fnum,nfreq+1,npts))
msv[:][:][:][:]=0.
msv_mean[:][:][:]=0.

data_sym=np.ndarray((nfreq,indx+1)) # two-side averaged stack CF
fmsv_mean=np.ndarray((fnum,nfreq+1,indx+1))

level=np.ndarray((fnum,nfreq,1))
ratio=3
twinbe=np.ndarray((fnum,nfreq,2))

for aa in range(fnum):
    for ncmp in  range(len(comp_arr)):
        ccomp=comp_arr[ncmp]
        msv[aa][ncmp][0]=MSE[aa][ncmp][0][:]
        for fb in range(nfreq):
            data=MSE[aa][ncmp][fb+1][:]
            fmin=freq[fb]
            fmax=freq[fb+1]

            para = { 'winlen':winlen[fb], 'dt':dt , 'npts': len(data)}
            msv[aa][ncmp][fb+1]=get_smooth(data, para)
            
            msv[aa][ncmp][fb+1]=msv[aa][ncmp][fb+1]/np.max(msv[aa][ncmp][fb+1])  # self-normalized 
    
            # self-normalized
            msv[aa][ncmp][fb+1]=msv[aa][ncmp][fb+1]/np.max(msv[aa][ncmp][fb+1])  

    msv_mean[aa][0]=msv[aa][0][0][:]
    for fb in range(nfreq):
        fmin=freq[fb]
        fmax=freq[fb+1]
        for ncmp in range(len(comp_arr)):
            msv_mean[aa][fb+1]+=msv[aa][ncmp][fb+1][:]
        msv_mean[aa][fb+1]=msv_mean[aa][fb+1]/len(comp_arr)
        # --- fold
        sym=get_symmetric(msv_mean[aa][fb+1],indx)
        data_sym[fb]=sym
        Val_mad=mad(sym)
        level[aa][fb]=Val_mad*ratio
        
        twinbe[aa][fb][0]=int((1/fmin)*4)
        for pt in range(len(sym)):
            if (sym[pt] < float(level[aa][fb])):
                twinbe[aa][fb][1]=int(msv[aa][0][0][indx+pt])
                print(aa,fb,pt,sym[pt],level[aa][fb],twinbe[aa][fb])
                break
            if ( pt >= int((twinbe[aa][fb][0]+(1/fmin)*20)/dt)):
                twinbe[aa][fb][1]=twinbe[aa][fb][0]+int((1/fmin)*20)
                print("*20 times criteria ",aa,fb,pt,sym[pt],level[aa][fb],twinbe[aa][fb])
                break

    fmsv_mean[aa]=[msv[aa][0][0][indx:],data_sym[0],data_sym[1],data_sym[2]]
    # -- plot
    plot_envelope(comp_arr,freq,msv[aa],msv_mean[aa],fname[aa],vdist[aa])
    plot_fmsv_waveforms(freq,fmsv_mean[aa],fname[aa],level[aa],twinbe[aa])

# --- measuring coda window
cvel=[2.6, 2.0, 1.8]    # Rayleigh wave velocities over the freqency bands
mfpx=np.zeros(1)        # mean_free_path search array
intby=np.zeros(30)      # intrinsic_b search array

SSR_final=np.ndarray((len(mfpx),len(intby)))
SSR=np.ndarray((nfreq,len(mfpx),len(intby)))
SSR[:][:][:]=0.

for fb in range(nfreq):
    fmin=freq[fb]
    fmax=freq[fb+1]
    c=cvel[fb]
    SSR_final[:][:]=0.
    vdist[:]=0.000001  # To avoid zero value at denominator
    
    # parameters for getting the sum of squared residuals (SSR) between Eobs and Esyn 
    para={ 'fb':fb, 'vdist':vdist, 'npts':npts, 'dt':dt, 'cvel':c, \
        'mfp':mfpx, 'intb':intby,'twin':twinbe, 'fmsv':fmsv_mean }
    # call function get_SSR
    SSR_final, mfpx, intby = get_SSR(fnum, para )
    SSR[fb]=SSR_final

# --- Final result 
result_intb=np.ndarray((nfreq))
result_mfp=np.ndarray((nfreq))

Eobs=np.ndarray((nfreq,npts//2+1))
Esyn=np.ndarray((nfreq,npts//2+1))
aa=0
r=np.take(vdist[aa],0) #+0.000001
print("Final: ",sta_pair)
for fb in range(nfreq):
    fmin=freq[fb]
    fmax=freq[fb+1]

    # parameters for getting optimal value from the sum of squared residuals (SSR) between Eobs and Esyn 
    para={ 'fb':fb, 'fmin':fmin, 'fmax':fmax, 'vdist':vdist, 'npts':npts, 'dt':dt, 'cvel':c, 'filenum':aa, \
        'mfp':mfpx, 'intb':intby,'twin':twinbe, 'fmsv':fmsv_mean, 'SSR':SSR , 'sta':sta_pair}
    # call function get_optimal
    result_intb[fb], result_mfp[fb], Eobs[fb], Esyn[fb] = get_optimal(fnum,para)
    
    plot_fitting_result(result_mfp[fb],result_intb[fb],fmsv_mean[aa][0][:],
                        Eobs[fb],Esyn[fb],sta_pair,vdist[aa],twinbe[aa][fb],fmin,fmax)
#plot_grid_searching(sta_pair,freq,SSR,x,y)

# --- move figures to dics
root=os.getcwd()
output_imgdir = root+"/../figures"
if not os.path.exists(output_imgdir):
        os.makedirs(output_imgdir)
images= [f for f in os.listdir() if '.png' in f.lower()]
for image in images:
    shutil.move(image, output_imgdir)
