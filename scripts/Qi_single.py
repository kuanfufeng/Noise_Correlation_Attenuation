import os 
import glob
import shutil
import obspy
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from Esyn_utils import *
from plotting import *
import scipy
import math

from obspy.signal.filter import bandpass


path="/home/kffeng/DVV_kura/COR/UU06-22/RCF_H5/"
sfiles = sorted(glob.glob(os.path.join(path, 'SPU*.h5')))
sta_pair=sfiles[0].split("H5/")[1].split(".h5")[0]
print("# Read in file: ",sfiles,"\n # Station Pair: ",sta_pair)

# compoent
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

for aa in range(fnum):
    for ncmp in  range(len(comp_arr)):
        ccomp=comp_arr[ncmp]
        msv[aa][ncmp][0]=MSE[aa][ncmp][0][:]
        for fb in range(nfreq):
            data=MSE[aa][ncmp][fb+1][:]
            fmin=freq[fb]
            fmax=freq[fb+1]

            # small window smoothing
            npt=int(winlen[fb]/dt)+1
            half_npt=int(npt/2)
            arr=np.zeros(int(npt))
            for jj in range(0, (npts)):
                if jj < half_npt:
                    arr=data[jj: jj+half_npt]
                elif jj > (npts)-half_npt:
                    arr=data[jj: npts]
                else:
                    arr=data[jj-half_npt : jj+half_npt]
                msv[aa][ncmp][fb+1][jj] = msValue(arr,len(arr))
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
        sym = 0.5 * msv_mean[aa][fb+1][indx:] + 0.5 * np.flip(msv_mean[aa][fb+1][: indx + 1], axis=0)
        data_sym[fb]=sym
    fmsv_mean[aa]=[msv[aa][0][0][indx:],data_sym[0],data_sym[1],data_sym[2]]
    # -- plot
    plot_envelope(comp_arr,freq,msv[aa],msv_mean[aa],fname[aa],vdist[aa])
    plot_fmsv_waveforms(freq,fmsv_mean[aa],fname[aa])


# --- measuring coda window
twindow=[5,6,7,14,15,16]
cvel=[2.6, 2.0, 1.8]  
# phase velocity --> should be vary

x=np.zeros(2)  # mean_free_path search array
y=np.zeros(30)  # intrinsic_b search array

Esyn_temp=np.ndarray((len(x),len(y),npts//2+1))
Eobs_temp=np.ndarray((len(x),len(y),npts//2+1))

SSR_temppp=np.ndarray((len(x),len(y),len(twindow)))

SSR_final=np.ndarray((len(x),len(y)))
SSR_final[:][:]=0.

SSR=np.ndarray((nfreq,len(x),len(y)))
SSR[:][:][:]=0.

for fb in range(nfreq):
    fmin=freq[fb]
    fmax=freq[fb+1]
    c=cvel[fb]
    for aa in range(fnum):
        r=float(vdist[aa]) 

        # grid search in combination of mean_free_path and intrinsic_b
        Esyn_temp[:][:][:]=0.
        Eobs_temp[:][:][:]=0.

        for nfree in range(len(x)):
            mean_free= 0.4 + 0.2 *nfree
            x[nfree]=mean_free
            for nb in range(len(y)):
                intrinsic_b=0.02*nb
                y[nb]=intrinsic_b

                # calculate the Esyn and SSR for combination of mean_free_path and intrinsic_b
                for twn in range(npts//2+1):
                    tm=dt*twn
                    Eobs_temp[nfree][nb][twn]= fmsv_mean[aa][fb+1][twn]

                    s0=c**2 * tm**2 -r**2
                    if s0 <= 0:
                        #print(twn,tm,s0,tm-r/c)
                        continue

                    tmp=ESYN_RadiaTrans_onesta(mean_free, tm , r, c)
                    Esyn_temp[nfree][nb][twn]= tmp * math.exp(-1* intrinsic_b * tm)
                # using scalar factor for further fitting processes --> shape matters more than amplitude

                #### specific window --> find the scaling factor in the specific window
                for tsn in range(len(twindow)):
                    tsb=int(twindow[tsn]//dt)
                    SSR_temppp[nfree][nb][tsn]=0.
                    SSR_temppp[nfree][nb][tsn]=(math.log10(Eobs_temp[nfree][nb][tsb]) - math.log10(Esyn_temp[nfree][nb][tsb]))

                crap=np.mean(SSR_temppp[nfree][nb])
                Esyn_temp[nfree][nb]*=(10**crap)  # scale the Esyn

            #### specific window
            #### Calculate the SSR in the specific window
                for tsn in range(len(twindow)):
                    tsb=int(twindow[tsn]//dt)
                    tse=int((twindow[tsn]+1)//dt)
                    SSR_temp=0.
                    for twn in range(tsb,tse):
                        SSR_temp+=(math.log10(Eobs_temp[nfree][nb][twn]) - math.log10(Esyn_temp[nfree][nb][twn]))**2
                SSR_final[nfree][nb]+=SSR_temp
            #print("mean_free: %.2f " % mean_free,", intri_b %.2f " %  intrinsic_b,
            #      "mean(Eobs): %.2e" % np.mean(Eobs_temp[nfree][nb]),"mean(Esyn): %.2e" % np.mean(Esyn_temp[nfree][nb]),

            plot_fitting_curves(mean_free,y,fmsv_mean[aa][0][:],Eobs_temp[nfree],Esyn_temp[nfree],fname[aa],vdist[aa],twindow,fmin,fmax)
        SSR_final= SSR_final / (np.min(SSR_final[:][:]))
    SSR[fb]=SSR_final

# --- Final result 
result_intb=np.ndarray((nfreq))
result_mfp=np.ndarray((nfreq))
Eobs=np.ndarray((nfreq,npts//2+1))
Esyn=np.ndarray((nfreq,npts//2+1))
temppp=np.ndarray((len(twindow)))
aa=0
r=np.take(vdist[aa],0) #+0.000001
print("Final: ",sta_pair)
for fb in range(nfreq):
    fmin=freq[fb]
    fmax=freq[fb+1]

    loc=np.where(SSR[fb].T == np.amin(SSR[fb].T))
    locx=list(zip(loc[0], loc[1]))
    print("%4.2f-%4.2f Hz " % (fmin,fmax),"loc ",loc)
    ymin=y[loc[0]]
    xmin=x[loc[1]]
    print(" intrinsic_b %.2f " % ymin,"mean_free: %.2f " % xmin)
    result_intb[fb]=np.take(ymin,0)
    result_mfp[fb]=np.take(xmin,0)
    for twn in range(npts//2+1):
        tm=dt*twn
        s0=c**2 * tm**2 -r**2
        if s0 <= 0:
            continue
        Eobs[fb][twn]=fmsv_mean[aa][fb+1][twn]
        tmp=ESYN_RadiaTrans_onesta(result_mfp[fb], tm ,r, c)
        Esyn[fb][twn]= tmp * math.exp(-1 * result_intb[fb] * tm )

    for tsn in range(len(twindow)):
        tsb=int(twindow[tsn]//dt)
        temppp[tsn]=0.
        temppp[tsn]=(math.log10(Eobs[fb][tsb]) - math.log10(Esyn[fb][tsb]))

    crap=np.mean(temppp)
    Esyn[fb]*=(10**crap)
    plot_fitting_result(result_mfp[fb],result_intb[fb],fmsv_mean[aa][0][:],
                        Eobs[fb],Esyn[fb],sta_pair,vdist[aa],twindow,fmin,fmax)
#plot_grid_searching(sta_pair,freq,SSR,x,y)

# --- move figures to dics
root=os.getcwd()
output_imgdir = root+"/../figures"
if not os.path.exists(output_imgdir):
        os.makedirs(output_imgdir)
images= [f for f in os.listdir() if '.png' in f.lower()]
for image in images:
    shutil.move(image, output_imgdir)
