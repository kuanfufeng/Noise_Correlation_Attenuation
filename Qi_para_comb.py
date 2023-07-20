import os
import sys
import glob
import obspy
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import pyasdf
import scipy
import math

from obspy.signal.filter import bandpass

## Select data 
path = "./"
data_path = os.path.join(path, "STACK_BH_Rainier")
sfiles = sorted(glob.glob(os.path.join(data_path, '*PANH*/*.h5')))
#sfiles = sorted(glob.glob(os.path.join(data_path, '*/*.h5')))
print(sfiles)

## 
comp_arr = ["ZZ", "ZR","ZT","RZ","TZ"]              # component
#comp_arr = ["ZZ"]
num_cmp=len(comp_arr)
fnum=len(sfiles)
print(fnum)

## sampling rate (samp) and max lag time (lag)
lag=200
samp=20
leng=int(lag*samp*2+1)
npts=leng
dt=1/samp

# ----- read files from asdf -----
stackf=np.ndarray((fnum,num_cmp,2,leng))
print(stackf.shape)
vdist=np.zeros((fnum,1))  # S-R distance array
fname=[]                  # file name array

aa=0
# loop through each station-pair
for sfile in sfiles:
    ncmp=0

    # skip auto-correlation
    if ( sfile.split("/")[3].split("_")[0] == sfile.split("/")[3].split("_")[1].split(".h5")[0] ):
        continue
    fname.append(sfile.split("/")[3].split(".h5")[0])
    for ccomp in comp_arr:
        print(aa, sfile, ccomp)
        # read stacked waveforms
        if ( read_pyasdf(sfile,ccomp) == None):
            continue
        dist, tvec,sdata = read_pyasdf(sfile,ccomp) # read waveform from pyasdf
        stackf[aa][ncmp]=[tvec,sdata]
        vdist[aa]=dist
        ncmp=ncmp+1
    plot_waveforms(num_cmp,stackf[aa],fname[aa],comp_arr)

    aa=aa+1
fnum=len(fname)

# ----- Waveform filtering -----
freq = [0.5, 1, 2,4]  # targeted frequency band for waveform monitoring
nfreq = len(freq) - 1

indx = npts // 2                      # half-side number of points 
MSE=np.ndarray((fnum,num_cmp,nfreq+1,indx+1)) # filtered two-side averaged stack CF
data_sym=np.ndarray((nfreq+1,indx+1)) # two-side averaged stack CF

for aa in range (fnum):
    dafbp=np.ndarray((nfreq,npts))
      
    for ncmp in range (len(comp_arr)):
        ccomp=comp_arr[ncmp]
        
        for fb in range(nfreq):
            fmin=freq[fb]
            fmax=freq[fb+1]
            tt = np.arange(0, npts) * dt
            data = stackf[aa][ncmp][1]
            dafbp[fb] = bandpass(data, fmin, fmax, int(1 / dt), corners=4, zerophase=True)
        
            # stack positive and negative lags  
            sym = 0.5 * dafbp[fb][indx:] + 0.5 * np.flip(dafbp[fb][: indx + 1], axis=0)
            data_sym[fb]=sym
        MSE[aa][ncmp]=[stackf[aa][ncmp][0][indx:],data_sym[0],data_sym[1],data_sym[2]] 
        print(MSE[aa][ncmp].shape)
        plot_filtered_waveforms(freq,stackf[aa][ncmp][0],dafbp,MSE[aa][ncmp],fname[aa])


# ----- mean-squared envelope -----
msv=np.ndarray((fnum,num_cmp,nfreq+1,indx+1))
msv_mean=np.ndarray((fnum,nfreq+1,indx+1))
msv[:][:][:][:]=0.
msv_mean[:][:][:]=0.
#print(msv.shape)

winlen=[4,2,1]   # smoothing window lengths corresponding to the frequency bands
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
            for jj in range(0, (npts//2)):
                if jj < half_npt:
                    arr=data[jj: jj+half_npt]
                elif jj > (npts//2)-half_npt:
                    arr=data[jj: npts]
                else:
                    arr=data[jj-half_npt : jj+half_npt]
                msv[aa][ncmp][fb+1][jj] = msValue(arr,len(arr))
        
            msv[aa][ncmp][fb+1]=msv[aa][ncmp][fb+1]/np.max(msv[aa][ncmp][fb+1])  # self-normalized 
    
    msv_mean[aa][0]=msv[aa][0][0][:]
    for fb in range(nfreq):
        fmin=freq[fb]
        fmax=freq[fb+1]
        for ncmp in range(len(comp_arr)):
            msv_mean[aa][fb+1]+=msv[aa][ncmp][fb+1][:]
        msv_mean[aa][fb+1]=msv_mean[aa][fb+1]/len(comp_arr)
        
    plot_envelope(comp_arr,freq,msv[aa],msv_mean[aa],fname[aa],vdist[aa])

# ----- window setting -----
#twindow=[5.24, 5.74, 6.24, 36.1, 37.1, 37.6] # time windows   --> should be determined depended on station pair
twindow=[20, 21, 22, 30, 31, 32]
cvel=[2.6, 2.0, 1.8]                                     # phase velocity --> should be vary

x=np.zeros(120)    # mean_free_path search array
y=np.zeros(16)     # intrinsic_b search array

Esyn_temp=np.ndarray((len(x),len(y),npts//2+1))
Eobs_temp=np.ndarray((len(x),len(y),npts//2+1))

# ----- grid-searching mean_free_path and intrinsic_b -----
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
    #for aa in range(0,1):
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
                    Eobs_temp[nfree][nb][twn]= msv_mean[aa][fb+1][twn]

                    s0=c**2 * tm**2 -r**2
                    if s0 <= 0:
                        #print(twn,tm,s0,tm-r/c)
                        continue

                    tmp=ESYN_RadiaTrans(mean_free, tm , r, c)
                    Esyn_temp[nfree][nb][twn]= tmp * math.exp(-1* intrinsic_b * tm)
                # using scalar factor for further fitting processes --> shape matters more than amplitude

                #### specific window --> find the scaling factor in the specific window
                for tsn in range(len(twindow)):
                    tsb=int(twindow[tsn]//dt)
                    SSR_temppp[nfree][nb][tsn]=0.
                    SSR_temppp[nfree][nb][tsn]=(math.log10(Eobs_temp[nfree][nb][tsb]) - math.log10(Esyn_temp[nfree][nb][tsb]))

                #print("tsn",tsn,tsb,tse,"Eobs %.2f" % Eobs_temp[nfree][nb][tsn],"Esyn %.2e" % Esyn_temp[nfree][nb][tsn],
                #    "SSR_temp %.2f" % SSR_temppp[nfree][nb][tsn],"crap",crap,res_max)
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
            #      "crap %.2f" % crap, "res_max %.2f" % res_max)
            
            ### --- plotting the Eobs and Esyn curves --> will take time to plot all out
            #plot_fitting_curves(nfree,y,msv_mean[aa][0][:],Eobs_temp[nfree],Esyn_temp[nfree],fname[aa])
        #print("mean_free: %.2f " % mean_free,", intrinsic_b %.2f " %  intrinsic_b,"SSR: %.4f" % SSR_temp)
        SSR_final= SSR_final / (np.min(SSR_final[:][:]))
    SSR[fb]=SSR_final

## optimal fit
plot_grid_searching(freq,SSR)

