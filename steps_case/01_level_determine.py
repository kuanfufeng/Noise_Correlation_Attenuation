#!/usr/bin/python
import os 
import sys
import glob
import numpy as np
from Esyn_func import *
from plotting import *
from obspy.signal.filter import bandpass


path = sys.argv[1]
sta_pair = sys.argv[2]

# compoent
# comp_arr = ["ZN", "ZE","NE", "EN","EZ","NZ"]
comp_arr = ["ZN_ave30-30", "ZE_ave30-30", "NE_ave30-30", "EN_ave30-30",
            "EZ_ave30-30", "NZ_ave30-30"]
num_cmp = len(comp_arr)

lag = 60
samp = 20
leng = int(lag*samp*2+1)
npts = leng
# print("\n# Lag-time (sec): ",lag,"\n
# # Sampling rate: ",samp,"\n
# # Data length (npts): ",leng)
# targeted frequency band for waveform monitoring
freq = [0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.6, 2, 4, 6, 8]
winlen = [8, 8, 6, 4, 4, 2, 2, 2, 1, 1]
# winlen=[8,4,2]
nfreq = len(freq) - 1
# print("# Target Frequency Bands: ", freq)
# print("# Smoothing window length: ", winlen)

# --- print results to file
ftxt = "NOISE_LEVEL_"+sta_pair+".txt"
f0 = open(ftxt, "w")

line = "ENV,Fband_(Hz),Noise_level,Amp(et),window_bt,window_et,file_name\n"
f0.write(line)

evd = np.loadtxt("evtlst", dtype=str)
for nev in range(len(evd)):
    ev = evd[nev]
    inpf = sta_pair+"_"+sta_pair+"_"+ev+"_ave30-30.h5"
    sfiles = sorted(glob.glob(os.path.join(path, inpf)))
    # print("# Read in file: ",sfiles,"\n # Station Pair: ",sta_pair)
    fnum = len(sfiles)
    stackf = np.ndarray((fnum, num_cmp, 2, leng))
    # print("Array shape: ",stackf.shape)
    vdist = np.zeros((fnum, 1))  # S-R distance array
    fname = []                  # file name array
    aa = 0
    # loop through each station-pair
    for sfile in sfiles:
        if not os.path.exists(sfile):
            continue
        ncmp = 0
        fname.append(sta_pair+"_"+ev)
        for ccomp in comp_arr:
            # print(aa, sfile, ccomp)
            # read stacked waveforms
            if (read_pyasdf(sfile, ccomp) == None):
                continue
            # read waveform from pyasdf
            dist, dt, tvec, sdata = read_pyasdf(sfile, ccomp)

            stackf[aa][ncmp] = [tvec, sdata]
            vdist[aa] = 0   # dist
            ncmp = ncmp+1
        aa = aa+1
    fnum = len(fname)
    # print("# total files number: ",fnum)
    if fnum == 0:
        continue
    # --- filtering and mean-squared values
    indx = npts // 2                      # half-side number of points
    MSE = np.ndarray((fnum, num_cmp, nfreq+1, npts))
    # filtered two-side averaged stack CF

    for aa in range(fnum):
        dafbp = np.ndarray((nfreq, npts))

        for ncmp in range(len(comp_arr)):
            ccomp = comp_arr[ncmp]
            # print(fname[aa],ccomp)

            MSE[aa][ncmp][0] = stackf[aa][ncmp][0]
            for fb in range(nfreq):
                fmin = freq[fb]
                fmax = freq[fb+1]
                tt = np.arange(0, npts) * dt
                data = stackf[aa][ncmp][1]
                dafbp[fb] = bandpass(data, fmin, fmax, int(1 / dt),
                                     corners=4, zerophase=True)
                MSE[aa][ncmp][fb+1] = dafbp[fb][:]

    # --- calculate average msv
    msv = np.zeros((fnum, num_cmp, nfreq+1, npts))
    msv_mean = np.zeros((fnum, nfreq+1, npts))

    data_sym = np.zeros((nfreq, indx+1))  # two-side averaged stack CF
    fmsv_mean = np.zeros((fnum, nfreq+1, indx+1))

    # --- noise level
    level = np.ndarray((fnum, nfreq, 1))
    twinbe = np.ndarray((fnum, nfreq, 2))

    for aa in range(fnum):

        for ncmp in range(len(comp_arr)):
            ccomp = comp_arr[ncmp]
            msv[aa][ncmp][0] = MSE[aa][ncmp][0][:]
            for fb in range(nfreq):
                data = MSE[aa][ncmp][fb+1][:]
                fmin = freq[fb]
                fmax = freq[fb+1]

                # small window smoothing
                para = {'winlen': winlen[fb], 'dt': dt, 'npts': len(data)}
                msv[aa][ncmp][fb+1] = get_smooth(data, para)
                # self-normalized
                msv[aa][ncmp][fb+1] = msv[aa][ncmp][fb+1]/np.max(msv[aa][ncmp][fb+1])

        msv_mean[aa][0] = msv[aa][0][0][:]
        for fb in range(nfreq):
            fmin = freq[fb]
            fmax = freq[fb+1]

            for ncmp in range(len(comp_arr)):
                msv_mean[aa][fb+1] += msv[aa][ncmp][fb+1][:]
            msv_mean[aa][fb+1] = msv_mean[aa][fb+1]/len(comp_arr)
            # --- fold
            sym = get_symmetric(msv_mean[aa][fb+1], indx)
            data_sym[fb] = sym

            twinbe[aa][fb][0] = (1/fmin)*4
            twinbe[aa][fb][1] = twinbe[aa][fb][0]+(1/fmin)*3
            ampend = sym[int(twinbe[aa][fb][1]/dt)]
            for pt in range(len(sym)):
                if (sym[pt] < float(ampend)):
                    level[aa][fb] = sym[pt]
                    break
            # print(fname[aa],"Fband_(Hz) %.1f-%.1f"%(fmin,fmax),"noise_level %.6f"%level[aa][fb],"ampend %.6f"%ampend, "coda_window %.1f %.1f "%(twinbe[aa][fb][0],twinbe[aa][fb][1]) )
            line = str(ev)+",%.1f-%.1f" % (fmin,fmax) \
                + ",%.6f" % (level[aa][fb]) \
                + ",%.6f" % (ampend) \
                + ",%.2f,%.2f" % (twinbe[aa][fb][0],twinbe[aa][fb][1]) \
                + ",%s" % sfiles[aa]+"\n"
            # print(line)
            f0.write(line)


# Closing the opened file
f0.close()
