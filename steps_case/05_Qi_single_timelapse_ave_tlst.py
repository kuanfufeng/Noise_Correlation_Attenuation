#!/usr/bin/python
import os
import sys
import glob
import numpy as np
import pandas as pd
from Esyn_func import *
from plotting import *
from obspy.signal.filter import bandpass


boots = "ave30-30"
path = sys.argv[1]
sta_pair = sys.argv[2]
fblst = sys.argv[3]
winfn = sys.argv[4]

# targeted frequency band for waveform monitoring
freq = np.loadtxt(fblst, dtype=float)
flen = len(freq)

twin_pd = pd.read_csv(winfn)
tbeg = twin_pd["tbeg"].tolist()
tend = twin_pd["tend"].tolist()
for k in range(len(freq)):
    if (tend[k] <= tbeg[k]):
        tend[k] = tbeg[k]+(1./freq[k])*2

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
# smoothing window lengths corresponding to the frequency bands
winlen = [8, 8, 6, 4, 4, 2, 2, 2, 1, 1]
nfreq = len(freq) - 1
# print("# Target Frequency Bands: ", freq)
# print("# Smoothing window length: ", winlen)

# --- print results to file
ftxt = "OUTPUT_Qi_"+sta_pair+".txt"
file = open(ftxt, "w")
line = "ENV,fband,fcen,intrinsic_b,intrinsic_Q,tbeg,tend,filename\n"
file.write(line)

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
            vdist[aa] = 0
            ncmp = ncmp+1
        # -- plot
        plot_waveforms(num_cmp, stackf[aa], fname[aa], comp_arr)
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
            # -- plot
            # plot_filtered_waveforms(freq, stackf[aa][ncmp][0],
            # dafbp, fname[aa], ccomp)

    # --- calculate average msv
    msv = np.zeros((fnum, num_cmp, nfreq+1, npts))
    msv_mean = np.zeros((fnum, nfreq+1, npts))

    data_sym = np.zeros((nfreq, indx+1))  # two-side averaged stack CF
    fmsv_mean = np.zeros((fnum, nfreq+1, indx+1))

    # --- noise level
    level = np.ndarray((fnum, nfreq, 1))
    ratio = 5
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
        fmsv_mean[aa][0] = msv[aa][0][0][indx:]
        for fb in range(nfreq):
            fmin = freq[fb]
            fmax = freq[fb+1]

            for ncmp in range(len(comp_arr)):
                msv_mean[aa][fb+1] += msv[aa][ncmp][fb+1][:]
            msv_mean[aa][fb+1] = msv_mean[aa][fb+1]/len(comp_arr)
            # --- fold
            sym = get_symmetric(msv_mean[aa][fb+1], indx)
            fmsv_mean[aa][fb+1] = get_symmetric(msv_mean[aa][fb+1], indx)

            Val_mad = mad(sym)
            level[aa][fb] = Val_mad*ratio

            twinbe[aa, fb, :] = tbeg[fb], tend[fb]
        # -- plot
        plot_envelope(comp_arr, freq, msv[aa], msv_mean[aa],
                      fname[aa], vdist[aa])
        plot_fmsv_waveforms(freq, fmsv_mean[aa], fname[aa],
                            level[aa], twinbe[aa])

    # --- measuring coda window
    cvel = 2.5    # Rayleigh wave velocities over the freqency bands
    mfpx = np.zeros(1)        # mean_free_path search array
    intby = np.zeros(60)      # intrinsic_b search array
    # getting the sum of squared residuals (SSR) between Eobs and Esyn
    SSR_final = np.zeros((len(mfpx), len(intby)))
    SSR = np.zeros((nfreq, len(mfpx), len(intby)))

    for fb in range(nfreq):
        fmin = freq[fb]
        fmax = freq[fb+1]
        c = cvel
        SSR_final[:][:] = 0.
        vdist[:] = 0.000001  # To avoid zero value at denominator

        # parameters for getting the sum of squared residuals (SSR)
        # between Eobs and Esyn
        para = {'fb': fb, 'vdist': vdist, 'npts': npts, 'dt': dt, 'cvel': c,
                'mfp': mfpx, 'intb': intby, 'twin': twinbe, 'fmsv': fmsv_mean}
        # call function get_SSR
        SSR_final, mfpx, intby = get_SSR(fnum, para)

        SSR[fb] = SSR_final

    # getting the optimal value from the SSR
    result_intb = np.ndarray((nfreq))
    result_mfp = np.ndarray((nfreq))

    Eobs = np.ndarray((nfreq, npts//2+1))
    Esyn = np.ndarray((nfreq, npts//2+1))
    aa = 0
    r = np.take(vdist[aa], 0)+0.000001
    for fb in range(nfreq):
        fmin = freq[fb]
        fmax = freq[fb+1]

        # parameters for getting optimal value from
        # the sum of squared residuals (SSR) between Eobs and Esyn
        para = {'fb': fb, 'fmin': fmin, 'fmax': fmax, 'vdist': vdist,
                'npts': npts, 'dt': dt, 'cvel': c, 'filenum': aa,
                'mfp': mfpx, 'intb': intby, 'twin': twinbe,
                'fmsv': fmsv_mean, 'SSR': SSR, 'sta': sta_pair}
        # call function get_optimal
        result_intb[fb], result_mfp[fb], Eobs[fb], Esyn[fb] = get_optimal(fnum, para)

        # --- print results to file
        REb = np.take(result_intb[fb], 0)
        Q = (2.0*np.pi*((fmax+fmin)/2.0))/REb
        # line = "ENV,fband,intrinsic_b,intrinsic_Q,tbeg,tend,file\n"
        line = str(ev)+",%.1f-%.1f" % (fmin, fmax)+",%.4f" % ((fmax+fmin)/2.0)\
                    + ",%.6f" % (REb) \
                    + ",%.6f" % (Q) \
                    + ",%.2f,%.2f" % (twinbe[aa][fb][0], twinbe[aa][fb][1]) \
                    + ",%s" % sfiles[aa]+"\n"
        # Writing the file
        file.write(line)

        plot_fitting_result(result_mfp[fb], result_intb[fb], Q,
                            fmsv_mean[aa][0][:], Eobs[fb], Esyn[fb],
                            sta_pair, str(ev),
                            vdist[aa], twinbe[aa][fb], fmin, fmax)

# Closing the opened file
file.close()
