#!/usr/bin/python
import sys
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import gridspec

noisefn=sys.argv[1]  #"WINDOW_ave30-30_2wvl-QT5_HF/"
sta=sys.argv[2]
fblst=sys.argv[3]


freq=np.loadtxt(fblst,dtype=str)
flen=len(freq)

ftxt=noisefn
fi = pd.read_csv(ftxt)

# --- print results to file
ftxt2="ENDTIME_QT5_"+sta+".txt"
f2=open(ftxt2, "w")
line="stan,freq,tbeg,tend\n"
f2.write(line)

plt.tight_layout()
fig2 = plt.figure(constrained_layout=True,figsize=(10,30))
gs = gridspec.GridSpec(nrows=int(flen), ncols=2, width_ratios=[2, 1],
                    height_ratios=np.linspace(1,1, num=flen).tolist(), figure=fig2)


for k in range(flen):

    pdk = fi[fi['Fband_(Hz)'] == freq[k]]
    ff = freq[k]
    dtime = pdk['ENV'].astype(str).tolist()
    level = pdk['Noise_qt5'].tolist()
    ampend = pdk['Amp(et)'].tolist()
    bt = pdk['window_bt'].tolist()
    et = pdk['window_et'].tolist()

    t = [datetime.datetime.strptime(d,'%Y%m%d').strftime('%Y-%m-%d') for d in dtime]
    t0 = pd.to_datetime(t)
    QT1 = 0.01
    QT5 = 0.05
    QTV1_F0 = np.percentile(et, QT1 * 100)
    QTV5_F0 = np.percentile(et, QT5 * 100)

    line = sta+",%s" % ff+",%.6f" % np.min(bt)+",%.6f" % QTV5_F0+"\n"
    f2.write(str(line))

    # plot to one figure
    f2_ax0 = fig2.add_subplot(gs[k,0])
    f2_ax0.set_xlim(np.datetime64('2006-01-01'), np.datetime64('2022-12-31'))
    f2_ax0.plot(np.array(t0), np.array(et), ls="-", c='blue', linewidth=0.5)
    f2_ax0.plot(np.array(t0), np.array(et), ls="", marker='.', c='blue')
    f2_ax0.plot([t0[0], t0[-1]], [QTV5_F0, QTV5_F0], '-', c='k',
                label="End time")
    f2_ax0.plot([t0[0], t0[-1]], [np.min(bt),np.min(bt)], '-', c='green',
                label="Begin time")
    f2_ax0.set_title('Window endtime of '+sta+' (quantile 0.05 of noise level)')
    f2_ax0.set_ylabel('Window endtime (sec)')
    f2_ax0.set_xlabel('Time')
    f2_ax0.grid(True)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    plt.gca().xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 7]))
    plt.xticks(rotation=50)

    f2_ax1 = fig2.add_subplot(gs[k, 1])
    f2_ax1.hist(et, bins=50, alpha=0.35, color='blue',
                label='End time', orientation="horizontal")
    f2_ax1.axhline(np.min(bt),   color='green', linestyle='solid',
                   linewidth=1, label='Begin Time')
    f2_ax1.axhline(QTV5_F0,   color='k', linestyle='solid',
                   linewidth=1, label='Quantile 0.05')

    f2_ax0.legend(loc='upper right')
    f2_ax1.legend(loc='upper right')
    plt.title('Window endtime of '+sta+' \n(quantile 0.05 of noise level)')
    plt.xlabel('Count')

plt.savefig('WindowEndtime_'+sta+'.png')
f2.close()
