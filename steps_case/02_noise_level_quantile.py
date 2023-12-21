#!/usr/bin/python
import sys
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import gridspec

fd = sys.argv[1]
sta = sys.argv[2]
fblst = sys.argv[3]

freq = np.loadtxt(fblst, dtype=str)
flen = len(freq)
ftxt = fd+"/NOISE_LEVEL_"+sta+".txt"
fi = pd.read_csv(ftxt)

plt.tight_layout()
fig2 = plt.figure(constrained_layout=True, figsize=(10, 30))
gs = gridspec.GridSpec(nrows=int(flen), ncols=2, width_ratios=[2, 1],
                    height_ratios=np.linspace(1, 1, num=flen).tolist(),
                    figure=fig2)

# --- print results to file
ftxt2 = "Noise_quantile_"+sta+".txt"
fq = open(ftxt2, 'w')
line = "stan,freq,qt001,qt005\n"
fq.write(str(line))

for k in range(flen):

    pdk = fi[fi['Fband_(Hz)'] == freq[k]]
    ff = freq[k]
    dtime = pdk['ENV'].astype(str).tolist()
    level = pdk['Noise_level'].tolist()
    ampend = pdk['Amp(et)'].tolist()
    bt = pdk['window_bt'].tolist()
    et = pdk['window_et'].tolist()

    t = [datetime.datetime.strptime(d, '%Y%m%d').strftime('%Y-%m-%d')
         for d in dtime]
    t0 = pd.to_datetime(t)
    QT1 = 0.01
    QT5 = 0.05
    QTV1_F0 = np.percentile(ampend, QT1 * 100)
    QTV5_F0 = np.percentile(ampend, QT5 * 100)

    line = sta+",%s" % ff+",%.6f" % QTV1_F0+",%.6f" % QTV5_F0+"\n"
    print(line)
    fq.write(str(line))

    f2_ax0 = fig2.add_subplot(gs[k,0])
    f2_ax0.set_xlim(np.datetime64('2006-01-01'), np.datetime64('2022-12-31'))
    f2_ax0.plot(np.array(t0), np.array(ampend), ls="-", c='blue', linewidth=0.5)
    f2_ax0.plot(np.array(t0), np.array(ampend), ls="", marker='.', c='blue')
    f2_ax0.set_title('Noise level of MSV at '+sta+' in '+ff+'Hz (endtime at 7 wvl)')
    f2_ax0.set_ylabel('Noise amplitude at 7*wvl')
    f2_ax0.set_xlabel('Time')
    f2_ax0.grid(True)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    plt.gca().xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 7]))
    plt.xticks(rotation=50)

    f2_ax1 = fig2.add_subplot(gs[k, 1])
    f2_ax1.hist(ampend, bins=50, alpha=0.35, color='blue',
                label='Endtime at '+ff+' Hz', orientation="horizontal")
    f2_ax1.axhline(QTV1_F0,   color='k', linestyle='solid', linewidth=1,
                   label='Quantile 0.01')
    f2_ax1.axhline(QTV5_F0,   color='k', linestyle='--', linewidth=1,
                   label='Quantile 0.05')

    f2_ax1.legend()
    plt.title('Station '+sta+' endtime amplitude')
    # plt.ylabel('Noise amplitude at 7*wvl')
    plt.xlabel('Count')
plt.savefig('NoiseQTsub_'+sta+'.png')
fq.close()
