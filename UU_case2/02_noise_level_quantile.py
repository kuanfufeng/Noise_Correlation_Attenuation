#!/usr/bin/python
import sys
import datetime
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import gridspec

fd=sys.argv[1] #"LEVEL_ave30-30_HighFreq/"
sta=sys.argv[2]
fb0=sys.argv[3]
fb1=sys.argv[4]
fb2=sys.argv[5]
# Read txt
ftxt=fd+"/NOISE_LEVEL_"+sta+"_"+fb0+".txt"
fi = pd.read_csv(ftxt)
f0=fi['Fband_(Hz)'].astype(str)[0]
dtime=fi['ENV'].astype(str).tolist()
level=fi['Noise_level'].tolist()
ampend=fi['Amp(et)'].tolist()
bt=fi['window_bt'].tolist()
et=fi['window_et'].tolist()

# Read txt
ftxt=fd+"/NOISE_LEVEL_"+sta+"_"+fb1+".txt"
fi = pd.read_csv(ftxt)
f1=fi['Fband_(Hz)'].astype(str)[0]
dtime1=fi['ENV'].astype(str).tolist()
level1=fi['Noise_level'].tolist()
ampend1=fi['Amp(et)'].tolist()
bt1=fi['window_bt'].tolist()
et1=fi['window_et'].tolist()

# Read txt
ftxt=fd+"/NOISE_LEVEL_"+sta+"_"+fb2+".txt"
fi = pd.read_csv(ftxt)
f2=fi['Fband_(Hz)'].astype(str)[0]
dtime2=fi['ENV'].astype(str).tolist()
level2=fi['Noise_level'].tolist()
ampend2=fi['Amp(et)'].tolist()
bt2=fi['window_bt'].tolist()
et2=fi['window_et'].tolist()


t=[datetime.datetime.strptime(d,'%Y%m%d').strftime('%Y-%m-%d') for d in dtime]
t0=pd.to_datetime(t)
t=[datetime.datetime.strptime(d,'%Y%m%d').strftime('%Y-%m-%d') for d in dtime1]
t1=pd.to_datetime(t)
t=[datetime.datetime.strptime(d,'%Y%m%d').strftime('%Y-%m-%d') for d in dtime2]
t2=pd.to_datetime(t)
# plot to one figure
'''
fig = plt.figure(figsize=(16, 5)) 
gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1]) 
ax0 = plt.subplot(gs[0])
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,7]))
plt.grid(True)
plt.xticks(rotation=90)
ax0.set_xlim(np.datetime64('2006-01-01'), np.datetime64('2022-12-31'))
ax0.plot(np.array(t0),np.array(ampend), '-',c='orange',label="Freq "+f0+" Hz")
ax0.plot(np.array(t1),np.array(ampend1), '-',c='r',label="Freq "+f1+" Hz")
ax0.plot(np.array(t2),np.array(ampend2), '-',c='b',label="Freq "+f2+" Hz")

ax0.set_title('Noise level of Mean-Squared vaule at '+sta+' (endtime at 7 wvl)')
ax0.set_ylabel('Noise Amplitude')
ax0.set_xlabel('Time')
plt.legend()


ax1 = plt.subplot(gs[1])
ax1.hist(ampend,bins = 100, alpha = 0.35, color = 'orange',label='Amp-et at '+f0, orientation="horizontal")
ax1.hist(ampend1,bins = 100, alpha = 0.35, color = 'red',label='Amp-et at '+f1, orientation="horizontal")
ax1.hist(ampend2,bins = 100, alpha = 0.25, color = 'blue',label='Amp-et at '+f2, orientation="horizontal")
plt.axhline(np.median(ampend), color='orange' , linestyle='dashed', linewidth=1,label='median of '+f0)
plt.axhline(np.mean(ampend),   color='orange' , linestyle='solid', linewidth=1,label='mean of '+f0)
plt.axhline(np.median(ampend1), color='red' , linestyle='dashed', linewidth=1,label='median of '+f1)
plt.axhline(np.mean(ampend1),   color='red' , linestyle='solid', linewidth=1,label='mean of '+f1)
plt.axhline(np.median(ampend2), color='cyan' , linestyle='dashed', linewidth=1,label='median of '+f2)
plt.axhline(np.mean(ampend2),   color='cyan' , linestyle='solid', linewidth=1,label='mean of '+f2)
plt.title('Station '+sta+' (endtime at 7 wvl)' )
plt.ylabel('Noise amplitude at 7*wvl')
plt.xlabel('Count')
plt.legend()

plt.savefig('Noiselevel_'+sta+'.png')
'''
QT1=0.01
QT5=0.05
QTV1_F0=np.percentile(ampend, QT1 * 100)
QTV5_F0=np.percentile(ampend, QT5 * 100)
QTV1_F1=np.percentile(ampend1, QT1 * 100)
QTV5_F1=np.percentile(ampend1, QT5 * 100)
QTV1_F2=np.percentile(ampend2, QT1 * 100)
QTV5_F2=np.percentile(ampend2, QT5 * 100)


# plot to one figure

plt.tight_layout()
fig2 = plt.figure(constrained_layout=True,figsize=(16,12))
gs = gridspec.GridSpec(nrows=3, ncols=2, width_ratios=[2, 1],height_ratios=[1,1,1], figure=fig2)

f2_ax0= fig2.add_subplot(gs[0,0])
f2_ax0.set_xlim(np.datetime64('2006-01-01'), np.datetime64('2022-12-31'))
f2_ax0.plot(np.array(t0),np.array(ampend), '-',c='orange',label="Freq "+f0+" Hz")
f2_ax0.set_title('Noise level of Mean-Squared vaule at '+sta+' (endtime at 7 wvl)')
f2_ax0.set_ylabel('Noise amplitude at 7*wvl')
f2_ax0.set_xlabel('Time')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,7]))
plt.xticks(rotation=50)

f2_ax1 = fig2.add_subplot(gs[0,1])
f2_ax1 .hist(ampend,bins = 100, alpha = 0.35, color = 'orange',label='Endtime at '+f0, orientation="horizontal")

f2_ax2 = plt.subplot(gs[1,0])
f2_ax2.set_xlim(np.datetime64('2006-01-01'), np.datetime64('2022-12-31'))
f2_ax2.plot(np.array(t1),np.array(ampend1), '-',c='red',label="Freq "+f1+" Hz")
f2_ax2.set_title('Noise level of Mean-Squared vaule at '+sta+' (endtime at 7 wvl)')
f2_ax2.set_ylabel('Noise amplitude at 7*wvl')
f2_ax2.set_xlabel('Time')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,7]))
plt.xticks(rotation=50)

f2_ax3 = fig2.add_subplot(gs[1,1])
f2_ax3.hist(ampend1,bins = 100, alpha = 0.35, color = 'red',label='Endtime at '+f1, orientation="horizontal")

f2_ax4 = fig2.add_subplot(gs[2,0])
f2_ax4.set_xlim(np.datetime64('2006-01-01'), np.datetime64('2022-12-31'))
f2_ax4.plot(np.array(t2),np.array(ampend2), '-',c='b',label="Freq "+f2+" Hz")
f2_ax4.set_title('Noise level of Mean-Squared vaule at '+sta+' (endtime at 7 wvl)')
f2_ax4.set_ylabel('Noise amplitude at 7*wvl')
f2_ax4.set_xlabel('Time')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,7]))
plt.xticks(rotation=50)

f2_ax5 = fig2.add_subplot(gs[2,1])
f2_ax5.hist(ampend2,bins = 100, alpha = 0.25, color = 'blue',label='Endtime at '+f2, orientation="horizontal")

f2_ax1.axhline(QTV1_F0,   color='k' , linestyle='solid', linewidth=1,label='Quantile 0.01 of '+f0)
f2_ax1.axhline(QTV5_F0,   color='k' , linestyle='--', linewidth=1,label='Quantile 0.05 of '+f0)

f2_ax3.axhline(QTV1_F1,   color='k' , linestyle='solid', linewidth=1,label='Quantile 0.01 of '+f1)
f2_ax3.axhline(QTV5_F1,   color='k' , linestyle='--', linewidth=1,label='Quantile 0.05 of '+f1)

f2_ax5.axhline(QTV1_F2,   color='k' , linestyle='solid', linewidth=1,label='Quantile 0.01 of '+f2)
f2_ax5.axhline(QTV5_F2,   color='k' , linestyle='--', linewidth=1,label='Quantile 0.05 of '+f2)

f2_ax0.grid(True)
f2_ax2.grid(True)
f2_ax4.grid(True)
f2_ax0.legend()
f2_ax1.legend()
f2_ax2.legend()
f2_ax3.legend()
f2_ax4.legend()
f2_ax5.legend()

plt.title('Station '+sta+' (endtime at Quantile 0.01 or 0.05)')
plt.ylabel('Noise amplitude at 7*wvl')
plt.xlabel('Count')

plt.savefig('NoiseQTsub_'+sta+'.png')

# --- print results to file
ftxt="Noise_quantile.txt"
f=open(ftxt, "a")
f.write('stan,qt001('+f0+')'+',qt005('+f0+')'+',qt001('+f1+')'+',qt005('+f1+')'+',qt001('+f2+')'+',qt005('+f2+')\n')
line=sta+",%.6f"%QTV1_F0+",%.6f"%QTV5_F0+",%.6f"%QTV1_F1+",%.6f"%QTV5_F1+",%.6f"%QTV1_F2+",%.6f"%QTV5_F2+"\n"
print(line)
f.write(str(line))
f.close()
