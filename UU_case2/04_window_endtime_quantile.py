#!/usr/bin/python
import sys
import datetime
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import gridspec


fd=sys.argv[1]  #"WINDOW_ave30-30_2wvl-QT5_HF/"
sta=sys.argv[2]
fb0=sys.argv[3]
fb1=sys.argv[4]
fb2=sys.argv[5]

# Read txt
ftxt=fd+"/NOISE_WINDOW_"+sta+"_"+fb0+".txt"
fi = pd.read_csv(ftxt)
f0=fi['Fband_(Hz)'].astype(str)[0]
dtime=fi['ENV'].astype(str).tolist()
level=fi['Noise_qt5'].tolist()
ampend=fi['Amp(et)'].tolist()
bt=fi['window_bt'].tolist()
et=fi['window_et'].tolist()

# Read txt
ftxt=fd+"/NOISE_WINDOW_"+sta+"_"+fb1+".txt"
fi = pd.read_csv(ftxt)
f1=fi['Fband_(Hz)'].astype(str)[0]
dtime1=fi['ENV'].astype(str).tolist()
level1=fi['Noise_qt5'].tolist()
ampend1=fi['Amp(et)'].tolist()
bt1=fi['window_bt'].tolist()
et1=fi['window_et'].tolist()

# Read txt
ftxt=fd+"/NOISE_WINDOW_"+sta+"_"+fb2+".txt"
fi = pd.read_csv(ftxt)
f2=fi['Fband_(Hz)'].astype(str)[0]
dtime2=fi['ENV'].astype(str).tolist()
level2=fi['Noise_qt5'].tolist()
ampend2=fi['Amp(et)'].tolist()
bt2=fi['window_bt'].tolist()
et2=fi['window_et'].tolist()

# Quantile part
QT1=0.01
QT5=0.05
QTV1_F0=np.percentile(et, QT1 * 100)
QTV5_F0=np.percentile(et, QT5 * 100)
QTV1_F1=np.percentile(et1, QT1 * 100)
QTV5_F1=np.percentile(et1, QT5 * 100)
QTV1_F2=np.percentile(et2, QT1 * 100)
QTV5_F2=np.percentile(et2, QT5 * 100)

# time rescale
t=[datetime.datetime.strptime(d,'%Y%m%d').strftime('%Y-%m-%d') for d in dtime]
t0=pd.to_datetime(t)
t=[datetime.datetime.strptime(d,'%Y%m%d').strftime('%Y-%m-%d') for d in dtime1]
t1=pd.to_datetime(t)
t=[datetime.datetime.strptime(d,'%Y%m%d').strftime('%Y-%m-%d') for d in dtime2]
t2=pd.to_datetime(t)

# plot to one figure

plt.tight_layout()
fig2 = plt.figure(constrained_layout=True,figsize=(16,12))
gs = gridspec.GridSpec(nrows=3, ncols=2, width_ratios=[2, 1],height_ratios=[1,1,1], figure=fig2) 


f2_ax0= fig2.add_subplot(gs[0,0])
f2_ax0.set_xlim(np.datetime64('2006-01-01'), np.datetime64('2022-12-31'))
f2_ax0.plot(np.array(t0),np.array(et), '-',c='orange',label="Freq "+f0+" Hz")
f2_ax0.plot([t0[0],t0[-1]],[QTV5_F0,QTV5_F0], '-',c='k',label="End time")
f2_ax0.plot([t0[0],t0[-1]],[np.min(bt),np.min(bt)], '-',c='green',label="Begin time")
f2_ax0.set_title('Wondow endtime of '+sta+' (quantile 0.05 of noise level)')
f2_ax0.set_ylabel('Window endtime (sec)')
f2_ax0.set_xlabel('Time')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,7]))
plt.xticks(rotation=50)


f2_ax1 = fig2.add_subplot(gs[0,1])
f2_ax1.hist(et,bins = 100, alpha = 0.35, color = 'orange',label='Endtime at '+f0, orientation="horizontal")

f2_ax2 = plt.subplot(gs[1,0])
f2_ax2.set_xlim(np.datetime64('2006-01-01'), np.datetime64('2022-12-31'))
f2_ax2.plot(np.array(t1),np.array(et1), '-',c='red', label="Freq "+f1+" Hz")
f2_ax2.plot([t1[0],t1[-1]],[QTV5_F1,QTV5_F1], '-',c='k',label="End time")
f2_ax2.plot([t1[0],t1[-1]],[np.min(bt1),np.min(bt1)], '-',c='green',label="Begin time")
f2_ax2.set_title('Wondow endtime of '+sta+' (quantile 0.05 of noise level)')
f2_ax2.set_ylabel('Window endtime (sec)')
f2_ax2.set_xlabel('Time')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,7]))
plt.xticks(rotation=50)


f2_ax3 = fig2.add_subplot(gs[1,1])
f2_ax3.hist(et1,bins = 100, alpha = 0.35, color = 'red',label='Endtime at '+f1, orientation="horizontal")

f2_ax4 = fig2.add_subplot(gs[2,0])
f2_ax4.set_xlim(np.datetime64('2006-01-01'), np.datetime64('2022-12-31'))
f2_ax4.plot(np.array(t2),np.array(et2), '-',c='b',label="Freq "+f2+" Hz")
f2_ax4.plot([t2[0],t2[-1]],[QTV5_F2,QTV5_F2], '-',c='k',label="End time")
f2_ax4.plot([t2[0],t2[-1]],[np.min(bt2),np.min(bt2)], '-',c='green',label="Begin time")
f2_ax4.set_title('Wondow endtime of '+sta+' (quantile 0.05 of noise level)')
f2_ax4.set_ylabel('Window endtime (sec)')
f2_ax4.set_xlabel('Time')
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
plt.gca().xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,7]))
plt.xticks(rotation=50)


f2_ax5 = fig2.add_subplot(gs[2,1])
f2_ax5.hist(et2,bins = 100, alpha = 0.25, color = 'blue',label='Endtime at '+f2, orientation="horizontal")

f2_ax1.axhline(np.min(bt),   color='green' , linestyle='solid', linewidth=1,label='Begin Time '+f0)
f2_ax3.axhline(np.min(bt1),  color='green' , linestyle='solid', linewidth=1,label='Begin Time '+f1)
f2_ax5.axhline(np.min(bt2),  color='green' , linestyle='solid', linewidth=1,label='Begin Time '+f2)
f2_ax1.axhline(QTV5_F0,   color='k' , linestyle='solid', linewidth=1,label='quantile 0.05 of EndTime '+f0)
f2_ax3.axhline(QTV5_F1,   color='k' , linestyle='solid', linewidth=1,label='quantile 0.05 of EndTime '+f1)
f2_ax5.axhline(QTV5_F2,   color='k' , linestyle='solid', linewidth=1,label='quantile 0.05 of EndTime '+f2)

f2_ax0.grid(True)
f2_ax2.grid(True)
f2_ax4.grid(True)
f2_ax0.legend(loc='upper right')
f2_ax1.legend(loc='upper right')
f2_ax2.legend(loc='upper right')
f2_ax3.legend(loc='upper right')
f2_ax4.legend(loc='upper right')
f2_ax5.legend(loc='upper right')

plt.title('Window endtime of '+sta+' (quantile 0.05 of noise level)')
plt.ylabel('Endtime from quantile 0.05 of noise level')
plt.xlabel('Count')

plt.savefig('WindowEndtime_'+sta+'.png')


# --- print results to file
ftxt="ENDTIME_QT5.txt"
f=open(ftxt, "a")
f.write('stan,tb('+f0+')'+',te('+f0+')'+',tb('+f1+')'+',te('+f1+')'+',tb('+f2+')'+',te('+f2+')\n')
line=sta+",%.6f"%np.min(bt)+",%.6f"%QTV5_F0+",%.6f"%np.min(bt1)+",%.6f"%QTV5_F1+",%.6f"%np.min(bt2)+",%.6f"%QTV5_F2+"\n"
f.write(str(line))
f.close()
