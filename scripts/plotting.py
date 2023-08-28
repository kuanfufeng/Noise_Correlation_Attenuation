import obspy
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt


### ----- 
def plot_waveforms(ncmp,wav,fname,comp_arr):
    #fig, ax = plt.subplots(1,ncmp, figsize=(16,2), sharex=False)
    fig, ax = plt.subplots(ncmp,1, figsize=(8,12), sharex=False)

    for n in range(ncmp):
        absy=max(wav[n][1], key=abs)
        ax[n].set_ylim(absy*-1,absy)
        ax[n].plot(wav[n][0],wav[n][1], linewidth=0.5)
        ax[n].set_xlabel("time [s]")
        ax[n].set_title(fname+" "+comp_arr[n])
    fig.tight_layout()
    #print("save figure as Waveform_readin_%s.png"%(fname))
    plt.savefig("Waveform_readin_%s.png"%(fname), format="png", dpi=100)
    plt.close(fig)

### ----- 
def plot_filtered_waveforms(freq,tt,wav,fname,ccomp):
    nfreq = len(freq) - 1
    fig, ax = plt.subplots(1,nfreq, figsize=(12,2), sharex=False)
    
    for fb in range(nfreq):
        fmin=freq[fb]
        fmax=freq[fb+1]
        absy=max(wav[fb], key=abs)
        #absx=max(tt, key=abs)
        #ax[fb].set_xlim(absx*-1,absx)
        ax[fb].set_ylim(absy*-1,absy)
        ax[fb].plot(tt,wav[fb], "k-", linewidth=0.2)
        #ax[fb].plot(wav_fold[0],wav_fold[fb+1], "b-", linewidth=1)
        ax[fb].set_xlabel("Time [s]")
        ax[fb].set_ylabel("Amplitude")
        ax[fb].set_title( "%s   %s   @%4.2f-%4.2f Hz" % ( fname,ccomp,fmin,fmax ) )
    fig.tight_layout()
    plt.savefig("Waveform_filtered_%s_%s_F%s-%s.png"%(fname,ccomp,fmin,fmax), format="png", dpi=100)
    plt.close(fig)


### ----- 
def plot_envelope(comp_arr,freq,msv,msv_mean,fname,vdist):
    nfreq = len(freq) - 1
    ncmp = len(comp_arr)
    
    fig, ax = plt.subplots(ncmp+1,nfreq, figsize=(16,10), sharex=False)   
    for n in  range(len(comp_arr)):
        for fb in range(nfreq):
            fmin=freq[fb]
            fmax=freq[fb+1]    
            ax[n,fb].plot(msv[n][0][:], msv[n][fb+1], "k-", linewidth=0.5)
            ax[n,fb].set_title("%s   %.2fkm  %s   @%4.2f-%4.2f Hz" % (fname,vdist,comp_arr[n],fmin,fmax))
            ax[n,fb].set_xlabel("Time [s]")
            ax[n,fb].set_ylabel("Amplitude")
            
    for fb in range(nfreq):
        fmin=freq[fb]
        fmax=freq[fb+1]
        ax[-1,fb].plot(msv_mean[0], msv_mean[fb+1], "b-", linewidth=1)
        ax[-1,fb].set_title(" Mean Squared Value %.2fkm  @%4.2f-%4.2f Hz" % (vdist,fmin,fmax))
        ax[-1,fb].set_xlabel("Time [s]")
        ax[-1,fb].set_ylabel("Amplitude")            
    plt.tight_layout() 
    plt.savefig("Waveform_envelope_%s_F%s-%s.png"%(fname,fmin,fmax), format="png", dpi=100)
    plt.close(fig)

### ----- 
def plot_fmsv_waveforms(freq,wav,fname,noise_level,twin):
    nfreq = len(freq) - 1
    fig, ax = plt.subplots(2,nfreq, figsize=(16,6), sharex=False)

    for fb in range(nfreq):
        fmin=freq[fb]
        fmax=freq[fb+1]
        absy=1 #max(wav[fb], key=abs)
        ax[0][fb].set_ylim(-0.1,absy)
        ax[0][fb].plot(wav[0],wav[fb+1], "k-", linewidth=0.5)
        ax[0][fb].set_xlabel("Time [s]")
        ax[0][fb].set_ylabel("Amplitude")
        ax[0][fb].set_title( "%s   @%4.2f-%4.2f Hz" % ( fname,fmin,fmax ) )

        ax[1][fb].plot([wav[0][0],wav[0][-1]],[noise_level[fb],noise_level[fb]],c='orange',marker='o',ls='--', linewidth=2)
        
        ax[1][fb].plot([twin[fb][0],twin[fb][0]],[-0.1,absy],c='orange',marker='o',ls='--', linewidth=2)
        ax[1][fb].plot([twin[fb][1],twin[fb][1]],[-0.1,absy],c='orange',marker='o',ls='--', linewidth=2)
        
        ax[1][fb].set_yscale('log', base=10)
        ax[1][fb].plot(wav[0],wav[fb+1], "k-", linewidth=0.5)
        ax[1][fb].set_xlabel("Time [s]")
        ax[1][fb].set_ylabel("Amplitude in log-scale")
        ax[1][fb].set_title( "%s   @%4.2f-%4.2f Hz" % ( fname,fmin,fmax ) )
    fig.tight_layout()
    plt.savefig("Waveform_fmsv_%s.png"%(fname), format="png", dpi=100)
    plt.close(fig)

### ----- 
def plot_fitting_curves(mean_free,intrinsic_b,tt,Eobs,Esyn,fname,dist,twin,fmin,fmax):
    numb=len(intrinsic_b)
    plt.figure(figsize=(8,2))
    for nb in range(numb):

        plt.yscale('log', base=10)
        pymin=np.min(Eobs[nb][:-2]/2)
        pymax=np.max(Eobs[nb][:-2]*2)
        plt.ylim( pymin , pymax )
        plt.plot( tt, Eobs[nb], "k-", linewidth=0.5)
        plt.plot( tt, Esyn[nb], "b-", linewidth=1)
        plt.plot([twin[0],twin[0],twin[-1],twin[-1],twin[0]],[pymin, pymax,pymax,pymin,pymin],"r", linewidth=2)
        #plt.plot([twin[-1],twin[-1]],[np.min(Eobs[nb][:-2]*2), np.max(Eobs[nb][:-2]/2)],"r", linewidth=1)

    plt.title("%s  %.2fkm   @%4.2f-%4.2f Hz, mean_free: %.2f  b: %.2f~%.2f"
            % ( fname,dist,fmin,fmax,mean_free,intrinsic_b[0],intrinsic_b[-1]))
    plt.xlabel("Time [s]")
    plt.ylabel("Energy density Amplitude")
    plt.tight_layout()
    plt.savefig("Fitting_fmsv_%s_F%s-%s_MFP%.2f.png"%(fname,fmin,fmax,mean_free), format="png", dpi=100)
    plt.close()

### ----- 
def plot_fitting_result(mean_free,intrinsic_b,tt,Eobs,Esyn,fname,dist,twin,fmin,fmax):
    plt.figure(figsize=(6,2))
    plt.yscale('log', base=10)
    pymax=np.max(Eobs[1:]*5)
    pymin=10**(-6)
    plt.ylim( pymin , pymax )
    plt.plot( tt, Eobs, "k-", linewidth=1)
    plt.plot( tt, Esyn, "b--", linewidth=1)
    plt.plot([twin[0],twin[0],twin[-1],twin[-1],twin[0]],[pymin, pymax,pymax,pymin,pymin],"r", linewidth=2)
    #plt.plot([twin[-1],twin[-1]],[np.min(Eobs[nb][:-2]*2), np.max(Eobs[nb][:-2]/2)],"r", linewidth=1)

    plt.title("%s  %.2fkm   @%4.2f-%4.2f Hz, mfp: %.2f  b: %.2f"
            % ( fname,dist,fmin,fmax,mean_free,intrinsic_b))
    plt.xlabel("Time [s]")
    plt.ylabel("Energy density Amp")
    plt.tight_layout()   
    plt.savefig("Final_fmsv_%s_F%s-%s.png"%(fname,fmin,fmax), format="png", dpi=100)

### ----- 
def plot_grid_searching(sta_pair,freq,SSR,x,y):
    nfreq=len(freq)-1

    fig, ax = plt.subplots(1,nfreq, figsize=(16,4), sharex=False)
    
    for fb in range(nfreq): 
        fmin=freq[fb]
        fmax=freq[fb+1]

        loc=np.where(SSR[fb].T == np.amin(SSR[fb].T))
        locx=list(zip(loc[0], loc[1]))
        ymin=y[loc[0]]
        xmin=x[loc[1]]
        print(" intrinsic_b %.2f " % ymin,"mean_free: %.2f " % xmin)

        grid = SSR[fb].T
        im=ax[fb].imshow(grid,extent=(x.min(), x.max(), y.max(), y.min()), aspect='auto',cmap = 'viridis_r',interpolation='spline16' )
        #im=ax[fb].imshow(grid,aspect='auto',cmap = 'viridis_r' )
        im.set_clim(1,1.3)    
        cb=plt.colorbar(im,extend='max')
        cb.set_label('SSR/SSR_min', rotation=90, labelpad=14)
        ax[fb].set_title("%s  SSR  @%4.2f-%4.2f Hz" % (sta_pair,fmin,fmax) )
        ax[fb].set_xlabel("mean free path")
        ax[fb].set_ylabel("intrinsic_b")
        ax[fb].invert_yaxis()
        ax[fb].plot(xmin,ymin,"+", markersize=20, color='red')
    plt.tight_layout() 
    plt.savefig("Final_SSR_%s.png"%(fname), format="png", dpi=100)
### ----- 
### ----- 
