import matplotlib.pyplot as plt
import math

def plot_waveforms(num_cmp,wav,fname,comp_arr):
    fig, ax = plt.subplots(1,num_cmp, figsize=(16,3), sharex=False)
    
    for n in range(ncmp):
        absy=max(wav[n][1], key=abs)
        ax[n].set_ylim(absy*-1,absy)
        ax[n].plot(wav[n][0],wav[n][1])
        ax[n].set_xlabel("time [s]")
        ax[n].set_title(fname+" "+comp_arr[n])
    fig.tight_layout()
    plt.show()


## plot filtered waveform and two-side average stack
def plot_filtered_waveforms(freq,tt,wav,wav_fold,fname):
    nfreq = len(freq) - 1
    fig, ax = plt.subplots(1,nfreq, figsize=(16,3), sharex=False)
    
    for fb in range(nfreq):
        fmin=freq[fb]
        fmax=freq[fb+1]
        absy=max(wav[fb], key=abs)
        #absx=max(tt, key=abs)
        #ax[fb].set_xlim(absx*-1,absx)
        ax[fb].set_ylim(absy*-1,absy)
        ax[fb].plot(tt,wav[fb], "k-", linewidth=0.1)
        ax[fb].plot(wav_fold[0],wav_fold[fb+1], "b-", linewidth=1)
        ax[fb].set_xlabel("Time [s]")
        ax[fb].set_ylabel("Amplitude")
        ax[fb].set_title( "%s   %s   @%4.2f-%4.2f Hz" % ( fname,ccomp,fmin,fmax ) )
    
    fig.tight_layout()
    plt.show()


## plot envelope and mean-squared envelope
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
        ax[-1,fb].plot(msv_mean[0], msv_mean[fb+1], "b-", linewidth=1)
        ax[-1,fb].set_title(" Mean Squared Value %.2fkm  @%4.2f-%4.2f Hz" % (vdist,fmin,fmax))
        ax[-1,fb].set_xlabel("Time [s]")
        ax[-1,fb].set_ylabel("Amplitude")            

    plt.tight_layout()   
    plt.show()       

## plot fitting curves
def plot_fitting_curves(mean_free,intrinsic_b,tt,Eobs,Esyn,fname):
    numb=len(intrinsic_b)
    plt.figure(figsize=(8,2))
    for nb in range(numb):

        plt.yscale('log', base=10)
        plt.xlim(0,120)
        plt.ylim( np.min(Eobs[nb][:-2]/2) , np.max(Eobs[nb][:-2]*2) )
        plt.plot( tt, Eobs[nb], "k-", linewidth=0.5)
        plt.plot( tt, Esyn[nb], "b-", linewidth=1)

    plt.title("%s     @%4.2f-%4.2f Hz, mean_free: %.2f  b: %.2f~%.2f"
            % ( fname,fmin,fmax,mean_free,y[0],y[-1]))
    plt.xlabel("Time [s]")
    plt.ylabel("Energy density Amplitude")

## plot grid-searching results

