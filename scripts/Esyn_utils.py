import obspy
import numpy as np
import pyasdf
import scipy
import math

from obspy.signal.filter import bandpass

### -----
def read_pyasdf(sfile,ccomp):
    # useful parameters from each asdf file
    with pyasdf.ASDFDataSet(sfile, mode="r") as ds:
        alist = ds.auxiliary_data.list()
        try:
            dt = ds.auxiliary_data[alist[0]][ccomp].parameters["dt"]
            dist = ds.auxiliary_data[alist[0]][ccomp].parameters["dist"]
            print("working on %s (comp: %s) that is %5.2fkm apart. dt: %.3f " % (sfile, ccomp, dist, dt))
            # read stacked data 
            sdata = ds.auxiliary_data[alist[0]][ccomp].data[:]
            para = ds.auxiliary_data[alist[0]][ccomp].parameters

            # time domain variables
            nwin = len(alist[1:])
            npts = sdata.size
            tvec = np.arange(-npts // 2 + 1, npts // 2 + 1) * dt
            #indx = np.where(np.abs(tvec) <= lag)[0]
            return dist,dt,tvec,sdata
        
        except Exception:
            print("continue! no %s component exist" % ccomp)
            return None
        
### -----
#Function that Calculate Mean Square
def msValue(arr, n):
    square = 0
    mean = 0.0

    #Calculate square
    for i in range(0,n):
        square += (arr[i]**2)

    #Calculate Mean
    mean = square / (float)(n)

    return mean

### -----
# Dirac delta function
def impulse(x):
    return 1 * (x == 0)

# Heaviside function (step function)
def step(x):
    return 1 * (x > 0)

# Esyn -->The 2-D radiative transfer equation for scalar waves
def ESYN_RadiaTrans_onesta(mean_free,tm, r, c):

    s0=c**2 * tm**2 -r**2
    if s0 > 0:
        # second term
        ind2=mean_free**(-1) * (math.sqrt(s0)-c*tm)
        a2up=math.exp(ind2)
        a2bot=2 * math.pi * mean_free * math.sqrt(s0)
        second=(a2up/a2bot)* step(tm-r/c)
        Esyn=  second

        return Esyn

# Esyn -->The 2-D radiative transfer equation for scalar waves
def ESYN_RadiaTrans_intersta(mean_free,tm, r, c):

    s0=c**2 * tm**2 -r**2
    if s0 > 0:
        # first term
        a1up=math.exp(-1 * c * tm * (mean_free**(-1)))
        a1bot= 2 * math.pi * c * r
        first= (a1up/a1bot)* impulse( tm-r/c )

        # second term
        ind2=mean_free**(-1) * (math.sqrt(s0)-c*tm)
        a2up=math.exp(ind2)
        a2bot=2 * math.pi * mean_free * math.sqrt(s0)
        second=(a2up/a2bot)* step(tm-r/c)

        Esyn= first + second

        return Esyn

### -----

### -----
### -----
### -----
### -----

