# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 09:40:07 2017

@author: Frederik Vardinghus
"""

import matplotlib.pyplot as plt
import numpy as np
import time

#==============================================================================
# Generér data
#==============================================================================
N = 2**15 # Antal samples og længde af FFT
f_s = 2**2 # Samplingsfrekvens
td = 1/float(f_s) # Samplingsperiode

x = np.linspace(0,N*td,N) # Samplingspunkter i tid
xf = np.linspace(0,1/float(2*td),N/float(2)) # Halvdelen af samplingspunkter i frekvens

def f(x):
    return np.sin(np.pi/3*x)
def g(x):
    return np.sin(2*np.pi/3 + np.pi/2*x)
def h(x):
    return np.sin(4*np.pi/3 + 3*np.pi/4*x)
def j(x):
    return np.sin(x)

y = f(x)+g(x)+h(x) # Funktion, som samples og transformeres

#==============================================================================
# DFT
#==============================================================================
def DFT(x,c):
    X = np.zeros(c,dtype=complex)
    for k in range(len(x)):
        a = 0+0*1j
        for n in range(c):
            a += x[n]*np.exp(-2*np.pi*1j*k*n/float(c))
            X[k] = a/float(np.sqrt(N))
    return X

#==============================================================================
# FFT
#==============================================================================
def FFT(x):
    N_new = len(x)
    if N % 2 > 0:
        raise ValueError('nej.') # Brug N = potenser af 2
    elif N_new <= 1:
        return DFT(x,N_new) # Returnerer DFT når data ikke kan deles mere op
    else:
        X_even = FFT(x[::2]) # Deler rekursivt input op - lige dele
        X_odd = FFT(x[1::2]) # Deler rekursivt input op - ulige dele
        factor = np.exp(-2j * np.pi * np.arange(N_new) / N_new) # Twiddlefaktor
        return np.concatenate([X_even + factor[:N_new / 2] * X_odd,
                               X_even + factor[N_new / 2:] * X_odd])


    
    
#start = time.time()
#Y_slow = DFT(y,N)
#end = time.time()
#DFT_time = end - start
#print'Seconds to evaluate DFT', DFT_time

start = time.time()
Y = FFT(y)
end = time.time()
FFT_time = end - start



#plt.plot(x,y)
plt.plot(xf,2/float(N)*np.abs(Y[:N/2]))

print 'Seconds to evaluate FFT', FFT_time



