# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 09:40:07 2017

@author: Frederik Vardinghus
"""

import matplotlib.pyplot as plt
import numpy as np
import time
from scipy import signal

#==============================================================================
# Generer data
#==============================================================================
N = 2**10 # Antal samples og laengde af FFT
f_s = 2**5   # Samplingsfrekvens
td = 1/float(f_s) # Samplingsperiode
t = td*N

x = np.linspace(0,N*td,N) # Samplingspunkter i tid
xf = 2*np.pi*np.linspace(0,1/float(2*td),N/float(2)) # Hoejre halvdel af samplingspunkter i frekvens

def f(x):
    return np.sin(np.pi/3*x)
def g(x):
    return np.sin(2*np.pi/3 + np.pi/2*x)
def h(x):
    return np.sin(4*np.pi/3 + 3*np.pi/4*x)
def j(x):
    return np.sin(x)+np.sin(3*x)

def window(x): # Hammingvindue
    return 0.54-0.46*np.cos((2*np.pi*x)/float(N-1))

w = window(np.linspace(0,N-1,N))

y = f(x)+g(x)+h(x)# Funktion, som samples og transformeres
#y = y*w # Windowing af signal
     
#==============================================================================
# DFT
#==============================================================================

def DFT(x,c):
    X = np.zeros(c,dtype=complex)
    for k in range(len(x)):
        a = 0+0*1j
        for n in range(c):
            a += x[n]*np.exp(-2*np.pi*1j*k*n/float(c))
            X[k] = a
    return X

#==============================================================================
# FFT
#==============================================================================

def is_power2(num): # Checks if a number is a power of 2
	return num != 0 and ((num & (num - 1)) == 0)

def FFT(x):
    N_new = len(x)
    if is_power2(N) == False:
        raise ValueError('N should be a power of 2.') # Brug N = potenser af 2
    elif N_new == 2:
        return DFT(x,N_new) # Returnerer DFT naar data ikke kan deles mere op
    else:
        X_even = FFT(x[::2]) # Deler rekursivt input op - lige dele
        X_odd = FFT(x[1::2]) # Deler rekursivt input op - ulige dele
        factor = np.exp(-2j * np.pi * np.arange(N_new) / N_new) # Twiddlefaktor
        return np.concatenate([X_even + factor[:N_new / 2] * X_odd,
                               X_even + factor[N_new / 2:] * X_odd])

#==============================================================================
# Egen DFT og hastighed
#==============================================================================
#start = time.time()
#Y_slow = DFT(y,N)
#end = time.time()
#DFT_time = end - start
#print'Seconds to evaluate DFT', DFT_time

#==============================================================================
# Egen FFT og hastighed
#==============================================================================
start = time.time()
Y = FFT(y)
end = time.time()
FFT_time = end - start
Y = 2/float(N)*np.abs(Y[:N/2])

#==============================================================================
# numpy.fft og hastighed
#==============================================================================
start2 = time.time()
Y2 = np.fft.fft(y)
end2 = time.time()
FFT2_time = end2 - start2
Y2 = 2/float(N)*np.abs(Y2[:N/2])

#==============================================================================
# Plots
#==============================================================================
#plt.style.use('ggplot')
plt.plot(x[:10],y[:10])
#plt.plot(xf,Y)
#plt.plot(xf,Y2)
plt.legend('Y')
plt.xlabel('Angular frequency')
plt.ylabel('Amplitude')
plt.show()

Y_sort = np.argpartition(Y2,-3)[-3:]

print 'Seconds to evaluate own FFT:  ', FFT_time
print 'Seconds to evaluate numpy.fft:', FFT2_time