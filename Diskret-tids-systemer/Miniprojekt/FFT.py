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
N = 2**8
f_s = 2**8
t = N/float(f_s)
f_res = f_s/float(N)
x = np.linspace(0,t,N)

def f(x):
    return np.sin(np.pi/3*x)

def g(x):
    return np.sin(2*np.pi/3 + np.pi/2*x)

def h(x):
    return np.sin(4*np.pi/3 + 3*np.pi/4*x)

def j(x):
    return np.sin(x)

y = f(x)+g(x)+h(x)

#==============================================================================
# DFT
#==============================================================================
def DFT(x,c):
    X = np.zeros(c,dtype=complex)
    for k in range(len(x)):
        a = 0+0*1j
        for n in range(c):
            a += x[n]*np.exp(-2*np.pi*1j*k*n/float(c))
            X[k] = a/np.sqrt(N)
    return X

#==============================================================================
# FFT
#==============================================================================
    
def FFT(x):
    """A recursive implementation of the 1D Cooley-Tukey FFT"""
    N_new = len(x)
    if N % 2 > 0:
        raise ValueError('Nej')
    elif N_new <= 1:  # this cutoff should be optimized
        return DFT(x,N_new)
    else:
        X_even = FFT(x[::2])
        X_odd = FFT(x[1::2])
        factor = np.exp(-2j * np.pi * np.arange(N_new) / N_new)
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
print 'Seconds to evaluate FFT', FFT_time



#plt.plot(x,y)
plt.plot(Y)
#plt.plot(Y_slow)
#plt.plot(np.fft.fft(y))

#print(np.allclose(DFT(y,N),np.fft.fft(y)))




