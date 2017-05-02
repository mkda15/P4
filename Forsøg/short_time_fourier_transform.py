# -*- coding: utf-8 -*-
"""
Created on Tue May 02 08:50:23 2017

@author: cht15
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sc

#==============================================================================
# STFT
#==============================================================================
def stft(x, fftsize = 101, overlap = 2):   
    hop = fftsize / overlap
    w = sc.hanning(fftsize+1)[:-1]      # better reconstruction with this trick +1)[:-1]  
    return np.array([np.fft.rfft(w*x[i:i+fftsize]) for i in range(0, len(x)-fftsize, hop)])

#==============================================================================
# DFT
#==============================================================================
def dft(x,c):
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
def fft(x):
    N_new = len(x)
    if N_new == 2:
        return dft(x) # Returnerer DFT naar data ikke kan deles mere op
    else:
        X_even = fft(x[::2]) # Deler rekursivt input op - lige dele
        X_odd = fft(x[1::2]) # Deler rekursivt input op - ulige dele
        factor = np.exp(-2j * np.pi * np.arange(N_new) / N_new) # Twiddlefaktor
        return np.concatenate([X_even + factor[:N_new / 2] * X_odd,
                               X_even + factor[N_new / 2:] * X_odd])

#==============================================================================
# Decibel calculation
#==============================================================================
def db(x):
    return 20*np.log10(x)