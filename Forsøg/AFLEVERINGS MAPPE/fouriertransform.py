# -*- coding: utf-8 -*-

import numpy as np

#==============================================================================
# STFT
#==============================================================================

def stft(signal, fftsize = 2**12, overlap = 2):   
    hop = fftsize / overlap
    w = np.blackman(fftsize+1)[:-1]      # Better reconstruction with this trick +1)[:-1] 
    x =  np.array([w*signal[i:i+fftsize] for i in range(0, len(signal)-fftsize, hop)])
    X =  np.array([np.fft.rfft(w*signal[i:i+fftsize]) for i in range(0, len(signal)-fftsize, hop)])
    return X,x,w

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
        return dft(x) # Returns DFT when the length of the data is 2
    else:
        X_even = fft(x[::2]) # Even indexes of the data
        X_odd = fft(x[1::2]) # Odd indexes of the data 
        factor = np.exp(-2j * np.pi * np.arange(N_new) / N_new) # Twiddle factor
        return np.concatenate([X_even + factor[:N_new / 2] * X_odd,
                               X_even + factor[N_new / 2:] * X_odd])

#==============================================================================
# Decibel calculation
#==============================================================================

def db(x):
    return 20*np.log10(x)