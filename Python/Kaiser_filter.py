# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 15:37:25 2017

@author: Martin Kamp Dalgaard
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
 
fs = 44100.
if not type(fs) == float:
    raise ValueError("The sampling frequency should be a float.")
ft = 50 / fs
ft1 = 500 / fs
ft2 = 1000 / fs
b = 0.1
d = 0.5
A = int(-20*np.log10(d))

# N = 201
N = int(np.ceil((A - 8) / (2.285 * 2 * np.pi * b))) + 1 # Length of the filter
if not N % 2: N += 1  # Make sure that N is odd.
M = N-1 # Order of the filter
M1 = 200

if A > 50:
    beta = 0.1102 * (A - 8.7)
elif A <= 50 and A >= 21:
    beta = 0.5842 * (A - 21) ** 0.4 + 0.07886 * (A - 21)
else:
    beta = 0

#alpha = 3 # Notice: the Kaiser window is a rectangular window for alpha = 0.
#beta = np.pi*alpha
n = np.arange(N)


#==============================================================================
# Windows
#==============================================================================

def Kaiser(n,M): # Kaiser window
    w = np.zeros(len(n))
    if beta == 0 :
        for i in range(len(n)):
            if n[i] >= 0 and n[i] <= M:
                w[i] = 1
            else:
                w[i] = 0
        return w
    else:
        for i in range(len(n)):
            sum_t = 0
            sum_n = 0
            for j in range(M1):
                sum_t += ((1/np.math.factorial(j))**2) * (((beta/2)*np.sqrt(1 - ((2*i)/(N-1) - 1)**2))**(2*j))
                sum_n += ((1/np.math.factorial(j))**2) * ((beta/2)**(2*j))
            w[i] = sum_t/sum_n
    return w

def ha(n,M,a): # Hann window, if a = 0.5. Hamming window, if a = 0.54.
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] >= 0 and n[i] <= M:
            w[i] = a - (1 - a)*np.cos((2*np.pi*n[i])/M)
        else:
            w[i] = 0
    return w

def rect(n,M): # Rectangular window
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] >= 0 and n[i] <= M:
            w[i] = 1
        else:
            w[i] = 0
    return w

#==============================================================================
# Filters
#==============================================================================

def lp(n,M,ft): # Lowpass filter
    h = np.zeros(len(n))
    for i in range(len(n)):
        if i == M/2.:
            h[i] = 2*ft
        else:
            h[i] = np.sin(2*np.pi*ft*(i - M/2.)) / (np.pi*(i - M/2.))
    return h

def bp(n,M,ft1,ft2): # Bandpass filter
    h = np.zeros(len(n))
    for i in range(len(n)):
        if i == M/2:
            h[i] = 2*(ft2 - ft1)
        else:
            h[i] = (1 / (np.pi*(i - M/2.)))*(np.sin(ft2*2*np.pi*(i - M/2.)) \
            - (np.sin(ft1*2*np.pi*(i - M/2.))))
    return h

def hp(n,M,ft): # Highpass filter
    h = np.zeros(len(n))
    for i in range(len(n)):
        if i == M/2.:            
            h[i] = 1 - 2*ft
        else:
            h[i] = - np.sin(2*np.pi*ft*(i - M/2.)) / (np.pi*(i - M/2.))            
    return h

h = bp(n,M,ft1,ft2)*Kaiser(n,M1)

# Normalize to get unity gain.
h = h / np.sum(h)

omega = np.linspace(0,np.pi,len(n))

H = np.fft.fft(h)

plt.plot(n,(np.abs(H)))