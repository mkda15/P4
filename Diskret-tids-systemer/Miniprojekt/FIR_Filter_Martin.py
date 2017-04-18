# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 21:50:53 2017

@author: Martin Kamp Dalgaard
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

M = 200
l = M+1
Fs = 0.9

n = np.linspace(0,l,l+1)
x = np.linspace(-np.pi,np.pi,len(n))
xf = np.linspace(0,10*np.pi,len(n))

f1 = (5*np.pi/12.) / (2*np.pi)
f2 = (5*np.pi/8.) / (2*np.pi)

def h(n,M,f1,f2):
    hd = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] == M/2:
            hd[i] = 1 - 2*(f2 - f1)
        else:
            hd[i] = (np.sin(2*np.pi*f1*(n[i] - M/2.)) / (np.pi*(n[i] - M/2.))) \
            - (np.sin(2*np.pi*f2*(n[i] - M/2.)) / (np.pi*(n[i] - M/2.)))
    return hd

def ha(n,M,a): # Hanning window if a = 0.5. Hamming window if a = 0.54.
    w = np.zeros(len(n))
    for i in range(len(n)):
        if abs(x[i]) <= M/2.:
            w[i] = a - (1 - a)*np.cos((2*np.pi*n[i])/M)
        else:
            w[i] = 0
    return w

def blackman(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        w[i] = 0.42 - 0.5*np.cos((2*np.pi*n[i])/M) + 0.8*np.cos((4*np.pi*n[i])/M)
    return w
    
w = ha(n,M,0.54)
hd = h(n,M,f1,f2)

def fft(x,n):
    return np.fft.fft(x)

h = hd * w

H = np.abs(fft(h,n))

plt.figure(2)
plt.plot(x, H)
plt.axis([0,np.pi,0,2])
plt.title('Our bandstop filter frequency response')
plt.xlabel('Frequency [radians / second]')
plt.ylabel('Amplitude')
plt.axvline(f1*(2*np.pi), color='yellow') # lower cutoff frequency
plt.axvline(f2*(2*np.pi), color='yellow') # upper cutoff frequency
plt.axvline(np.pi/2, color='red') # frequency to be eliminated
plt.axvline(np.pi/3, color='green') # frequency to keep
plt.axvline(3*np.pi/4, color='green') # frequency to keep

#==============================================================================
# Scipy 
#==============================================================================

omega1_scp = (5*np.pi/12)
omega2_scp = (5*np.pi/8)

N = [omega1_scp,omega2_scp]
plt.figure(3)
b, a = signal.butter(15, N, 'bandstop', analog=True)
w, h = signal.freqs(b, a)
plt.plot(w, abs(h), "b-", label = "Bandstop filter")
plt.title('Scipys bandstop filter frequency response')
plt.xlabel('Frequency [radians / second]')
plt.ylabel('Amplitude')
plt.legend(loc = "lower left")
plt.margins(0, 0.1)
plt.axis([0,np.pi,0,2])
plt.grid(which='both', axis='both')
plt.axvline(omega1_scp, color='yellow') # lower cutoff frequency
plt.axvline(omega2_scp, color='yellow') # upper cutoff frequency
plt.axvline(np.pi/2, color='red') # frequency to be eliminated
plt.axvline(np.pi/3, color='green') # frequency to keep
plt.axvline(3*np.pi/4, color='green') # frequency to keep
plt.show()