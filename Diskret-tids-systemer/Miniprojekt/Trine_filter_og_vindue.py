# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

def impulse_response(n,M,omega1,omega2):
    hd = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] == M/2.:
            hd[i] = 1 - ((omega2 - omega1)/np.pi)
        else:
            hd[i] = (1./(2*np.pi*(n[i]-(M/2.)))) * (np.sin(2*np.pi*omega1*(n[i]-(M/2))) - np.sin(np.pi*omega1*(n[i]-(M/2))))
    return hd
    
delta = 0.5
M = 8
length = M+1
samples = 100

n = np.linspace(0,length,samples)
x = np.linspace(-np.pi, np.pi, len(n))

omega1 = (np.pi)/2. - delta
omega2 = (np.pi)/2. + delta


hd = impulse_response(n,M,omega1,omega2)

plt.stem(n,hd)
plt.show()
#==============================================================================
# Windows
#==============================================================================

def hamming(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] <= M:
            w[i]=0.54-0.46*np.cos((2*np.pi*n[i])/M)
    return w
    
w = hamming(n,M)

plt.plot(n,w)
plt.show()


h = w * hd
plt.stem(n,h)
plt.axis([-2,15,-0.6,0.8])
plt.show()

#==============================================================================
# Fourier tansformation
#==============================================================================

def fft(x,n):
    return np.fft.fft(x)/len(n)

H = fft(h,n)
W = fft(w,n)
plt.plot(x,np.abs(H))
plt.show()

