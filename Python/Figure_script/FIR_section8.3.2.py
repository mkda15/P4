# -*- coding: utf-8 -*-
"""
Created on Fri May 19 14:40:53 2017

@author: Martin Kamp Dalgaard
"""

import numpy as np
import matplotlib.pyplot as plt

def rectangular(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if np.abs(n[i]) <= M:     
            w[i] = 1
        else:
            w[i] = 0
    return w

def hd(n,M,f):
    hd = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] == M/2:
            hd[i] = f/np.pi
        else:
            hd[i] = (np.sin(f*(n[i] - M/2)))/(np.pi*(n[i]-M/2))
    return hd

M = 50
l = M+1

n = np.linspace(0,l,l+1)
x = np.linspace(0,2*np.pi,len(n))

lp = hd(n,M,np.pi/2)

h = lp*rectangular(n,M)
H = np.abs(np.fft.fft(h))

plt.figure(2)
plt.plot(x,H)
plt.axis([0,2*np.pi,0,1.2])
plt.xlabel("Frequency [rad/s]")
plt.ylabel("Amplitude")
plt.savefig("lowpass_real.pdf")