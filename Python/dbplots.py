# -*- coding: utf-8 -*-
"""
Created on Fri May 19 11:40:10 2017

@author: Frederik Vardinghus
"""




import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

N = 2**20
M = 100
p = N-M
beta = 6
yaxis_max = 1.2
bins = 2*np.pi*np.linspace(0,0.5,N/2)

def rect(M):
    return np.array([1 for i in range(M)])

def fft(x):
    return np.abs(np.fft.fft(x))/np.max(np.abs(np.fft.fft(x)))

def db(x):
    return 20*np.log10(x)

def pad(x,p):
    return np.pad(x,(0,p),'constant',constant_values=0)

def k(M,beta):
    return pad(sc.kaiser(M,beta),p)

hamming = pad(sc.hamming(M),p)
hann = pad(sc.hanning(M),p)
bartlett= pad(sc.bartlett(M),p)
blackman = pad(sc.blackman(M),p)
kaiser = k(M,beta)
rectangular = pad(rect(M),p)

i = 0
while bins[i] < yaxis_max:
    inter = i
    i += 1

plt.plot(bins[:inter],db(fft(rectangular)[:inter]))
plt.xlabel('Frequency [rad/s]',fontsize=13)
plt.ylabel('Amplitude [dB]',fontsize=13)
plt.axis([0,bins[inter-1],-100,0])















