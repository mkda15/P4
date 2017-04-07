# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 14:34:11 2017

@author: Trine
"""
#==============================================================================
# Plot of windos for report 
#==============================================================================

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


w_range = np.linspace(-np.pi, np.pi, 1000)
n = np.linspace(0,11,len(w_range))
 
#windows

def w1(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] <= M:     
            w[i]=1
    return w

def w2(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] <= M/2:
            w[i]=(2*n[i])/M
        if n[i] > M/2 and n[i] <= M:
            w[i]=2-(2*n[i]/M)
    return w
    
def w3(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] <= M:
            w[i]=0.5-0.5*np.cos((2*np.pi*n[i])/M)
    return w
 
def w4(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] <= M:
            w[i]=0.54-0.46*np.cos((2*np.pi*n[i])/M)
    return w

def w5(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] <= M:
            w[i]=0.42-0.5*np.cos((2*np.pi*n[i])/M)+0.08*np.cos(4*np.pi*n[i]/M)
    return w
   
w1=w1(n,8)
w2=w2(n,8)
w3=w3(n,8)
w4=w4(n,8)
w5=w5(n,8)


W_fft=np.fft.fft(w1)

dB = 20*np.log(np.abs(W_fft))

plt.plot(n,w1, label = "Rectangular")
plt.plot(n,w2, label = "Bartlett")
plt.plot(n,w3, label = "Hanning")
plt.plot(n,w4, label = "Hamming")
plt.plot(n,w5, label = "Blackman")
plt.axis([0,11,0,1.2])
plt.legend(loc= "upper right", bbox_to_anchor=(1.12, 1))
plt.xlabel("n")
plt.ylabel("w [n]")
#plt.savefig("figures/window_types.pdf")
plt.show()