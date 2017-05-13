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
from scipy import signal


 
#windows

def w1(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if np.abs(n[i]) <= M:     
            w[i]=1
    return w

def w2(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if np.abs(n[i]) <= M/2:
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
M=100
N=M+1

sampels = 10000
n = np.linspace(0,M,N) 
freq_ax = np.linspace(0,np.pi,sampels/2)
w = w1(n,M)




plt.plot(n,w1(n,M), label = "Rectangular")
plt.plot(n,w2(n,M), label = "Bartlett")
plt.plot(n,w3(n,M), label = "Hanning")
plt.plot(n,w4(n,M), label = "Hamming")
plt.plot(n,w5(n,M), label = "Blackman")
plt.axis([0,M,0,1.2])
plt.legend(loc= "upper right", bbox_to_anchor=(1., 1))
plt.xlabel("n")
plt.ylabel("w [n]")
plt.savefig("figures/window_types.pdf")
plt.show()

#W = np.abs(np.fft.fft(np.pad(w,(0,sampels-N),'constant',constant_values=0)))
#W_dB= 20 * np.log10(np.abs(W))
#
#plt.plot(freq_ax,W[:len(freq_ax)])
#plt.xlabel("Frequency [rad/s.]")
#plt.ylabel("|W($\omega$)|")
##plt.savefig("figures/W_rect.pdf")
#plt.show()
#plt.plot(freq_ax,W_dB[:len(freq_ax)])
#plt.show()


#plt.gca().xaxis.set_major_locator(plt.NullLocator())
#plt.gca().yaxis.set_major_locator(plt.NullLocator())