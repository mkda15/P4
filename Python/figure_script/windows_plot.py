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
import scipy as sc
from scipy import signal



 
#windows

def w1(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if np.abs(n[i]) <= M and n[i] >= 0:  
            w[i]=1
        else:
            w[i] = 0
    return w

def w2(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if np.abs(n[i]) <= M/2 and n[i] >= 0:
            w[i]=(2*n[i])/M
        elif n[i] > M/2 and n[i] <= M:
            w[i]=2-(2*n[i]/M)
        else:
            w[i] = 0
    return w
    
def w3(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] <= M and n[i] >= 0:
            w[i]=0.5-0.5*np.cos((2*np.pi*n[i])/M)
        else:
            w[i] = 0
    return w
 
def w4(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] <= M and n[i] >= 0:
            w[i]=0.54-0.46*np.cos((2*np.pi*n[i])/M)
        else:
            w[i] = 0
    return w

def w5(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] <= M and n[i] >= 0:
            w[i]=0.42-0.5*np.cos((2*np.pi*n[i])/M)+0.08*np.cos(4*np.pi*n[i]/M)
        else:
            w[i] = 0
    return w

def Kaiser( d1, tw):                                        # Kaiservindue
    M = 0
    beta = 0
    # defining Beta 
    A = int(-20*np.log10(d1))
    
    if A > 50:
        beta = 0.1102 * (A - 8.7)
    elif A <= 50 and A >= 21:
        beta = 0.5842 * (A - 21) ** 0.4 + 0.07886 * (A - 21)
    else:
        beta = 0
        
    
    M = int(np.ceil((A - 8) / (2.285 * tw)))
    print('The filter order calculated by Kaiser window: \n M = %.0f' %M)
    
    n = np.linspace(0,M,M+1) 
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] >= 0 and n[i] < M:
            variabel = beta*np.sqrt(1-((n[i]-(M/2))/(M/2))**2)
            w[i] = sc.special.i0(variabel)/sc.special.i0(beta)
        else:
            w[i] = 0
    if (M%2 == 1):
        M=M+1
    return w,M,n,beta

def ImpulsresponsLP(n, M, cut):     # Den ønskværdige impulsrespons for et lavpass filter
    hd = np.zeros(len(n))
    for i in range(len(hd)):
        if n[i] == M/2:
            hd[i] = cut/np.pi
        else:
            hd[i] = (np.sin(cut*(n[i] - (M/2))))/(np.pi * (n[i] - (M/2)))
    return hd 



M=100
N=M+1

sampels = 10000
n = np.linspace(0,N,sampels)
freq_ax = np.linspace(0,np.pi,sampels/2)
w,M,n,beta = Kaiser(0.01, 0.5*np.pi-0.4*np.pi)
#
#plt.plot(n,w1(n,M), "b-", label = "Rectangular")
#plt.plot(n,w2(n,M), "g-", label = "Bartlett")
#plt.plot(n,w3(n,M), "r-", label = "Hanning")
#plt.plot(n,w4(n,M), "c-", label = "Hamming")
#plt.plot(n,w5(n,M), "m-", label = "Blackman")
#plt.axis([-M/10,M+M/10,-0.2,1.2])
#plt.legend(loc= "best")
#plt.xlabel("$n$")
#plt.ylabel("$w[n]$")
#plt.savefig("figures/window_types.pdf")
#plt.show()
#W = np.abs(np.fft.fft(w,len(n)))
W = np.abs(np.fft.fft(np.pad(w,(0,sampels-M+1),'constant',constant_values=0)))

W_dB= 20 * np.log10(np.abs(W))

hd = ImpulsresponsLP(n,M,np.pi/2.)
h = hd*w

x = np.linspace(0,np.pi,M)
H = np.abs(np.fft.fft(h,len(n)))

plt.plot(x[:len(H)/2.],H[:len(H)/2.])
plt.xlabel("Frequency [rad/s]")
plt.ylabel("Amplitude")
plt.savefig("figures/kaiser_H.png")
plt.show()

#plt.plot(freq_ax,W[:len(freq_ax)])
#plt.xlabel("Frequency [rad/s.]")
#plt.ylabel("|W($\omega$)|")
##plt.savefig("figures/W_rect.pdf")
#plt.show()
#plt.plot(freq_ax,W_dB[:len(freq_ax)])
#plt.show()


#plt.gca().xaxis.set_major_locator(plt.NullLocator())
#plt.gca().yaxis.set_major_locator(plt.NullLocator())