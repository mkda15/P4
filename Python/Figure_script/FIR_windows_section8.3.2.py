# -*- coding: utf-8 -*-
"""
Created on Tue May 23 22:12:33 2017

@author: Martin Kamp Dalgaard
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def w1(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if np.abs(n[i]) <= M and n[i] >= 0:  
            w[i] = 1
        else:
            w[i] = 0
    return w

def w2(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if np.abs(n[i]) <= M/2 and n[i] >= 0:
            w[i]=(2*n[i])/M
        elif n[i] > M/2 and n[i] <= M:
            w[i] = 2 - (2*n[i]/M)
        else:
            w[i] = 0
    return w
    
def w3(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] <= M and n[i] >= 0:
            w[i] = 0.5 - 0.5*np.cos((2*np.pi*n[i])/M)
        else:
            w[i] = 0
    return w
 
def w4(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] <= M and n[i] >= 0:
            w[i] = 0.54 - 0.46*np.cos((2*np.pi*n[i])/M)
        else:
            w[i] = 0
    return w

def w5(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] <= M and n[i] >= 0:
            w[i] = 0.42 - 0.5*np.cos((2*np.pi*n[i])/M) + 0.08*np.cos(4*np.pi*n[i]/M)
        else:
            w[i] = 0
    return w

M = 100
N = M+1

samples = 10000
n = np.linspace(-10,N+10,samples)

plt.plot(n,w1(n,M), "b-", label = "Rectangular")
plt.plot(n,w2(n,M), "g-", label = "Bartlett")
plt.plot(n,w3(n,M), "r-", label = "Hann")
plt.plot(n,w4(n,M), "c-", label = "Hamming")
plt.plot(n,w5(n,M), "m-", label = "Blackman")
plt.axis([-M/10,M+M/10,-0.2,1.2])
plt.legend(loc= "best")
plt.xlabel(r"$n$")
plt.ylabel(r"$w[n]$")
plt.savefig("figures/window_types.pdf")
plt.show()