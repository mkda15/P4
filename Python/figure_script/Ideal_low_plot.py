# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 11:53:02 2017

@author: Trine
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


w_range = np.linspace(-np.pi, np.pi, 8)
n = np.linspace(-10,10,(len(w_range)))
M = len(n)-1

# plot at ideal lavpas filter w_range(10000), og n(11)
def h(n,w_c):
    h = np.zeros(len(n))
    for i in range(len(h)):
        if n[i] == 0:
            h[i]= w_c/np.pi
        else:
            h[i]=(np.sin(n[i]*w_c))/((np.pi*n[i]))
    return h 

#n = range(11)
cutoff = np.pi/2 

def H(w):
    H = np.zeros(len(w))
    for i in range(len(w)):
        if np.abs(w[i])< np.pi/2:
            H[i] = 1
        else :
            H[i] = 0
    return H

def w1(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] <= M:     
            w[i]=1
    return w



plt.plot(freq,W)
plt.axis([-0.1,0.1,0,2.5])
plt.show()


#plt.plot(20*np.log10(np.abs(np.fft.fft(np.pad(w,(0,10000),'constant',constant_values=0))))[:400])
#plt.show()
plt.plot(np.abs(np.fft.fft(np.pad(w,(0,10000),'constant',constant_values=0)))[:400])
plt.show()

#plt.stem(n,h(n,cutoff))
#plt.xlabel("n", fontsize=15)
#plt.ylabel("h[n]", fontsize=15)
##plt.savefig("figures/ideal_low1.pdf")
#plt.show()

plt.plot(w_range,H(w_range))
plt.axis([-np.pi, np.pi, 0, 1.5])
plt.xlabel("$\omega$", fontsize=20)
plt.ylabel("|H($\omega$)|", fontsize=15)
#plt.savefig("figures/ideal_low2.pdf")
plt.show()

