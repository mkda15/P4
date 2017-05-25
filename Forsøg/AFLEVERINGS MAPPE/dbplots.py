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
xaxis_max = 1.2
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
while bins[i] < xaxis_max:
    inter = i
    i += 1

# dB plot of window as specified above

plt.figure(1)
plt.plot(bins[:inter],db(fft(k(M,beta))[:inter]))
plt.xlabel('Frequency [rad/s]',fontsize=13)
plt.ylabel('Amplitude [dB]',fontsize=13)
plt.axis([0,bins[inter-1],-100,0])
plt.legend(loc = "best")
plt.savefig("kaiser.pdf")

# dB plots of the Kaiser window with different values of beta and M

M = [50, 100, 150]
beta = [0, 2, 4]
color = ["b-", "r-", "g-"]

for i in range(len(M)):
    plt.figure(2)
    plt.plot(bins[:inter],db(fft(k(100,beta[i]))[:inter]), color[i], label = r"$\beta \ = \ %d$" %(beta[i]))
    plt.xlabel('Frequency [rad/s]',fontsize=13)
    plt.ylabel('Amplitude [dB]',fontsize=13)
    plt.axis([0,bins[inter-1],-100,0])
    plt.legend(loc = "best")
    plt.savefig("kaiser_beta.pdf")

for i in range(len(M)):
    plt.figure(3)
    plt.plot(bins[:inter],db(fft(k(M[i],6))[:inter]), color[i], label = r"$M \ = \ %d$" %(M[i]))
    plt.xlabel('Frequency [rad/s]',fontsize=13)
    plt.ylabel('Amplitude [dB]',fontsize=13)
    plt.axis([0,bins[inter-1],-100,0])
    plt.legend(loc = "best")
    plt.savefig("kaiser_order.pdf")