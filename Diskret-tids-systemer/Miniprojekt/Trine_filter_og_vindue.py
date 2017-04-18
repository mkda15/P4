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
            hd[i] = (1./(np.pi*(n[i]-(M/2.)))) * (np.sin(omega2 * (n[i]-(M/2.))) - np.sin( omega1 *(n[i]-(M/2.))))
    return hd

    
M = 40

length = M
samples = 100+length



n = np.linspace(0,length,(samples)+1)
x = np.linspace(0, length, len(n))

omega1 =   5*np.pi/12. #(np.pi)/2. - delta 
omega2 =   5*np.pi/8. # (np.pi)/2. + delta


hd = impulse_response(n,M,omega1,omega2)



plt.stem(n,hd)
plt.axis([0,20,-1.5,1.5])
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
plt.axis([-2,22,-0.8,0.8])
plt.show()

#==============================================================================
# Fourier tansformation
#==============================================================================


def fft(x):
    return np.fft.fft(x)

def dB(X):
    return 20*np.log(np.abs(X))

H = fft(h)
W = fft(w)

plt.plot(x,np.abs(W))
plt.axis([0,1,-10,21])
plt.show()

plt.plot(x,dB(W))
plt.axis([0,1,-100,100])
plt.show()

plt.plot(x,H)

plt.show()

plt.plot(x,fft(hd))
plt.show()