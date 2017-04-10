# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

def response(n,length,omega1,omega2):
    hd = np.zeros(length)
    for i in range(length):
#        if i == 0:
#            h[i] = 1 - ((omega2 - omega1)/np.pi)
#        else:
            hd[i] = (1./(np.pi*n[i])) * (np.sin(omega1*n[i])-np.sin(omega2*n[i]))
    return hd
delta = 0.5
M = 12
length = 100
x = np.linspace(-np.pi, np.pi,length)
n = np.linspace(0,(M/2), length)
omega1 = (np.pi)/2. - delta
omega2 = (np.pi)/2. + delta


hd = response(n,length,omega1,omega2)

plt.stem(n,hd)
plt.show()
#==============================================================================
# Windows
#==============================================================================
trine = np.linspace(0,M,length)
def hamming(n,length,M):
    w = np.zeros(length)
    for i in range(length):
        if n[i] <= M:
            w[i]=0.54-0.46*np.cos((2*np.pi*n[i])/M)
    return w
w = hamming(trine,length,M)

plt.plot(trine,w)
plt.show()
h = w * hd
W = np.fft.fft(w)
plt.plot(x,h)

