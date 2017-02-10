# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 09:08:05 2017

@author: Frederik Vardinghus
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

j = np.complex(0,1)
N = 500

x = np.zeros(N)
#x[3] = 1
lin = np.linspace(0,2*np.pi,N)
#lin = np.linspace(-np.pi,np.pi,N)
for i in range(len(lin)):
    x[i] = np.sin(lin[i])

def X(x,f):
    a = 0
    for n in range(len(x)):
        a += x[n]*np.exp(-j*2*np.pi*f*n)
    return a

#lin2 = np.linspace(-0.5,0.5,N)
lin2 = np.linspace(0,1,N)

amp = np.zeros(len(lin),dtype="complex64")
for i in range(len(amp)):
    amp[i] = X(x,lin2[i])
    
plt.plot(lin,np.real(amp),".")
