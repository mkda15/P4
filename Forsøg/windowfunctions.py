# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 08:46:36 2017

@author: cht15
"""
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

#==============================================================================
# Window function, for implementing filters or the STFT.
#==============================================================================
""" Rectangular vindue """
def Rectangular(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] >= 0 and n[i] <= M:     
            w[i]=1
        else:
            w[i] = 0
    return w

""" Bartlett vindue """
def Bartlett(n, M):                                         # Bartlett vindue
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] <= M/2 and n[i] >= 0:
            w[i] = ((2*n[i])/M)
        elif n[i] > M/2 and n[i] <= M:
            w[i] = 2-(2*n[i]/M)
        else:
            w[i] = 0
    return w
  
""" Hanning vindue """    
def Hanning(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] >= 0 and n[i] <= M:
            w[i]=0.5-0.5*np.cos((2*np.pi*n[i])/M)
        else:
            w[i] = 0
    return w

""" Hamming vindue """ 
def Hamming(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] >= 0 and n[i] <= M:
            w[i]=0.54-0.46*np.cos((2*np.pi*n[i])/M)
        else:
            w[i] = 0
    return w

""" Blackman vindue """
def Blackman(n,M):
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] >= 0 and n[i] <= M:
            w[i]=0.42-0.5*np.cos((2*np.pi*n[i])/M)+0.08*np.cos(4*np.pi*n[i]/M)
        else:
            w[i] = 0
    return w

""" Kaiser vindue"""
def Kaiser( d1, d2, fs):                                        # Kaiservindue
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
        
    tw = ((2*d2)/fs)*2*np.pi
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
    return w,M,beta,A,n

#plt.style.use('classic')
#
#M=50
#n = np.linspace(0,M+1,M+2)
#
#w = Rectangular(n,M)
#plt.plot(n,w)
#plt.show()
#
#W = np.fft.fft(w,len(n))
#plt.plot(n,W)