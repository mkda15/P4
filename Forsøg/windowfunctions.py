# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 08:46:36 2017

@author: cht15
"""
import numpy as np
import scipy as sc

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

""" BArtlett vindue """
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
def Kaiser(n, M, B):                                        # Kaiservindue
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] >= 0 and n[i] <= M:
            variabel = B*np.sqrt(1-((n[i]-(M/2))/(M/2))**2)
            w[i] = sc.special.i0(variabel)/sc.special.i0(B)
        else:
            w[i] = 0
    return w