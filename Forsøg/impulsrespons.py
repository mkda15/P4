# -*- coding: utf-8 -*-
"""
Created on Mon May 01 13:22:16 2017

@author: cht15
"""

import numpy as np
#==============================================================================
# 
#==============================================================================
""" Impulsrespons Band Stop """
def ImpulsresponsBS(n,M,cut1,cut2): # Den ønskværdige impulsrespons for et båndstop filter
    hd = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] == M/2:
            hd[i] = 1 - (cut2 - cut1) * 2
        else:
            hd[i] = (np.sin(2 * np.pi * cut1 * (n[i] - M/2.)) / (np.pi*(n[i] - M/2.))) \
            - (np.sin(2 * np.pi * cut2 * (n[i] - M/2.)) / (np.pi*(n[i] - M/2.)))
    return hd

""" Impulsrespons High Pass """
def ImpulsresponsHP(n, M, f):     # Den ønskværdige impulsrespons for et højpass filter
    hd = np.zeros(len(n))
    for i in range(len(hd)):
        if n[i] == M/2:
            hd[i] = 1-f * 2
        else:
            hd[i] = (-np.sin(2 * np.pi * f * (n[i] - (M/2))))/(np.pi * (n[i] - (M/2)))
    return hd       

""" Impulserespons Low Pass """
def ImpulsresponsLP(n, M, f):     # Den ønskværdige impulsrespons for et lavpass filter
    hd = np.zeros(len(n))
    for i in range(len(hd)):
        if n[i] == M/2:
            hd[i] = f * 2
        else:
            hd[i] = (np.sin(2 * np.pi * f*(n[i] - (M/2))))/(np.pi * (n[i] - (M/2)))
    return hd 
""" Impulserespons Band pass """
def ImpulsresponsBP(n,M,ft1,ft2): # Den ønskværdige impulsrespons for et bandpass filter
    h = np.zeros(len(n))
    for i in range(len(n)):
        if i == M/2:
            h[i] = 2*(ft2 - ft1)
        else:
            h[i] = (1 / (np.pi*(i - M/2.)))*(np.sin(ft2*2*np.pi*(i - M/2.)) \
            - (np.sin(ft1*2*np.pi*(i - M/2.))))
    return h

#==============================================================================
# Noise adding
#==============================================================================
""" Noise adding """
def add_noise(data,noise,c = 0.5): #kilde side 229 i DTSP, tager to signaler og addere de to 
#    signal = np.convolve(data,noise)
#    return signal
    signal=np.zeros(len(data))    
    for i in range(len(data)):
        signal[i]=data[i]+float(c)*noise[i]
    return signal