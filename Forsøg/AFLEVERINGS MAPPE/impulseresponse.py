# -*- coding: utf-8 -*-

import numpy as np

""" Impulse response for bandstop filter """
def ImpulseresponseBS(n,M,cut1,cut2):
    hd = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] == M/2:
            hd[i] = 1 - (cut2 - cut1) * 2
        else:
            hd[i] = (np.sin(2 * np.pi * cut1 * (n[i] - M/2.)) / (np.pi*(n[i] - M/2.))) \
            - (np.sin(2 * np.pi * cut2 * (n[i] - M/2.)) / (np.pi*(n[i] - M/2.)))
    return hd

""" Impulse response for high pass filter """
def ImpulseresponseHP(n, M, f):
    hd = np.zeros(len(n))
    for i in range(len(hd)):
        if n[i] == M/2:
            hd[i] = 1-f * 2
        else:
            hd[i] = (-np.sin(2 * np.pi * f * (n[i] - (M/2))))/(np.pi * (n[i] - (M/2)))
    return hd       

""" Impulse response for lowpass filter """
def ImpulseresponseLP(n, M, f): 
    hd = np.zeros(len(n))
    for i in range(len(hd)):
        if n[i] == M/2:
            hd[i] = f * 2
        else:
            hd[i] = (np.sin(2 * np.pi * f*(n[i] - (M/2))))/(np.pi * (n[i] - (M/2)))
    return hd

""" Impulse response for bandpass filter """
def ImpulseresponseBP(n,M,ft1,ft2):
    h = np.zeros(len(n))
    for i in range(len(n)):
        if i == M/2:
            h[i] = 2*(ft2 - ft1)
        else:
            h[i] = (1 / (np.pi*(i - M/2.)))*(np.sin(ft2*2*np.pi*(i - M/2.)) \
            - (np.sin(ft1*2*np.pi*(i - M/2.))))
    return h

#==============================================================================
# Adding noise
#==============================================================================

def add_noise(data,noise,c = 0.5):
    signal=np.zeros(len(data))    
    for i in range(len(data)):
        signal[i]=data[i]+float(c)*noise[i]
    return signal