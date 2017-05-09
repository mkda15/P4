# -*- coding: utf-8 -*-
"""
Created on Tue May  9 09:37:39 2017

@author: Trine
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io.wavfile as siw
import scipy as sc

def bp_filter (d1, d2, fs, ft1,ft2):
    
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
        
        n = np.linspace(0,M,M+1) 
        w = np.zeros(len(n))
        for i in range(len(n)):
            if n[i] >= 0 and n[i] <= M:
                variabel = beta*np.sqrt(1-((n[i]-(M/2))/(M/2))**2)
                w[i] = sc.special.i0(variabel)/sc.special.i0(beta)
            else:
                w[i] = 0
       
        if (M%2 == 1):
            M=M+1
        
        return w,M,n
    
    w,M,n = Kaiser(d1, d2, fs)
    
    
    def bp(n,M,ft1,ft2): # Bandpass filter
        h = np.zeros(len(n))
        for i in range(len(n)):
            if i == M/2:
                h[i] = 2*(ft2 - ft1)
            else:
                h[i] = (1 / (np.pi*(i - M/2.)))*(np.sin(ft2*2*np.pi*(i - M/2.)) \
                - (np.sin(ft1*2*np.pi*(i - M/2.))))
        return h
        
    h_d = bp(n,M,ft1,ft2)
    
    h = h_d*w

    return h
    
 


#==============================================================================
# Eksempel
#==============================================================================
fs = 44100.         # sampling frequency
         
if not type(fs) == float:
    raise ValueError("The sampling frequency should be a float.")
    
ft1 = 75 / fs      # cut off 1 
ft2 = 500 / fs     # cut off 2 

delta_1 = 0.05 # peak approximation error in amplitude 
delta_2 = 10  # max transition width in Hz is 2*delta_2
   
h = bp_filter(delta_1,delta_2,fs,ft1,ft2)
