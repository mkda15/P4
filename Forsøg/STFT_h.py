# -*- coding: utf-8 -*-
"""
Created on Sun May 14 17:39:40 2017

@author: Trine
"""
import numpy as np
import scipy.io.wavfile as siw


freq , data  = siw.read('Lydfiler/forsoeg_nopeak/enkelt_tone/forsoeg_enkelt_dyb.wav')

#==============================================================================
# STFT - optize lenght by Heisenberg
#==============================================================================
def stft_h(signal, overlap = 2):   
    fftsize = 500
    wlist = [0]
    tlist = [0]
    for k in range(5000):
        hop = fftsize / overlap
        w = np.kaiser(fftsize+1,4)[:-1]      # better reconstruction with this trick +1)[:-1]  
    
        x =  np.array([w*signal[i:i+fftsize] for i in range(0, len(signal)-fftsize, hop)])
        X =  np.array([np.fft.rfft(w*signal[i:i+fftsize]) for i in range(0, len(signal)-fftsize, hop)])    
        
        v_w = np.var(X)
        v_t = np.var(x) 
        
        wlist = np.append(wlist,[v_w]) 
        tlist = np.append(tlist,[v_t])
        if v_t*v_w < 1/4.:
            X = [0]
            
        elif v_t*v_w > 1/4.:
            fftsize += 10
   
    print fftsize 
    print wlist        
    return X
 
X = stft_h(data[:500],overlap=2)  