# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 08:31:00 2017

@author: cht15
"""

import numpy as np
from windowfunctions import Hamming, Hanning
import matplotlib.pyplot as plt
import scipy.signal as ss
import scipy.io.wavfile as siw
import winsound
#==============================================================================
# Variable
#==============================================================================
cut = np.pi/35.
cut1 = np.pi/50.
cut2 = np.pi/45.

freq , data  = siw.read('Lydfiler/forsoeg_nopeak/enkelt_tone/forsoeg_enkelt_lys.wav')
freq2, noise = siw.read('Lydfiler/forsoeg_nopeak/stoej/klap_random_1.wav')
#freq3, signal = siw.read('Lydfiler/noise_pc.wav')

#==============================================================================
# Filter funktion defineres
#==============================================================================
""" Filte koefficenter med filter orden M = 92 genereres"""
M = 1000

n = np.linspace(0,M,M+1)
#x = np.linspace(-np.pi,np.pi,len(n))

def ImpulsresponsBS(n,M,cut1,cut2): 
    hd = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] == M/2:
            hd[i] = 1 - (cut2 - cut1)/np.pi
        else:
            hd[i] = (np.sin(cut1*(n[i] - M/2.)) / (np.pi*(n[i] - M/2.))) \
            - (np.sin(cut2*(n[i] - M/2.)) / (np.pi*(n[i] - M/2.)))
    return hd

def ImpulsresponsHP(n, M, f):     # Den ønskværdige impulsrespons
    hd = np.zeros(len(n))
    for i in range(len(hd)):
        if n[i] == M/2:
            hd[i] = 1-f/np.pi
        else:
            hd[i] = (-np.sin(f*(n[i] - (M/2))))/(np.pi*(n[i] - (M/2)))
    return hd       

def ImpulsresponsLP(n, M, f):     # Den ønskværdige impulsrespons
    hd = np.zeros(len(n))
    for i in range(len(hd)):
        if n[i] == M/2:
            hd[i] = f/np.pi
        else:
            hd[i] = (np.sin(f*(n[i] - (M/2))))/(np.pi*(n[i] - (M/2)))
    return hd 

def add_noise(data,noise,c = 0.5): #kilde side 229 i DTSP
#    signal = np.convolve(data,noise)
#    return signal
    signal=np.zeros(len(data))    
    for i in range(len(data)):
        signal[i]=data[i]+float(c)*noise[i]
    return signal
#==============================================================================
# Filter koefficenter udregnes
#==============================================================================
noise = noise[:len(data)]
signal = add_noise(data,noise,c = 0.50)

w = Hanning(n,M) #Hanning eller Hamming for nu
#hd = ImpulsresponsBS(n,M,cut1,cut2)
#hd = ImpulsresponsHP(n,M,cut)
hd = ImpulsresponsLP(n,M,cut)
h = hd * w
H = np.fft.fft(h,(len(signal)))


DATA = np.fft.fft(data)     # Pure signal in fourier
NOISE = np.fft.fft(noise)   # Noise in fourier
SIGNAL = np.fft.fft(signal) # Signal with noise in fourier

SIGNAL_FILT = H * DATA      # Convolution between the filter H and the noise SIGNAL
signal_filt = np.fft.ifft(SIGNAL_FILT) # Filtered data
signal_filt = np.real(signal_filt)
#==============================================================================
# Plt
#==============================================================================

plt.plot(np.abs(DATA)[:5500])
plt.show()
plt.plot(np.abs(H)[:5500])
plt.show()
plt.plot(np.abs(NOISE)[:6500])
plt.show()
plt.plot(data)
plt.show()
plt.plot(noise)
plt.show()
plt.plot(signal)
plt.show()
plt.plot(signal_filt)
plt.show()
plt.plot(np.abs(SIGNAL_FILT[:6500]))

#winsound.PlaySound('Lydfiler/forsoeg_nopeak/skala/forsoeg_skala_hurtig.wav', winsound.SND_FILENAME)

siw.write('Lydfiler/forsoeg_nopeak/output/out_signal_filt.wav',freq,signal_filt)
siw.write('Lydfiler/forsoeg_nopeak/output/out_data.wav',freq,data)
siw.write('Lydfiler/forsoeg_nopeak/output/out_noise.wav',freq,noise)
siw.write('Lydfiler/forsoeg_nopeak/output/out_signal.wav',freq,signal)



