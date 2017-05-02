# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 08:31:00 2017

@author: cht15
"""

import numpy as np
from windowfunctions import Hamming, Hanning
import impulsrespons as impuls
import matplotlib.pyplot as plt
import scipy.signal as ss
import scipy.io.wavfile as siw
#import winsound
#==============================================================================
# Variable
#==============================================================================
freq , data  = siw.read('Lydfiler/forsoeg_nopeak/enkelt_tone/forsoeg_enkelt_lys.wav')
freq2, noise = siw.read('Lydfiler/forsoeg_nopeak/stoej/klap_random_1.wav')
#freq3, signal = siw.read('Lydfiler/noise_pc.wav')
noise = noise
data = data
cut  = 400./freq
cut1 = 50./freq
cut2 = 500./freq
sampels = len(data)
#==============================================================================
# Filter funktion defineres
#==============================================================================
""" Filte koefficenter med filter orden M = 1000 genereres"""
M = 1000
freq_axis = np.linspace(0,freq/2,sampels/2)
n = np.linspace(0,M,M+1)
#x = np.linspace(-np.pi,np.pi,len(n))
#==============================================================================
# Filter koefficenter udregnes
#==============================================================================
noise = noise[:len(data)]
signal = impuls.add_noise(data,noise,c = 0.50)

w = Hamming(n,M) #Hanning eller Hamming for nu
#hd = impuls.ImpulsresponsBS(n,M,cut1,cut2)
#hd = impuls.ImpulsresponsHP(n,M,cut)
hd = impuls.ImpulsresponsLP(n,M,cut)
h = hd * w
H = np.fft.fft(h,(len(signal)))

signal = signal / float((np.max(signal)))

DATA = np.fft.fft(data)     # Pure signal in fourier
NOISE = np.fft.fft(noise)   # Noise in fourier
SIGNAL = np.fft.fft(signal) # Signal with noise in fourier

                   
SIGNAL_FILT = H * SIGNAL      # Convolution between the filter H and the noise SIGNALx
signal_filt = np.fft.ifft(SIGNAL_FILT) # Filtered data
signal_filt = np.real(signal_filt)
#==============================================================================
# Plt
#==============================================================================

plt.plot(freq_axis[:5500],np.abs(DATA)[:5500])
plt.show()
plt.plot(freq_axis[:5500],np.abs(H)[:5500])
plt.show()
plt.plot(freq_axis[:6500],np.abs(NOISE)[:6500])
plt.show()
plt.plot(data)
plt.show()
plt.plot(noise)
plt.show()
plt.plot(signal)
plt.show()
plt.plot(signal_filt)
plt.show()
plt.plot(freq_axis[:6500],np.abs(SIGNAL_FILT[:6500]))

#winsound.PlaySound('Lydfiler/forsoeg_nopeak/skala/forsoeg_skala_hurtig.wav', winsound.SND_FILENAME)

siw.write('Lydfiler/forsoeg_nopeak/output/out_signal_filt.wav',freq,signal_filt)
siw.write('Lydfiler/forsoeg_nopeak/output/out_data.wav',freq,data)
siw.write('Lydfiler/forsoeg_nopeak/output/out_noise.wav',freq,noise)
siw.write('Lydfiler/forsoeg_nopeak/output/out_signal.wav',freq,signal)



