# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 08:31:00 2017

@author: cht15
"""

#==============================================================================
# Imports
#==============================================================================

from short_time_fourier_transform import stft , db
from windowfunctions import Hamming, Hanning, Kaiser
import numpy as np
import impulsrespons as impuls
import matplotlib.pyplot as plt
#import scipy.signal as ss
import scipy.io.wavfile as siw

#==============================================================================
# Variable og data import 
#==============================================================================
""" Data import """
freq , data  = siw.read('Lydfiler/forsoeg_nopeak/enkelt_tone/forsoeg_enkel_dyb.wav')  # Data signal
freq2, noise = siw.read('Lydfiler/forsoeg_nopeak/stoej/klap_takt_2.wav')                  # Noise signal
                    #freq3, signal = siw.read('Lydfiler/noise_pc.wav')                      # Noise and data as a single file

""" Length of data and noise alings"""
if len(data) < len(noise):
    noise = noise[:len(data)]
elif len(data) > len(noise):
    data = data[:len(noise)]

""" Variabler til filter """
window = Kaiser     # The wanted window is named (Has to be capitalised and has to be imported under windowfunctions)
M    = 1000.        # Filter order
cut  = 1000./freq   # Cut off frequency
cut1 = 180./freq     # Cut off frequency for band
cut2 = 500./freq    # Cut off frequency for band 
sampels = len(data) # Amount of sampels in the signal (data points)
plotlength = int(sampels/2) # Length for plotting (arbitrary)


""" Til Kaiser vinduet """
delta_1 = 0.05 # peak approximation error in amplitude 
delta_2 = 4. # max transition width is 2*delta_2

""" Aksis og linspaces """
t   = sampels/float(freq)                   # The time for howlong the system runs (for making the time axis)
tid = np.linspace(0,t,sampels)              # Axis for time domain
freq_axis = np.linspace(0,freq/2,sampels/2) # Axis for frequency domain
n   = np.linspace(0,M,M+1)                  # Integer numbers for making window and impulsrespons

""" Variabler til spektogram """
freq_inter1 = 0    
freq_inter2 = 100

fontsize = 13
                 
print("variabler og data importeret 1/9")

signal = impuls.add_noise(data,noise,c = 1.0)   # Noise and data conjoined

print('støj adderet 2/9')

#==============================================================================
# Filter koefficenter udregnes og filtere anvendes
#==============================================================================
""" Vindue funktion og impuls respons udregnes """

if window == Kaiser:                                # What to be returned for different windowfunctions
    w,M,beta,A,n = window(delta_1,delta_2,freq)
    w = w[:-1]
    n = n[:-1]
else:
    w = window(n,M)

print('vindue generet 3/9')

hd = impuls.ImpulsresponsBP(n,M,cut1,cut2)
#hd = impuls.ImpulsresponsHP(n,M,cut)
#hd = impuls.ImpulsresponsLP(n,M,cut)    # Desired impulsrespons

h = hd * w                              # The final impulsrespons
H = np.fft.fft(h,(len(signal)))         # The fourier transformed of the final impulsrespons zero padded to fit the signal

print('impuls respons udregnet 4/9')

signal = signal / float((np.max(signal))) # Reduktion of amplitude

""" Dataen fourier transformeres """
DATA = np.fft.fft(data)     # Pure signal in fourier
NOISE = np.fft.fft(noise)   # Noise in fourier
SIGNAL = np.fft.fft(signal) # Signal with noise in fourier

print('Data fourier transformeret 5/9')
                   
SIGNAL_FILT = H * SIGNAL                # Convolution between the filter H and the noise SIGNALx
signal_filt = np.fft.ifft(SIGNAL_FILT)  # Filtered data
signal_filt = np.real(signal_filt)      # Cast to real, to remove the 0j from ifft.

print('Data filtreret 6/9')


#==============================================================================
# Plt plots af alt det intresante og data gemmes
#==============================================================================

plt.plot(freq_axis[:plotlength],np.abs(DATA)[:plotlength])          # FFT of clean data
plt.xlabel('Ren Data freq')
plt.show()

plt.plot(freq_axis[:plotlength],np.abs(H)[:plotlength])             # FFT of impuls respons   
plt.xlabel('Ren stoej freq')
plt.show()

plt.plot(freq_axis[:2000],np.abs(NOISE)[:2000])         # FFT of the noise
plt.xlabel('Stoej og data freq')
plt.show()

plt.plot(tid,data)                                                  # Original data
plt.xlabel('Ren Data tid')
plt.show()

plt.plot(tid,noise)                                                 # Original noise
plt.xlabel('Ren Stoej tid')
plt.show()

plt.plot(tid,signal)                                                # Original data with noise added 
plt.xlabel('Stoej og data tid')
plt.show()

plt.plot(tid,signal_filt)                                           # The filtered data
plt.xlabel('Filtreret signal')
plt.show()

plt.plot(freq_axis[:plotlength],np.abs(SIGNAL_FILT[:plotlength]))   # FFT of the filtered data
plt.xlabel('Freq af Filtreret signal')
plt.show()

print('plot plotteret 7/9')


""" Data gemmes """
#siw.write('Lydfiler/forsoeg_nopeak/output/out_signal_filt.wav',freq,signal_filt)    # The filtered data is saved
#siw.write('Lydfiler/forsoeg_nopeak/output/out_data.wav',freq,data)                  # Original noise is saved with same length as data
#siw.write('Lydfiler/forsoeg_nopeak/output/out_noise.wav',freq,noise)                # Original data is saved with same length as noise
#siw.write('Lydfiler/forsoeg_nopeak/output/out_signal.wav',freq,signal)              # The signal with noise is saved

print('Data gemt 8/9')

#==============================================================================
# Spectrogram
#==============================================================================

X = stft(signal_filt,fftsize = 2500,overlap = 2)     # STFT calculated

print('stft udregnet 9/9')

X = db(np.abs(X).T)                             # Calculated to dB

x = np.linspace(0,tid[-1],np.shape(X)[1])       
y = np.linspace(0,freq_axis[-1],np.shape(X)[0]) 


spec = plt.pcolormesh(x,y[freq_inter1:freq_inter2],X[freq_inter1:freq_inter2],cmap='hot')
cb   = plt.colorbar(spec)
cb.set_label(label = 'Amplitude (dB)', fontsize=fontsize)
plt.xlabel('Time (sec)', fontsize = fontsize)
plt.ylabel('Frequency (Hz)', fontsize = fontsize)
plt.show()

X = stft(signal,fftsize = 2500,overlap = 2)     # STFT calculated


X = db(np.abs(X).T)                             # Calculated to dB

x = np.linspace(0,tid[-1],np.shape(X)[1])       
y = np.linspace(0,freq_axis[-1],np.shape(X)[0]) 


spec = plt.pcolormesh(x,y[freq_inter1:freq_inter2],X[freq_inter1:freq_inter2],cmap='hot')
cb   = plt.colorbar(spec)
cb.set_label(label = 'Amplitude (dB)', fontsize=fontsize)
plt.xlabel('Time (sec)', fontsize = fontsize)
plt.ylabel('Frequency (Hz)', fontsize = fontsize)
plt.show()
X = stft(data,fftsize = 2500,overlap = 2)     # STFT calculated


X = db(np.abs(X).T)                             # Calculated to dB

x = np.linspace(0,tid[-1],np.shape(X)[1])       
y = np.linspace(0,freq_axis[-1],np.shape(X)[0]) 


spec = plt.pcolormesh(x,y[freq_inter1:freq_inter2],X[freq_inter1:freq_inter2],cmap='hot')
cb   = plt.colorbar(spec)
cb.set_label(label = 'Amplitude (dB)', fontsize=fontsize)
plt.xlabel('Time (sec)', fontsize = fontsize)
plt.ylabel('Frequency (Hz)', fontsize = fontsize)
