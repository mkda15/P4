# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 08:31:00 2017

@author: cht15
"""

#==============================================================================
# Imports
#==============================================================================

from short_time_fourier_transform import stft , db
from windowfunctions import Kaiser, Hamming, Hanning
import numpy as np
import impulsrespons as impuls
import matplotlib.pyplot as plt
#import scipy.signal as ss
import scipy.io.wavfile as siw

#==============================================================================
# Variable og data import 
#==============================================================================
""" Data import """
freq , data  = siw.read('Lydfiler/forsoeg_nopeak/enkelt_tone/forsoeg_enkelt_dyb.wav')  # Data signal
freq2, noise = siw.read('Lydfiler/forsoeg_nopeak/stoej/kroelle_stoej.wav')                  # Noise signal
                    #freq3, signal = siw.read('Lydfiler/noise_pc.wav')                      # Noise and data as a single file

""" Length of data and noise alings"""
if len(data) < len(noise):
    noise = noise[:len(data)]
elif len(data) > len(noise):
    data = data[:len(noise)]

""" Variabler til filter """
M    = 8000         # Filter order
cut  = 1000./freq   # Cut off frequency
cut1 = 75./freq     # Cut off frequency for band
cut2 = 1000./freq    # Cut off frequency for band 
sampels = len(data) # Amount of sampels in the signal (data points)
plotlength = int(sampels/20) # Length for plotting (arbitrary)

""" Aksis og linspaces """
t   = sampels/float(freq)                   # The time for howlong the system runs (for making the time axis)
tid = np.linspace(0,t,sampels)              # Axis for time domain
freq_axis = np.linspace(0,freq/2,sampels/2) # Axis for frequency domain
n   = np.linspace(0,M,M+1)                  # Integer numbers for making window and impulsrespons

""" Variabler til spektogram """
freq_inter1 = 0    
freq_inter2 = 100

#fontsize = 13
                 
print("variabler og data importeret 1/9")

signal = impuls.add_noise(data,noise,c = 1.0)   # Noise and data conjoined

print('støj adderet 2/9')

#==============================================================================
# Filter koefficenter udregnes og filtere anvendes
#==============================================================================
""" Vindue funktion og impuls respons udregnes """
w = Hamming(n,M)        #Hanning eller Hamming for nu andre kan vælge impoteret i starten
hd = impuls.ImpulseresponseBP(n,M,cut1,cut2)
##########    kaiser window   ##########
delta_1= 0.01
delta_2 = float(0.5)         # half transition width in Hz

w1,M,beta,A,n = Kaiser(delta_1,delta_2,freq)            # window

########################################

print('vindue generet 3/9')

hd1 = impuls.ImpulseresponseBP(n,M,cut1,cut2)
#hd = impuls.ImpulsresponsHP(n,M,cut)
#hd = impuls.ImpulsresponsLP(n,M,cut)    # Desired impulsrespons

h = hd * w 
h = h[:len(h)-1] 
h1 = hd1 * w1 
h1 = h1[:len(h1)-1]
                           # The final impulsrespons
H = np.fft.fft(h,(len(signal)))         # The fourier transformed of the final impulsrespons zero padded to fit the signal
H1 = np.fft.fft(h1,(len(signal)))


signal = signal / float((np.max(signal))) # Reduktion of amplitude

""" Dataen fourier transformeres """
DATA = np.fft.fft(data)     # Pure signal in fourier
NOISE = np.fft.fft(noise)   # Noise in fourier
SIGNAL = np.fft.fft(signal) # Signal with noise in fourier

                   
SIGNAL_FILT = H * SIGNAL                # Convolution between the filter H and the noise SIGNALx
signal_filt = np.fft.ifft(SIGNAL_FILT)  # Filtered data
signal_filt = np.real(signal_filt)      # Cast to real, to remove the 0j from ifft.


#==============================================================================
# Plt plots af alt det intresante og data gemmes
#==============================================================================

plt.style.use('classic')

plt.plot(freq_axis[:plotlength],np.abs(H1)[:plotlength],'r') 
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
#plt.axis([0,5000,0,1.05])
plt.savefig("figures/filter_test/freq_response1.pdf")
plt.show()

f, axarr = plt.subplots(2, sharex=True)

axarr[0].plot(freq_axis[:plotlength], np.abs(H1)[:plotlength],'r')
axarr[0].axis([70,80,0.975,1.025])
axarr[0].set_ylabel('Amplitude')
axarr[0].axvline((75-delta_2), color='green') # Nedre transitionsgrænse
axarr[0].axvline((75+delta_2), color='green') # Øvre transitionsgrænse
axarr[0].axhline((1+delta_1), color='green') # Nedre knækfrekvens
axarr[0].axhline((1-delta_1), color='green') # Nedre knækfrekvens 
 
plt.plot(freq_axis[:plotlength],np.abs(H1)[:plotlength],'r')  
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
plt.axvline((75-delta_2), color='green') # Nedre transitionsgrænse
plt.axvline((75+delta_2), color='green') # Øvre transitionsgrænse
plt.axhline((0+delta_1), color='green') # Nedre knækfrekvens
plt.axis([70,80,0,0.05])
plt.savefig("figures/filter_test/freq_response2.pdf")
plt.show()

plt.plot(freq_axis[:plotlength],np.abs(SIGNAL)[:plotlength],'b')          # FFT of clean data
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
plt.axis([60,80, 0,200])
#plt.savefig("figures/filter_test/SIGNAL.pdf")
plt.show()

plt.plot(freq_axis[:plotlength],np.abs(SIGNAL_FILT[:plotlength]),'b')   # FFT of the filtered data
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
plt.axis([60,80, 0,200])
#plt.savefig("figures/filter_test/filt_SIGNAL.pdf")
plt.show()

#plt.plot(freq_axis[:plotlength],np.abs(NOISE)[:plotlength])         # FFT of the noise
#plt.title('frequency response of noise')           # FFT of impuls respons   
#plt.xlabel('frequency [Hz]')
#plt.show()

#plt.plot(tid,data)                                                  # Original data
#plt.title('Original data over time')           # FFT of impuls respons   
#plt.xlabel('Time [sec.]')
#plt.show()

#plt.plot(tid,noise)                                                 # Original noise
#plt.title('Original noise over time')           # FFT of impuls respons   
#plt.xlabel('Time [sec.]')
#plt.show()

#plt.plot(tid,signal)                                                # Original data with noise added 
#plt.title('Noisy signal over time')           # FFT of impuls respons   
#plt.xlabel('Time[Hz]')
#plt.show()
#
#plt.plot(tid,signal_filt)                                           # The filtered data
#plt.title('Filtreret signal over time ')           # FFT of impuls respons   
#plt.xlabel('Time [sec.]')
#plt.show()


#""" Data gemmes """
#siw.write('Lydfiler/forsoeg_nopeak/output/out_signal_filt.wav',freq,signal_filt)    # The filtered data is saved
#siw.write('Lydfiler/forsoeg_nopeak/output/out_data.wav',freq,data)                  # Original noise is saved with same length as data
#siw.write('Lydfiler/forsoeg_nopeak/output/out_noise.wav',freq,noise)                # Original data is saved with same length as noise
#siw.write('Lydfiler/forsoeg_nopeak/output/out_signal.wav',freq,signal)              # The signal with noise is saved
#
#print('Data gemt 8/9')

#==============================================================================
# Spectrogram
#==============================================================================
#
#X = stft(signal_filt,fftsize = 2500,overlap = 2)     # STFT calculated
#
#print('stft udregnet 9/9')
#
#X = db(np.abs(X).T)                             # Calculated to dB
#
#x = np.linspace(0,tid[-1],np.shape(X)[1])       
#y = np.linspace(0,freq_axis[-1],np.shape(X)[0]) 
#
#
#spec = plt.pcolormesh(x,y[freq_inter1:freq_inter2],X[freq_inter1:freq_inter2],cmap='hot')
#cb   = plt.colorbar(spec)
#cb.set_label(label = 'Amplitude (dB)', fontsize=fontsize)
#plt.xlabel('Time (sec)', fontsize = fontsize)
#plt.ylabel('Frequency (Hz)', fontsize = fontsize)
