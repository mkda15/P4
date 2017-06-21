# -*- coding: utf-8 -*-

#==============================================================================
# Imports
#==============================================================================

from fouriertransform import stft , db 
from windowfunctions import Kaiser
from peak_detection import peak_dec
import numpy as np
import impulseresponse as impuls
import matplotlib.pyplot as plt
import scipy.io.wavfile as siw

#==============================================================================
# Variable and data import 
#==============================================================================

""" Data import """
# Single tone with clap
freq , data  = siw.read('Lydfiler/enkelt_tone/forsoeg_enkelt_dyb.wav')   # Data signal
freq2, noise = siw.read('Lydfiler/stoej/klap_takt_2.wav')                # Noise signal

# TEST 1 
#freq , data  = siw.read('Lydfiler/skala/forsoeg_skala_hurtig.wav')     # Data signal
#freq2, noise = siw.read('Lydfiler/stoej/klap_takt_2.wav')              # Noise signal

# TEST 2 
#freq , data  = siw.read('Lydfiler/melodi/alene/forsoeg_lillepeteredderkop_langsom.wav')     # Data signal
#freq2, noise = siw.read('Lydfiler/stoej/klap_takt_2.wav')                                   # Noise signal

# TEST 3                                  
#freq , data  = siw.read('Lydfiler/akkorder/forsoeg_akkord_dyb2.wav')   # Data signal
#freq2, noise = siw.read('Lydfiler/stoej/klap_takt_2.wav')              # Noise signal

down_data = np.zeros(len(data/2))

down = 5

for i in range(len(data)/down):
    down_data[i] = data[down*i]

data = down_data

freq = freq/down

""" Lengths of data and noise aligned"""
if len(data) > len(noise):
    while len(data) > len(noise):
        noise = np.append(noise,noise)
if len(data) < len(noise):
        noise = noise[:len(data)]

print("Functions and data imported (1/8)")

""" Generate signal with noise """
signal = impuls.add_noise(data,noise,c = 1.0)   # Noise and data conjoined

print('Noise added (2/8)')

""" Variables for filter """
window = Kaiser     # The wanted window is named (has to be capitalised and has to be imported under windowfunctions)
cut1 = 70./freq     # First cut-off frequency
cut2 = 1000./freq   # Second cut-off frequency
samples = len(data) # Amount of samples in the signal (data points)

plotlength = int(samples/20) # Length for plotting (arbitrary)

""" Specifications for Kaiser window """
delta_1 = 0.05   # Peak approximation error in amplitude 
delta_2 = 10.    # Max transition width is 2*delta_2

""" Axis og linspaces """
t   = samples/float(freq)                   # The time for how long the system runs (for making the time axis)
time = np.linspace(0,t,samples)             # Axis for time domain
freq_axis = np.linspace(0,freq/2,samples/2) # Axis for frequency domain
freq_axis_norm = np.linspace(0,1,samples/2) 

""" Variables for spektogram """
freq_inter1 = 0    
freq_inter2 = 200

fontsize = 13

#==============================================================================
# Filter coefficients calculated and filter applied
#==============================================================================

""" Window function and impulse response """
if window == Kaiser:                            # If Kaiser window is chosen do this
    w,M,n,beta = window(delta_1,delta_2,freq)   # Retuns window funcion(time), order M and linspace n of length M+1  

else:                                           # If other window is chosen
    M = 1000.                                   # Set order of filter 
    n = np.linspace(0,M,M+1)
    w = window(n,M)

print('Window generated (3/8)')

hd = impuls.ImpulseresponseBP(n,M,cut1,cut2)    # Ideal impulse reponse 

h = hd * w                                      # Actual impulse response

print('Impulse response calculated (4/8)')


#==============================================================================
# Fourier transform of filter and data
#==============================================================================

""" Fourier transform of data """
H = np.fft.fft(h,(len(signal)))         # The Fourier transform of the final impulse response zero-padded to fit the signal
Hdb =  db(np.abs(H).T)                  # Frequency response of filter in dB
DATA = np.fft.fft(data)/len(signal)     # The Fourier transform of pure signal
NOISE = np.fft.fft(noise)               # The Fourier transform of noise
SIGNAL = np.fft.fft(signal)             # The Fourier transform of signal with noise

print('Fourier transform of data completed (5/8)')

""" Filtering of signal in the frequency domain """                  

SIGNAL_FILT = H * SIGNAL                # Convolution between the filter H and the noise SIGNALx
signal_filt = np.fft.ifft(SIGNAL_FILT)  # Filtered data
signal_filt = np.real(signal_filt)      # Cast to real, to remove the 0j from ifft.

print('Filtering of data completed (6/8)')

#==============================================================================
# Plt plots af alt det intresante og data gemmes
#==============================================================================

""" Impulse respose of filter """
#plt.plot(h[:sampels/2])
#plt.xlabel('Samples [n]')
#plt.ylabel('Amplitude')
#plt.axis([0,M+100,-0.015,0.046])
#plt.savefig("figures/filter_test/impulse.png")
#plt.show()

""" Frequency respose of filter """
plt.plot(freq_axis,np.abs(H)[:samples/2],'r')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
plt.axis([0,2500,0,1.1])
#plt.savefig("figures/filter_test/freq_response1.png")
plt.show()
#close up
f, axarr = plt.subplots(2, sharex=True)

axarr[0].plot(freq_axis[:plotlength], np.abs(H)[:plotlength],'r')
axarr[0].axis([20,130,0.9,1.1])
axarr[0].set_ylabel('Amplitude')
axarr[0].axvline((70-delta_2), color='green')   # Lower transition bound
axarr[0].axvline((70+delta_2), color='green')   # Upper transition bound
axarr[0].axhline((1+delta_1), color='green')    # Upper cut-off frequency bound
axarr[0].axhline((1-delta_1), color='green')    # Lower cut-off frequency bound 
 
plt.plot(freq_axis[:plotlength],np.abs(H)[:plotlength],'r')  
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
plt.axvline((70-delta_2), color='green')    # Lower transition bound
plt.axvline((70+delta_2), color='green')    # Upper transition bound
plt.axhline((0+delta_1), color='green')     # Lower cut-off frequency bound
plt.axis([20,130,0,0.2])
#plt.savefig("figures/filter_test/freq_response2.png")
plt.show()

""" dB representation of frequency response of filter """
plt.plot(freq_axis,Hdb[:samples/2])
plt.xlabel("Frequency [Hz]")
plt.ylabel("Amplitude [dB]")
plt.title('Frequency response of filter')
plt.axis([0,1500,-100,2])
plt.show() 

""" Signal with noise in the time domain """
plt.plot(time,signal)  # Original data with noise added 
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.title('Signal with noise')
#plt.axis([0,6,-1,1])
#plt.savefig("figures/integrationstest/signal.png")
plt.show()

""" Filtered signal in the time domain """
plt.plot(time,signal_filt)  # Original data with noise added 
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
plt.title('Filtered signal')
plt.axis([0,6,-6000,6000])
#plt.savefig("figures/integrationstest/f_signal.png")
plt.show()

""" Close-up on all data in the time domain """
plt.plot(time,signal, 'r-', label = "Noisy signal")
plt.plot(time,signal_filt, 'b-', label = "Filtered signal") 
plt.plot(time,data, 'g-',label = "Pure signal")                       
plt.legend(loc = 'upper right')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')
#plt.axis([1.03,1.04,-0.1,0.1])
plt.show()                                    

""" Pure signal with noise in the frequency domain """
plt.plot(freq_axis[:plotlength],np.abs(SIGNAL)[:plotlength])          
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
#plt.savefig("figures/filter_test/SIGNAL.png")
#plt.savefig("figures/integrationstest/FSIGNAL.png")
plt.show()

""" Filtered signal in the frequency domain """
plt.plot(freq_axis[:plotlength],np.abs(SIGNAL_FILT[:plotlength]))   
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
#plt.savefig("figures/filter_test/f_SIGNAL.png")
#plt.savefig("figures/integrationstest/f_FSIGNAL.png")
plt.show()

print('Signals plotted (7/8)')

#==============================================================================
# Computing Spectrogram
#==============================================================================

X,o,ws = stft(signal_filt)      # return STFT, stft and used window(time) 

print('STFT calculated (8/8)')

X = db(np.abs(X).T)             # Calculated to dB

x = np.linspace(0,time[-1],np.shape(X)[1])       
y = np.linspace(0,freq_axis[-1],np.shape(X)[0]) 

spec = plt.pcolormesh(x,y[freq_inter1:freq_inter2],X[freq_inter1:freq_inter2],cmap='jet')
cb   = plt.colorbar(spec)
cb.set_label(label = 'Amplitude [dB]', fontsize=fontsize)
plt.xlabel('Time [s]', fontsize = fontsize)
plt.ylabel('Frequency [Hz]', fontsize = fontsize)
plt.axis([0,time[-1],0,2000])
#plt.savefig("figures/skala.png")
#plt.savefig("figures/integrationstest/spectrogram.png")
#plt.savefig("figures/systemtest/final_spec1.png")
plt.show()

#==============================================================================
# Peak Dectection
#==============================================================================

max_freq_t = peak_dec(X,0.75,y) #limit is given as percentages of max amplitude in STFT

""" Plot peak dection """
plt.plot(x,max_freq_t,'o')
plt.xlabel('Time [s]')
plt.ylabel('Frequency [Hz]')
#plt.savefig("figures/peak/peak_lim4.png")
#plt.savefig("figures/integrationstest/peak_dec.png")
#plt.savefig("figures/systemtest/final_peak1.png")

#==============================================================================
# SNR
#==============================================================================

def RMS(x):
    x = x**2
    RMS = np.sqrt((np.sum(x))) / (len(x))
    return RMS

def SNR(signal,noise):
    SNR = ((RMS(signal))**2 / (RMS(noise))**2)
    return SNR

SNR = SNR(data,noise)
SNRdB = 10*np.log10(SNR)

print('SNR = %.0f dB' %SNRdB)