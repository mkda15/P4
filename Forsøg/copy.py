# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 08:31:00 2017

@author: cht15
"""

#==============================================================================
# Imports
#==============================================================================

from short_time_fourier_transform import stft , db , stft_h, variance_t
from windowfunctions import Kaiser
import numpy as np
import impulsrespons as impuls
import matplotlib.pyplot as plt
import scipy.io.wavfile as siw

#==============================================================================
# Variable og data import 
#==============================================================================
""" Data import """
freq , data  = siw.read('Lydfiler/forsoeg_nopeak/enkelt_tone/forsoeg_enkelt_dyb.wav')     # Data signal
freq2, noise = siw.read('Lydfiler/forsoeg_nopeak/stoej/klap_takt_2.wav')                  # Noise signal
                                  

""" Length of data and noise alings"""
if len(data) > len(noise):
    while len(data) > len(noise):
        noise = np.append(noise,noise)
if len(data) < len(noise):
        noise = noise[:len(data)]

print("variabler og data importeret 1/9")

""" Generate signal with noise """
signal = impuls.add_noise(data,noise,c = 1.0)   # Noise and data conjoined

print('støj adderet 2/9')

""" Variabler til filter """
window = Kaiser     # The wanted window is named (Has to be capitalised and has to be imported under windowfunctions)
cut1 = 75./freq     # Cut off frequency for band
cut2 = 1000./freq   # Cut off frequency for band 
sampels = len(data) # Amount of sampels in the signal (data points)
M    = 1000.        # Filter order

plotlength = int(sampels/30) # Length for plotting (arbitrary)


""" Til Kaiser vinduet """
delta_1 = 0.05 # peak approximation error in amplitude 
delta_2 = 10. # max transition width is 2*delta_2

""" Aksis og linspaces """
t   = sampels/float(freq)                   # The time for howlong the system runs (for making the time axis)
tid = np.linspace(0,t,sampels)              # Axis for time domain
freq_axis = np.linspace(0,freq/2,sampels/2) # Axis for frequency domain
freq_axis_norm = np.linspace(0,1,sampels/2)
n   = np.linspace(0,M,M+1) 

""" Variabler til spektogram """
freq_inter1 = 0    
freq_inter2 = 150

fontsize = 13
dataType = "Tabs" #Variable to peak detection, if the file is with chords dataType == Chords if its tabs dataType should be == Tabs
                 

#==============================================================================
# Filter koefficenter udregnes og filtere anvendes
#==============================================================================
""" Vindue funktion og impuls respons udregnes """

if window == Kaiser:                                # What to be returned for different windowfunctions
    w,M,n = window(delta_1,delta_2,freq)
#    w = w[:-1] ?
#    n = n[:-1]
else:
    w = window(n,M)

print('vindue generet 3/9')

hd = impuls.ImpulsresponsBP(n,M,cut1,cut2)

h = hd * w                              # The final impulsrespons


print('impuls respons udregnet 4/9')

data = data / float((np.max(signal)))       # Reduktion of amplitude
signal = signal / float((np.max(signal))) 

#==============================================================================
# Fourier Transformation of Filter and Data
#==============================================================================
""" Dataen fourier transformeres """
H = np.fft.fft(h,(len(signal)))         # The fourier transformed of the final impulsrespons zero padded to fit the signal
DATA = np.fft.fft(data)     # Pure signal in fourier
NOISE = np.fft.fft(noise)   # Noise in fourier
SIGNAL = np.fft.fft(signal) # Signal with noise in fourier

print('Data fourier transformeret 5/9')

""" Filtering af signal i frekvens domænet """                  
SIGNAL_FILT = H * SIGNAL                # Convolution between the filter H and the noise SIGNALx
signal_filt = np.fft.ifft(SIGNAL_FILT)  # Filtered data
signal_filt = np.real(signal_filt)      # Cast to real, to remove the 0j from ifft.

print('Data filtreret 6/9')


#==============================================================================
# Plt plots af alt det intresante og data gemmes
#==============================================================================

""" Close-up all data in time domain  """
##plt.plot(tid,signal, 'r-', label = "signal")  
##plt.plot(tid,signal_filt, 'b-', label = "filt signal")
##plt.plot(tid,data, 'g-',label = "ren signal")  # Original data with noise added 
##plt.legend(loc = 'upper right')
##plt.xlabel('Time [sec.]')
##plt.axis([1,1.02,-0.1,0.1])
##plt.show()

""" Signal with noise in time domain """
#plt.plot(tid,signal)  # Original data with noise added 
#plt.xlabel('Time [sec.]')
#plt.ylabel('Amplitude')
##plt.savefig("figures/integrationstest/signal.pdf")
#plt.show()

""" Filtered signal in time domain """
#plt.plot(tid,signal_filt)  # Original data with noise added 
#plt.xlabel('Time [sec.]')
#plt.ylabel('Amplitude')
##plt.axis([0,6,-1.5,1])
##plt.savefig("figures/integrationstest/signal_filt.pdf")
#plt.show()


#plt.plot(freq_axis[:plotlength],np.abs(SIGNAL)[:plotlength])          # FFT of clean data
#plt.xlabel('Frequency [Hz]')
#plt.ylabel('Amplitude')
##plt.axis([200,270,0,1000])
##plt.savefig("figures/integrationstest/f_signal.pdf")
#plt.show()
#
#plt.plot(freq_axis[:plotlength],np.abs(SIGNAL_FILT[:plotlength]))   # FFT of the filtered data
#plt.xlabel('Frequency [Hz]')
#plt.ylabel('Amplitude')
##plt.savefig("figures/integrationstest/f_signal_filt.pdf")
#plt.show()
##
##plt.plot(tid,signal, 'r-', label = "signal")
##plt.plot(tid,signal_filt, 'b-', label = "filt") 
##plt.plot(tid,data, 'g-',label = "ren signal")                                                 # Original data with noise added 
##plt.legend(loc = 'upper right')
##plt.xlabel('Time [sec.]')
##plt.axis([1.03,1.06,-0.1,0.1])
##plt.show()
#
#                                        

print('plot plotteret 7/9')


#==============================================================================
# Spectrogram
#==============================================================================

X,o,ws = stft(signal_filt,fftsize = 2048 ,overlap = 2)     # STFT calculated

W = np.abs(np.fft.fft(ws))
#X,v_w,v_t,t = stft_h(signal_filt,overlap = 2)


print('stft udregnet 9/9') 

X = db(np.abs(X).T) 
#G = db(np.abs(G).T)                            # Calculated to dB

x = np.linspace(0,tid[-1],np.shape(X)[1])       
y = np.linspace(0,freq_axis[-1],np.shape(X)[0]) 


spec = plt.pcolormesh(x,y[freq_inter1:freq_inter2],X[freq_inter1:freq_inter2],cmap='hot')
cb   = plt.colorbar(spec)
cb.set_label(label = 'Amplitude (dB)', fontsize=fontsize)
plt.xlabel('Time (sec.)', fontsize = fontsize)
plt.ylabel('Frequency (Hz)', fontsize = fontsize)
#plt.axis([0,25,0,3000])
#plt.savefig("figures/skala.png")
#plt.savefig("figures/integrationstest/spectrogram.pdf")
#plt.savefig("figures/systemtest/final_spec.pdf")
plt.show()

plt.plot(freq_axis,np.angle(H)[:sampels/2])
plt.axis([0,1075,-4,4])
plt.show()

Hdb =  db(np.abs(H).T)

plt.plot(freq_axis,Hdb[:sampels/2])
plt.axis([0,1300,-100,2])
plt.show()


#==============================================================================
# Peak Dectection
#==============================================================================

X = X.T
p = 20 # lower limit for amplitude to be detected, below p -> 0 
# kan laves til en definition og placeres i et andet dokument.
sortedX = np.zeros(len(X),dtype = object)
for i in range(len(X)):

    sortedX[i] = np.sort(X[i])
if dataType == "Tabs": #Tjeck if data is in single tabs or chords
    max_freq_pos = np.zeros(len(X))
    for i in range(len(X)):
        if np.max(X[i]) > p:
            a = np.where(X[i][:] == np.max(X[i]))
            max_freq_pos[i] = a[0][0]
        else:
            max_freq_pos[i] = 0
    max_freq_t = np.zeros(len(X))
    for i in range(len(X)):
        if max_freq_pos[i] == 0:
            max_freq_t[i] = 0
        else:
            max_freq_t[i] = y[int(max_freq_pos[i])]
            
    plt.stem(x,max_freq_t)
    plt.xlabel('Time (sec.)')
    plt.ylabel('Frequency (Hz)')
    #plt.savefig("figures/integrationstest/peak_dec.pdf")
    #plt.savefig("figures/systemtest/final_peak.pdf")
  #  print(max_freq_t[6])
elif dataType == "Chords":
    max_freq_pos1 = np.zeros(len(X))
    max_freq_pos2 = np.zeros(len(X))
    max_freq_pos3 = np.zeros(len(X))

    for i in range(len(X)):
        if sortedX[i][-1] > 20:
            a = np.where(X[i][:] == sortedX[i][-1])
        else:
            a = [[0]]
        if sortedX[i][-2] > 20:
            b = np.where(X[i][:] == sortedX[i][-2])
        else:
            b = [[0]]
        if sortedX[i][-3] > 20:
            c = np.where(X[i][:] == sortedX[i][-3])
        else:
            c = [[0]]
        max_freq_pos1[i] = a[0][0]
        max_freq_pos2[i] = b[0][0]
        max_freq_pos3[i] = c[0][0]

    max_freq_t1 = np.zeros(len(X))
    max_freq_t2 = np.zeros(len(X))
    max_freq_t3 = np.zeros(len(X))
    for i in range(len(X)):
        if max_freq_pos1[i] == 0:
            max_freq_t1[i] = 0
        else:
            max_freq_t1[i] = y[int(max_freq_pos1[i])]
        if max_freq_pos2[i] == 0:
            max_freq_t2[i] = 0
        else:
            max_freq_t2[i] = y[int(max_freq_pos2[i])]
        if max_freq_pos3[i] == 0:
            max_freq_t3[i] = 0
        else:
            max_freq_t3[i] = y[int(max_freq_pos3[i])]
    plt.plot(max_freq_t1)
    plt.plot(max_freq_t2)
    plt.plot(max_freq_t3)

    sted = 75
    print(max_freq_t1[sted])
    print(max_freq_t2[sted])
    print(max_freq_t3[sted])

##
#k = max_freq_t
#l =np.zeros(len(k))
#l[0]=k[0]
#for i in range(len(k)-1):
#    if k[i+1] == k[i] :
#        l[i+1]= k[i+1]
#    elif k[i+1] == k[i+2]:
#        l[i+1]= k[i+2]
#    else :
#        l[i+1]=0
        
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
print(SNRdB)   

#==============================================================================
# Heisenberg 
#==============================================================================



            
