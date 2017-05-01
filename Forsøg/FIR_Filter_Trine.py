# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 21:50:53 2017

@author: Trine Nyholm Jensen
"""

import numpy as np
import matplotlib.pyplot as plt
#from scipy import signal
import scipy.io.wavfile as siw

#==============================================================================
# impotere data 
#==============================================================================
signal = siw.read('Lydfiler/forsoeg_nopeak/enkelt_tone/forsoeg_enkelt_dyb.wav')
noise = siw.read('Lydfiler/forsoeg_nopeak/stoej/klap_takt_2.wav')

freq = signal[0] #same freq for signal and noise
signal = signal[1]

noise = noise[1]

# fits lenght of noise to lenght of signal
if len(noise) < len(signal):
    noise = np.lib.pad(noise,(0,len(signal)-len(noise)),'constant',constant_values=0)
else:
    noise = noise[0:len(signal)]
    #data1 = np.lib.pad(data1,(0,len(data2)-len(data1)),'constant',constant_values=0)

signal = signal/1000  #normalized data
noise = noise/1000

#==============================================================================
# Filter specifications
#==============================================================================
#N = len(signal) #lenght of input 
#if N % 2 > 0:
#    raise ValueError('Brug N = potenser af 2')
N = 1000  #tjek om ovenstående er bedre altså længde af signal

fs = 44100.
#fs = 2000.         # sampling frequency 
if not type(fs) == float:
    raise ValueError("The sampling frequency should be a float.")
    
ft1 = 75 / fs      # cut off 1 
ft2 = 500 / fs     # cut off 2 

delta = 0.1 # peak approximation error in amplitude 

M1 = N-1 #selvvalgt filter orden 


#n = np.linspace(0,N,N)     
sampels = len(signal)        # evt.  sampels som data længden med virker ikke med lp
n = np.linspace(0,N,sampels)
freq_ax = np.linspace(0,fs/2,sampels/2)   # normalised frequency axis
#==============================================================================
# Ideal impulse response
#==============================================================================
def bp(n,M,ft1,ft2): # Bandpass filter
    h = np.zeros(len(n))
    for i in range(len(n)):
        if i == M/2:
            h[i] = 2*(ft2 - ft1)
        else:
            h[i] = (1 / (np.pi*(i - M/2.)))*(np.sin(ft2*2*np.pi*(i - M/2.)) \
            - (np.sin(ft1*2*np.pi*(i - M/2.))))
    return h

def lp(n,M,ft): # Lowpass filter
    h = np.zeros(len(n))
    for i in range(len(n)):
        if i == M/2.:
            h[i] = 2*ft
        else:
            h[i] = np.sin(2*np.pi*ft*(i - M/2.)) / (np.pi*(i - M/2.))
    return h


#==============================================================================
# Vinduer
#==============================================================================

def rect(n,M): # Det rektangulaere vindue
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] >= 0 and n[i] <= M:
            w[i] = 1
        else:
            w[i] = 0
    return w

def ha(n,M,a): # Hann-vindue, hvis a = 0.5. Hamming-vindue, hvis a = 0.54.
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] >= 0 and n[i] <= M:
            w[i] = a - (1 - a)*np.cos((2*np.pi*n[i])/M)
        else:
            w[i] = 0
    return w

def blackman(n,M): # Blackman-vinduet
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] >= 0 and n[i] <= M:
            w[i] = 0.42 - 0.5*np.cos((2*np.pi*n[i])/M) + 0.8*np.cos((4*np.pi*n[i])/M)
        else:
            w[i] = 0
    return w

#==============================================================================
# impulse and frequency response
#==============================================================================
freq_inter1 = 100
freq_inter2 = 20000

h_d = bp(n,M1,ft1,ft2)  #   ideal impulse response
w = rect(n,M1)      #   window

h = h_d * w         #   windowed impulse response  

H = np.abs(np.fft.fft(h))[0:sampels/2]
plt.plot(freq_ax,np.abs(H))
plt.axis([0,4000,0,1.2])
plt.xlabel("frequency [Hz]")
plt.show()

# Pure signal
SIGNAL = 2/float(len(signal))*np.abs(np.fft.fft(signal)[:len(signal)/2])
#plt.plot(freq_ax[freq_inter1:freq_inter2],SIGNAL[freq_inter1:freq_inter2])
#plt.show()

def add_noise(data,noise,c = 1.5): #kilde side 229 i DTSP
    signal=np.zeros(len(data))    
    for i in range(len(data)):
        signal[i]=data[i]+float(c)*noise[i]
    return signal

# signal + noise
data = add_noise(signal,noise)
DATA = 2/float(len(data))*np.abs(np.fft.fft(data)[:len(data)/2])


#==============================================================================
# Filtereing 
#==============================================================================
filt_DATA = H * DATA

#plt.plot(freq_ax[freq_inter1:freq_inter2],DATA[freq_inter1:freq_inter2])
#plt.show()
#
#plt.plot(freq_ax[freq_inter1:freq_inter2],filt_DATA[freq_inter1:freq_inter2])
#plt.axis([0,4000,0,0.35])
#plt.show()

plt.plot(signal)
plt.show()

filt_data = np.fft.ifft(filt_DATA)
plt.plot(n[:len(n)/2],filt_data)
#==============================================================================
# save signal
#==============================================================================
#siw.write('Lydfiler/forsoeg_nopeak/output/out_signal_.wav',freq,)
#siw.write('Lydfiler/forsoeg_nopeak/output/out_data.wav',freq,data)





