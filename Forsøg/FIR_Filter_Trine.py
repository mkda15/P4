# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 21:50:53 2017

@author: Trine Nyholm Jensen
"""

import numpy as np
import matplotlib.pyplot as plt
#from scipy import signal
import scipy.io.wavfile as siw
import scipy as sc
from windowfunctions import Hamming, Hanning, Blackman

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

signal = signal/1000.  #normalized data
noise = noise/1000.


t = len(signal)/freq
lin = np.linspace(0,t,len(signal)) # time axis 

#==============================================================================
#  Add noise to signal
#==============================================================================

def add_noise(signal,noise,c = 1): #kilde side 229 i DTSP
    data = np.zeros(len(signal))    
    for i in range(len(signal)):
        data[i]=signal[i]+float(c)*noise[i]
    return data

# signal + noise = data
data = add_noise(signal,noise)

#==============================================================================
# Plot in time domain   
#==============================================================================
#start = 0 # start udsnit
#slut = start+len(signal) # slut udsnit
#plt.plot(lin[start:slut],data[start:slut], label='data')
##plt.plot(lin[start:slut],signal[start:slut], label='signal')
##plt.plot(lin[start:slut],noise[start:slut], label='noise')
#
#plt.xlabel('Time[Seconds]')
#plt.ylabel('Amplitude [Voltages]')
#plt.legend(loc= "lower right")
#plt.show()


#==============================================================================
# Filter specifications
#==============================================================================
#N = len(signal) #lenght of input 
#if N % 2 > 0:
#    raise ValueError('Brug N = potenser af 2')
N = 1001  #tjek om ovenstående er bedre altså længde af signal

fs = 44100.         # sampling frequency
         
if not type(fs) == float:
    raise ValueError("The sampling frequency should be a float.")
    
ft1 = 75 / fs      # cut off 1 
ft2 = 500 / fs     # cut off 2 

delta_1 = 0.05 # peak approximation error in amplitude 
delta_2 = 1  # max transition width in Hz is 2*delta_2

M1 = N-1 #selvvalgt filter orden 
    
sampels = len(signal)
n = np.linspace(0,M1,N) 
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
    
    return w,M,beta,A,n,tw


#==============================================================================
# impulse and frequency response
#==============================================================================
freq_inter1 = 100           # frequency inverval
freq_inter2 = 20000
#
#h_d = bp(n,M1,ft1,ft2)
#w = rect(n,M1)
#h = h_d * w


##########    kaiser window   ##########
w = Kaiser(delta_1,delta_2,fs)[0]             # window
M = Kaiser(delta_1,delta_2,fs)[1]             # order 
beta = Kaiser(delta_1,delta_2,fs)[2]
A = Kaiser(delta_1,delta_2,fs)[3]
n = Kaiser(delta_1,delta_2,fs)[4]
tw = Kaiser(delta_1,delta_2,fs)[5]

w1 = np.kaiser(M,beta)

h_d = bp(n,M,ft1,ft2)    #   ideal impulse response
h = h_d * w
h1 = h_d[:len(w1)] * w1              #   windowed impulse response  

h = h[:len(h)-1]

########################################
plt.plot(n,h_d)
plt.show()

plt.plot(n,w)
#plt.axis([0,1500,-0.5,1.2])
plt.show()

H = np.abs(np.fft.fft(h,len(signal)))
H1 = np.abs(np.fft.fft(h1,len(signal)))

plt.plot(freq_ax,np.abs(H)[:sampels/2],'r')
plt.plot(freq_ax,np.abs(H1)[:sampels/2],'b')
plt.axis([0,100,0,1.1])
plt.xlabel("frequency [Hz]")
plt.show()
plt.plot(freq_ax,np.unwrap(np.angle(np.fft.fft(h1,len(signal))[:sampels/2])))
#plt.axis([0,100,-100,100])
plt.show()

# Pure signal in frequency domain
SIGNAL = 2/float(len(signal))*np.abs(np.fft.fft(signal)) #normaliseret
#plt.plot(freq_ax[freq_inter1:freq_inter2],SIGNAL[freq_inter1:freq_inter2])
#plt.show()


DATA = 2/float(len(data))*np.abs(np.fft.fft(data))


#==============================================================================
# Filtereing 
#==============================================================================
filt_DATA = H * DATA

plt.plot(freq_ax[freq_inter1:freq_inter2],DATA[freq_inter1:freq_inter2], label='Data')
plt.show()

plt.plot(freq_ax[freq_inter1:freq_inter2],filt_DATA[freq_inter1:freq_inter2], label='Filteret data')
plt.axis([0,4000,0,0.12])
plt.show()

#==============================================================================
# plots til rapport
#==============================================================================

#W = np.abs(np.fft.fft(np.pad(w,(0,sampels-N),'constant',constant_values=0)))
#W_dB= 20 * np.log10(np.abs(W))

#W = np.abs(np.fft.fft(w,len(signal)))
#W_dB= (20 * np.log10(np.abs(W)))-60
#
#omega = np.linspace(0,np.pi,len(freq_ax)/2)
#
#plt.plot(omega,W[:len(omega)])
##plt.axis([0,2000,0,1100])
#plt.show()
#
#plt.plot(omega,W_dB[:len(omega)])
#plt.axis([0,0.5,-100,1])





