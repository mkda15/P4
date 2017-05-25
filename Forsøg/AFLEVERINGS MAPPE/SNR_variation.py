# -*- coding: utf-8 -*-

#==============================================================================
# Imports
#==============================================================================

from foriertransform import stft , db
from windowfunctions import Hamming, Hann, Kaiser
import numpy as np
import impulsrespons as impuls
import matplotlib.pyplot as plt
import scipy.io.wavfile as siw

#==============================================================================
# Variable and data import
#==============================================================================

""" Data import """

freq , signal = siw.read('Lydfiler/enkelt_tone/forsoeg_enkelt_dyb.wav')  # Data signal
freq2, noise = siw.read('Lydfiler/stoej/kroelle_stoej.wav')              # Noise signal
freq3, noise2 = siw.read('Lydfiler/stoej/klap_takt_2.wav')
freq4, noise3 = siw.read('Lydfiler/stoej/sang_2_lille_peter_edderkop_2.wav')
freq5, noise4 = siw.read('Lydfiler/stoej/stoej_fra_omraadet_rent.wav')

if len(noise) > len(noise2):
    noise = noise[:len(noise2)]+noise2
else:
    noise = noise+noise2[:len(noise)]
    
if len(noise) > len(noise3):
    noise = noise[:len(noise3)]+noise3
else:
    noise = noise+noise3[:len(noise)]

if len(noise) > len(noise4):
    noise = noise[:len(noise4)]+noise
else:
    noise = noise+noise4[:len(noise)]

signal = np.array([signal[i] for i in range(len(signal))],dtype=np.int64)
noise = np.array([noise[i] for i in range(len(noise))],dtype=np.int64)
                      
                         
""" Length of data and noise aligned"""

if len(signal) > len(noise):
    while len(signal) > len(noise):
        noise = np.append(noise,noise)
if len(signal) < len(noise):
        noise = noise[:len(signal)]

""" Variables for filter """
window = Kaiser     # The desired window is named (has to be capitalised and has to be imported under windowfunctions)
M    = 1000.        # Filter order
cut  = 1000./freq   # Cut-off frequency
cut1 = 70./freq     # Lower cut-off frequency
cut2 = 600./freq    # Upper cut-off frequency
sampels = len(signal) # Amount of sampels in the signal (data points)
plotlength = int(sampels/2) # Length for plotting (arbitrary)

""" For the Kaiser window """
delta_1 = 0.05 # peak approximation error in amplitude
delta_2 = 10. # max transition width is 2*delta_2

""" Axis and linspaces """
t   = sampels/float(freq)                   # The time for howlong the system runs (for making the time axis)
tid = np.linspace(0,t,sampels)              # Axis for time domain
freq_axis = np.linspace(0,freq/2,sampels/2) # Axis for frequency domain
freq_axis_norm = np.linspace(0,1,sampels/2)
n  = np.linspace(0,M,M+1)                   # Integer numbers for making window and impulse response

""" Variables for spektogram """
freq_inter1 = 0
freq_inter2 = 10000
fftsize = 2**12
overlap = 2
beta = 4

fontsize = 13
dataType = "Tabs" #Variable to peak detection: if the file is with chords dataType == "Chords"; if its tabs the datatype should be == "Tabs"

def noise_variation(signal,noise,gain):
    noise = noise*gain
    return impuls.add_noise(signal,noise,1), noise

#signal, noise = noise_variation(data,noise,0)

#==============================================================================
# Filter coefficients calculated and filters applied
#==============================================================================

""" Window and calculation of impulse response """

if window == Kaiser:   # What to be returned for different windowfunctions
    w,M,beta,A,n = window(delta_1,delta_2,freq)
    w = w[:-1]
    n = n[:-1]
else:
    w = window(n,M)

hd = impuls.ImpulseresponseBP(n,M,cut1,cut2)
#hd = impuls.ImpulsresponsHP(n,M,cut)
#hd = impuls.ImpulsresponsLP(n,M,cut)   # Desired impulsrespons

h = hd * w                              # The final impulse response
H = np.fft.fft(h,(len(signal)))         # The fourier transformed of the final impulsrespons zero padded to fit the signal

def filtrering(signal,data,noise,H):
    signal = signal / float((np.max(signal))) # Reduction of amplitude
    
    """ Fourier transform of data """
    DATA = np.fft.fft(data)     # Fourier transform of pure signal
    NOISE = np.fft.fft(noise)   # Fourier transform of noise
    SIGNAL = np.fft.fft(signal) # Fourier transform of signal with noise
    
    SIGNAL_FILT = SIGNAL*H                  # Convolution between the filter H and the noise SIGNAL
    signal_filt = np.fft.ifft(SIGNAL_FILT)  # Filtered data
    signal_filt = np.real(signal_filt)      # Cast to real, to remove the 0j from ifft.
    
    return signal_filt

#==============================================================================
# Spectrogram
#==============================================================================

def stft_plot(signal_filt,fftsize,overlap,freq_inter1,freq_inter2,beta):
    X = stft(signal_filt,fftsize,overlap,beta)     # STFT calculated

    #print('stft udregnet 9/9')
    
    X = db(np.abs(X).T)                             # Calculated to dB
    
    x = np.linspace(0,tid[-1],np.shape(X)[1])
    y = np.linspace(0,freq_axis[-1],np.shape(X)[0])
    
    return X,x,y




#==============================================================================
# Peak analysis
#==============================================================================

def sort(X):
    X = X.T
    sortedX = np.zeros(len(X),dtype = object)
    for i in range(len(X)):
        sortedX[i] = np.sort(X[i])

def tabs(X,x):
    max_freq_pos = np.zeros(len(X))
    for i in range(len(X)):
        if np.max(X[i]) > -3:
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
            
#    plt.stem(x,max_freq_t)
#    plt.xlabel('Time (sec)', fontsize = fontsize)
#    plt.ylabel('Frequency (Hz)', fontsize = fontsize)
#    plt.show()
    u,indices = np.unique(max_freq_t,return_inverse=True)
    return u[np.argmax(np.bincount(indices))]

#==============================================================================
# SNR
#==============================================================================

def RMS(x):
    x = x**2
    RMS = np.sqrt((np.sum(x))) / (len(x))
    return RMS

def SNR(signal,noise):
    if RMS(noise) == 0:
        return 'Infinite'
    else:
        return 20*np.log10(((RMS(signal))**2 / (RMS(noise))**2))

#==============================================================================
# Increase of SNR
#==============================================================================

N = 50 # Number of loops

SNR_a = np.zeros(N)
max_freq = np.zeros(N)
#noise = np.array([np.random.randn() for i in range(len(signal))])
#noise = np.array([np.sin(2*np.pi*300*tid[i]) for i in range(len(signal))])

for i in range(N):
#    noise = np.array([np.random.randn() for k in range(len(signal))])
    gain = (i+1)*0.1
    signal_new, noise_new = noise_variation(signal,noise,gain)
    signal_filt = filtrering(signal_new,signal,noise,H)
    X,x,y = stft_plot(signal_filt,fftsize,overlap,freq_inter1,freq_inter2,beta)
    X = X.T
    max_freq[i] = tabs(X,x)
    
    SNR_a[i] = SNR(signal,noise_new)
    print i
    
plt.stem(SNR_a,max_freq)
plt.xlabel('SNR [dB]',fontsize=13)
plt.ylabel('Most significant frequency',fontsize=13)
for i in range(N-1):
    if max_freq[i+1]-max_freq[i] != 0:
        print 'At SNR = %s the right frequency is no longer detected.' % SNR_a[i+1]

inter = 10000

S = np.abs(np.fft.fft(signal))
N = np.abs(np.fft.fft(noise))
