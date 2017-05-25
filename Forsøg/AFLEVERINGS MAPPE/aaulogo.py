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

freq , signal = siw.read('Lydfiler/forsoeg_nopeak/enkelt_tone/forsoeg_enkelt_dyb.wav')  # Data signal
freq2, noise = siw.read('Lydfiler/forsoeg_nopeak/stoej/kroelle_stoej.wav')                  # Noise signal
freq3, noise2 = siw.read('Lydfiler/forsoeg_nopeak/stoej/klap_takt_2.wav')
freq4, noise3 = siw.read('Lydfiler/forsoeg_nopeak/stoej/sang_2_lille_peter_edderkop_2.wav')
freq5, noise4 = siw.read('Lydfiler/forsoeg_nopeak/stoej/stoej_fra_omraadet_rent.wav')

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
                      
                         
""" Length of data and noise alings"""


if len(signal) > len(noise):
    while len(signal) > len(noise):
        noise = np.append(noise,noise)
if len(signal) < len(noise):
        noise = noise[:len(signal)]

""" Variabler til filter """
window = Kaiser     # The wanted window is named (Has to be capitalised and has to be imported under windowfunctions)
M    = 1000.        # Filter order
cut  = 1000./freq   # Cut off frequency
cut1 = 70./freq     # Cut off frequency for band
cut2 = 600./freq    # Cut off frequency for band
sampels = len(signal) # Amount of sampels in the signal (data points)
plotlength = int(sampels/2) # Length for plotting (arbitrary)


""" Til Kaiser vinduet """
delta_1 = 0.05 # peak approximation error in amplitude
delta_2 = 10. # max transition width is 2*delta_2

""" Aksis og linspaces """
t   = sampels/float(freq)                   # The time for howlong the system runs (for making the time axis)
tid = np.linspace(0,t,sampels)              # Axis for time domain
freq_axis = np.linspace(0,freq/2,sampels/2) # Axis for frequency domain
freq_axis_norm = np.linspace(0,1,sampels/2)
n   = np.linspace(0,M,M+1)                  # Integer numbers for making window and impulsrespons

""" Variabler til spektogram """
freq_inter1 = 0
freq_inter2 = 10000
fftsize = 2**12
overlap = 2
beta = 4

fontsize = 13
dataType = "Tabs" #Variable to peak detection, if the file is with chords dataType == Chords if its tabs dataType should be == Tabs



#print("variabler og data importeret 1/9")

def noise_variation(signal,noise,gain):
    noise = noise*gain
    return impuls.add_noise(signal,noise,1), noise

#signal, noise = noise_variation(data,noise,0)
                         
#print('støj adderet 2/9')

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

#print('vindue generet 3/9')

hd = impuls.ImpulsresponsBP(n,M,cut1,cut2)
#hd = impuls.ImpulsresponsHP(n,M,cut)
#hd = impuls.ImpulsresponsLP(n,M,cut)    # Desired impulsrespons

h = hd * w                              # The final impulsrespons
H = np.fft.fft(h,(len(signal)))         # The fourier transformed of the final impulsrespons zero padded to fit the signal

#print('impuls respons udregnet 4/9')

def filtrering(signal,data,noise,H):
    signal = signal / float((np.max(signal))) # Reduktion of amplitude
    
    """ Dataen fourier transformeres """
    DATA = np.fft.fft(data)     # Pure signal in fourier
    NOISE = np.fft.fft(noise)   # Noise in fourier
    SIGNAL = np.fft.fft(signal) # Signal with noise in fourier
    
    #print('Data fourier transformeret 5/9')
    
    SIGNAL_FILT = SIGNAL*H                # Convolution between the filter H and the noise SIGNALx
    signal_filt = np.fft.ifft(SIGNAL_FILT)  # Filtered data
    signal_filt = np.real(signal_filt)      # Cast to real, to remove the 0j from ifft.
    
    #print('Data filtreret 6/9')
    
    return signal_filt

#signal_filt = filtrering(signal,data,noise,H)
#signal_filt = signal

#==============================================================================
# Plt plots af alt det intresante og data gemmes
#==============================================================================

#plt.plot(freq_axis[:plotlength-50000],np.abs(DATA)[:plotlength-50000])          # FFT of clean data
#plt.xlabel('Ren Data freq')
#plt.show()
#
#plt.plot(freq_axis[:plotlength],np.abs(H)[:plotlength])             # FFT of impuls respons
#plt.xlabel('Ren stoej freq')
#plt.show()
#
#plt.plot(freq_axis[:2000],np.abs(NOISE)[:2000])         # FFT of the noise
#plt.xlabel('Stoej og data freq')
#plt.show()
#
#plt.plot(tid,data)                                                  # Original data
#plt.xlabel('Ren Data tid')
#plt.show()
#
#plt.plot(tid,noise)                                                 # Original noise
#plt.xlabel('Ren Stoej tid')
#plt.show()
#
#plt.plot(tid,signal)                                                # Original data with noise added
#plt.xlabel('Stoej og data tid')
#plt.show()
#
#plt.plot(tid,signal_filt)                                           # The filtered data
#plt.xlabel('Filtreret signal')
#plt.show()
#
#plt.plot(freq_axis[:plotlength],np.abs(SIGNAL_FILT[:plotlength]))   # FFT of the filtered data
#plt.xlabel('Freq af Filtreret signal')
#plt.show()

#print('plot plotteret 7/9')


""" Data gemmes """
#siw.write('Lydfiler/forsoeg_nopeak/output/out_signal_filt.wav',freq,signal_filt)    # The filtered data is saved
#siw.write('Lydfiler/forsoeg_nopeak/output/out_data.wav',freq,data)                  # Original noise is saved with same length as data
#siw.write('Lydfiler/forsoeg_nopeak/output/out_noise.wav',freq,noise)                # Original data is saved with same length as noise
#siw.write('Lydfiler/forsoeg_nopeak/output/out_signal.wav',freq,signal)              # The signal with noise is saved

#print('Data gemt 8/9')

#==============================================================================
# Spectrogram
#==============================================================================

def stft_plot(signal_filt,fftsize,overlap,freq_inter1,freq_inter2,beta):
    X = stft(signal_filt,fftsize,overlap,beta)     # STFT calculated

    #print('stft udregnet 9/9')
    
    X = db(np.abs(X).T)                             # Calculated to dB
    
    x = np.linspace(0,tid[-1],np.shape(X)[1])
    y = np.linspace(0,freq_axis[-1],np.shape(X)[0])
    
    
#    spec = plt.pcolormesh(x,y[freq_inter1:freq_inter2],X[freq_inter1:freq_inter2],cmap='hot')
#    cb   = plt.colorbar(spec)
#    cb.set_label(label = 'Amplitude (dB)', fontsize=fontsize)
#    plt.xlabel('Time (sec)', fontsize = fontsize)
#    plt.ylabel('Frequency (Hz)', fontsize = fontsize)
#    plt.show()
    
    return X,x,y

#X,x,y = stft_plot(signal_filt,fftsize,overlap,freq_inter1,freq_inter2)


#plt.plot(freq_axis,np.angle(H)[:sampels/2])
#plt.axis([0,500,-4,4])
#plt.show()
#
#Hdb =  db(np.abs(H).T)
#
#plt.plot(freq_axis,Hdb[:sampels/2])
#plt.axis([0,500,-100,2])
#plt.show()


#==============================================================================
# Peak analysis
#==============================================================================

def sort(X):
    X = X.T
    # kan laves til en definition og placeres i et andet dokument.
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
#    print 'Hyppigste frekvens', u[np.argmax(np.bincount(indices))]
    return u[np.argmax(np.bincount(indices))]


#if dataType == "Tabs": #Tjeck if data is in single tabs or chords
#    tabs(X,x)

#elif dataType == "Chords":
#    max_freq_pos1 = np.zeros(len(X))
#    max_freq_pos2 = np.zeros(len(X))
#    max_freq_pos3 = np.zeros(len(X))
#
#    for i in range(len(X)):
#        if sortedX[i][-1] > 20:
#            a = np.where(X[i][:] == sortedX[i][-1])
#        else:
#            a = [[0]]
#        if sortedX[i][-2] > 20:
#            b = np.where(X[i][:] == sortedX[i][-2])
#        else:
#            b = [[0]]
#        if sortedX[i][-3] > 20:
#            c = np.where(X[i][:] == sortedX[i][-3])
#        else:
#            c = [[0]]
#        max_freq_pos1[i] = a[0][0]
#        max_freq_pos2[i] = b[0][0]
#        max_freq_pos3[i] = c[0][0]
#
#    max_freq_t1 = np.zeros(len(X))
#    max_freq_t2 = np.zeros(len(X))
#    max_freq_t3 = np.zeros(len(X))
#    for i in range(len(X)):
#        if max_freq_pos1[i] == 0:
#            max_freq_t1[i] = 0
#        else:
#            max_freq_t1[i] = y[int(max_freq_pos1[i])]
#        if max_freq_pos2[i] == 0:
#            max_freq_t2[i] = 0
#        else:
#            max_freq_t2[i] = y[int(max_freq_pos2[i])]
#        if max_freq_pos3[i] == 0:
#            max_freq_t3[i] = 0
#        else:
#            max_freq_t3[i] = y[int(max_freq_pos3[i])]
#    plt.plot(max_freq_t1)
#    plt.plot(max_freq_t2)
#    plt.plot(max_freq_t3)
#
#    sted = 75
#    print(max_freq_t1[sted])
#    print(max_freq_t2[sted])
#    print(max_freq_t3[sted])

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

N = 50 # Antal loops

SNR_a = np.zeros(N)
max_freq = np.zeros(N)
#noise = np.array([np.random.randn() for i in range(len(signal))])
#noise = np.array([np.sin(2*np.pi*300*tid[i]) for i in range(len(signal))])

for i in range(N):
#    noise = np.array([np.random.randn() for k in range(len(signal))])
    gain = (i+1)*0.1
#    beta = i*0.5
    signal_new, noise_new = noise_variation(signal,noise,gain)
    signal_filt = filtrering(signal_new,signal,noise,H)
#    signal_filt = signal_new
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
        print 'Ved SNR = %s detekteres den rigtige frekvens ikke længere' % SNR_a[i+1]


inter = 10000

S = np.abs(np.fft.fft(signal))
N = np.abs(np.fft.fft(noise))

#plt.plot(freq_axis[:inter],S[:inter]/np.max(S[:inter]))
#plt.plot(freq_axis[:inter],np.abs(H)[:inter])
#plt.plot(freq_axis[:inter],N[:inter]/np.max(N[:inter]))
#plt.xlabel('Frequency [Hz]',fontsize=13)
#plt.ylabel('Amplitude',fontsize=13)

















