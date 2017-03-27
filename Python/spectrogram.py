# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 09:40:07 2017

@author: Frederik Vardinghus
"""

from __future__ import division
import scipy.io.wavfile as siw
import scipy.signal as ss
import matplotlib.pyplot as plt
import numpy as np

#==============================================================================
# Hent data
#==============================================================================
#wav = siw.read('4chords2.wav')
#freq = wav[0]
#data = wav[1]
#
#data1 = np.zeros(len(data))
#data2 = np.zeros(len(data))
#for i in range(len(data)):
#    data1[i] = data[i][0]
#    data2[i] = data[i][1]
#    
#time = len(data)/freq
#lin = np.linspace(0,time,len(data))

#plt.plot(lin,data1)

#==============================================================================
# Generér egne data
#==============================================================================
freq = 1000
N = 10000
time = N/freq
lin = np.linspace(0,time,N)

def f(x):
    return np.sin(2*np.pi*(x**3))+np.sin(100*np.pi*x)

data = np.zeros(N)
for i in range(N):
    data[i] = f(lin[i])

#plt.plot(lin,data)

#==============================================================================
# Udtag og plot afsnit af data
#==============================================================================
#N = 2000 # Længde af udsnit
#data = data[0:N]
#lin = lin[0:N]
#plt.plot(lin,data)

#==============================================================================
# Foretag FFT
#==============================================================================
#F = np.fft.fft(data)

#==============================================================================
# Udtag og plot afsnit af F
#==============================================================================
#NF = 22050 # Længde af udsnit
#dataF = np.fft.fft(data[0:NF])
#linF = np.linspace(0,NF-1,NF)
    
#plt.plot(linF,dataF)


#==============================================================================
# Spectrogram genereret med scipy.signal
#==============================================================================
plt.specgram(data,Fs=freq)

#plt.specgram(data1,Fs=freq)
#plt.specgram(data2,Fs=freq)












