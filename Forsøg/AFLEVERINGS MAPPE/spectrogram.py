# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 12:51:39 2017

@author: Frederik Vardinghus
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import matplotlib.colors as color
import scipy.io.wavfile as siw

#==============================================================================
# Load data
#==============================================================================
wav = siw.read('Lydfiler/skala/forsoeg_skala_hurtig.wav')
freq = wav[0]
data = wav[1][:-150000]

t = len(data)/float(freq)
tid = np.linspace(0,t,len(data))
bins = np.linspace(0,freq/2.,len(data)/2.)
N = len(data)

#==============================================================================
# Parametre
#==============================================================================
#N = 2**18 # Datalængde
#freq = 2**12 # Samplingsfrekvens i Hz
#td = 1/float(freq) # Samplingsperiode
#t = td*N # Samplingsinterval

M = 100 # Længde af vindue

fft_size = 2**12
overlap = 2
beta = 0

save = 0

#==============================================================================
# Genereret signal
#==============================================================================
def sig(x):
    if x < t/3.:
        return np.sin(500*2*np.pi*x)
    elif t/3. < x < 2*t/3.:
        return np.sin(1500*2*np.pi*x)
    else:
        return np.sin(1500*2*np.pi*x)+np.sin(1625*2*np.pi*x)

#def sig(x):
#    if x < t/12.:
#        return 0
#    elif t/12 <= x < t/5.:
#        return np.sin(x**2)
#    else:
#        return np.sin(np.sqrt(x))

#def sig(x):
#    return np.sin(500*2*np.pi*x)

#start = time.time()
#S = np.fft.fft(s)
#fft_time = time.time() - start

#==============================================================================
# Linspaces
#==============================================================================
#tid = np.linspace(0,t,N) # Samplingstidspunkter
#bins = np.linspace(0,freq/2.,N/2.) # Frekvensintervaller
#data = np.array([sig(i) for i in tid]) # Signal samplet i tid
               
#==============================================================================
# Functions                 
#==============================================================================
def db(x):
    return 20*np.log10(x)           

#==============================================================================
# STFT
#==============================================================================
def stft(x, fftsize, overlap, beta):  
    if beta == 0:    
        hop = fftsize / overlap
        w = sc.blackman(fftsize+1)[:-1]      # better reconstruction with this trick +1)[:-1]  
        X = np.array([np.fft.rfft(w*x[i:i+fftsize]) for i in range(0, len(x)-fftsize, hop)])
        y = np.linspace(0,bins[-1],np.shape(X)[1])
        x = np.linspace(0,tid[-1],np.shape(X)[0])
        return X,x,y
    if beta != 0:    
        hop = fftsize / overlap
        w = sc.hanning(fftsize+1)[:-1]      # better reconstruction with this trick +1)[:-1]  
        X = np.array([np.fft.rfft(w*x[i:i+fftsize]) for i in range(0, len(x)-fftsize, hop)])
        y = np.linspace(0,bins[-1],np.shape(X)[1])
        x = np.linspace(0,tid[-1],np.shape(X)[0])
        return X,x,y

X,x,y = stft(data,fft_size,overlap,beta)

#==============================================================================
# Peak detection
#==============================================================================
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
            
    plt.plot(x,max_freq_t)
    plt.axis([0,x[-1],0,np.max(max_freq_t)*1.1])
    plt.xlabel('Time [s]', fontsize = 13)
    plt.ylabel('Frequency [Hz]', fontsize = 13)
    plt.show()
    u,indices = np.unique(max_freq_t,return_inverse=True)
#    print 'Hyppigste frekvens', u[np.argmax(np.bincount(indices))]
    return u[np.argmax(np.bincount(indices))],max_freq_t

#==============================================================================
# Spectrogram
#==============================================================================

#a,b = tabs(X,x)

X = db(np.abs(X).T)
freq_inter1 = 0
freq_inter2 = 450
fontsize = 13
vmin = np.min(X)
vmax = np.max(X)
spec = plt.pcolormesh(x,y[freq_inter1:freq_inter2],X[freq_inter1:freq_inter2],cmap='jet',vmin=vmin,vmax=vmax)
cb = plt.colorbar(spec)
cb.set_label(label='Amplitude [dB]',fontsize=fontsize)
plt.xlabel('Time [s]',fontsize=fontsize)
plt.ylabel('Frequency [Hz]',fontsize=fontsize)

#print 'Max dB =',np.max(X)
#print 'Min dB =',np.min(X)

#if save == 1:
#    siw.write('forsoeg_nopeak/melodi/akkorder/forsoeg_lillepeteredderkop_langsom.wav',freq,data)




#==============================================================================
# Varying fftsize
#==============================================================================
#iters = 10
#vari = np.zeros(iters)
#max_freq = np.zeros(iters)

#for i in xrange(iters):
#    overlap = 2
#    fft_size = 2**12
#    beta = (i+1)*0.5
#    X,x,y = stft(data,fft_size,overlap,beta)
#    X = X.T
#
#    max_freq[i] = tabs(X,x)
#    vari[i] = beta
#
#    plt.stem(vari,max_freq)



#    freq_inter1 = 0
#    freq_inter2 = len(y)
#    
#    fontsize = 13
#    
#    spec = plt.pcolormesh(x,y[freq_inter1:freq_inter2],X[freq_inter1:freq_inter2],beta,cmap='jet')
#    cb = plt.colorbar(spec)
#    cb.set_label(label='Amplitude [dB]',fontsize=fontsize)
#    plt.xlabel('Time [s]',fontsize=fontsize)
#    plt.ylabel('Frequency [rad/s]',fontsize=fontsize)

#==============================================================================
# Varying window and beta
#==============================================================================













