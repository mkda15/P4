# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 09:40:07 2017

@author: Frederik Vardinghus
"""

from __future__ import division
import scipy.io.wavfile as siw
import matplotlib.pyplot as plt
import numpy as np
import time
#import winsound
plt.style.use('ggplot')
#==============================================================================
# Hent data
#==============================================================================
wav1 = siw.read('Lydfiler/forsoeg_nopeak/enkelt_tone/forsoeg_enkelt_dyb.wav')
#wav1 = siw.read('Lydfiler/forsoeg_nopeak/enkelt_tone/forsoeg_enkelt_lys.wav')
#wav2 = siw.read('Lydfiler/forsoeg_nopeak/stoej/kroelle_stoej.wav')
wav2 = siw.read('Lydfiler/forsoeg_nopeak/stoej/klap_takt_2.wav')
#wav2 = siw.read('Lydfiler/forsoeg_nopeak/stoej/tale_1.wav')
#wav2 = siw.read('Lydfiler/forsoeg_nopeak/stoej/klap_random_huj_1.wav')
#wav = siw.read('Lydfiler/forsoeg_nopeak/melodi/alene/forsoeg_lillepeteredderkop_langsom.wav')
#wav1 = siw.read('Lydfiler/forsoeg_nopeak/melodi/akkorder/forsoeg_lillepeteredderkop_langsom.wav')
#wav1 = siw.read('Lydfiler/forsoeg_nopeak/akkorder/forsoeg_akkord_lys.wav')
#wav1 = siw.read('Lydfiler/forsoeg_nopeak/skala/forsoeg_skala_hurtig.wav')
#wav2 = siw.read('Lydfiler/forsoeg_nopeak/stoej/sang_3_mester_jakob.wav')
freq1 = wav1[0] #signal
data1 = wav1[1]
freq2 = wav2[0] #noise
data2 = wav2[1]

# forlængelse af data, ligger nul i enden så de er lige lange #
if len(data2) < len(data1):
    data2 = np.lib.pad(data2,(0,len(data1)-len(data2)),'constant',constant_values=0)
else:
    data2 = data2[0:len(data1)] # forkorter støjen
    #data1 = np.lib.pad(data1,(0,len(data2)-len(data1)),'constant',constant_values=0)

freq = freq1
data1 = data1/1000  #normalized data
data2 = (data2/1000)*10
data = data1
N = len(data)
if N % 2 > 0:
    raise ValueError('Brug N = potenser af 2')

plot = 1  # ændre til 2 for plot af spectrogram
time_inter1 = 0
time_inter2 = len(data1)
freq_inter1 = 100
freq_inter2 = 20000

save = 0


#data1 = np.zeros(len(data))
#data2 = np.zeros(len(data))
#for i in range(len(data)):
#    data1[i] = data[i][0]
#    data2[i] = data[i][1]
    
t = len(data)/freq
lin = np.linspace(0,t,len(data))

bins = np.linspace(0,freq/2,len(data)/2)
#plt.plot(lin,data1)

#==============================================================================
# dB funktion
#==============================================================================
def db(x):
    return 20*np.log10(x)


#==============================================================================
# Funktion som foretager FFT, generere Frekvensplot og spektrogram med scipy.signal
#==============================================================================

def Freq_plot(lin,data):
    
    start = time.time()
    F = 2/float(len(data))*np.abs(np.fft.fft(data)[:len(data)/2])
    end = time.time()
    FFT_time = end - start

    if plot == 0:
        plt.plot(lin[time_inter1:time_inter2],data[time_inter1:time_inter2])
        plt.legend('Signal')
        plt.ylabel('Amplitude')
        plt.xlabel('Time (sec)')
    elif plot == 1:
        plt.plot(bins[freq_inter1:freq_inter2],F[freq_inter1:freq_inter2])
        plt.legend('Fourier transform of signal')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude')
    elif plot == 2:
        plt.specgram(data,Fs=freq)
        plt.xlabel('Time (sec)')
        plt.ylabel('Frequency (Hz)')

    return plt.show(),bins[np.argmax(F)], FFT_time, len(data), freq

#if save == 1:
#    siw.write('Lydfiler/forsoeg_nopeak/melodi/akkorder/forsoeg_lillepeteredderkop_langsom.wav',freq,data)
#winsound.PlaySound('Lydfiler/forsoeg_nopeak/skala/forsoeg_skala_hurtig.wav', winsound.SND_FILENAME)

#==============================================================================
# Udtag og plot afsnit af F (frekvens)
#==============================================================================
#dataF = np.fft.fft(data1)
#linF = np.linspace(0,len(data1)-1,len(data1))
#
##startF = 1200 # Længde af udsnit
##slutF = startF+2000
##dataF = np.fft.fft(data2[startF:slutF])
##linF = linF[startF:slutF]
#
#plt.plot(linF,dataF)
#plt.show()


#==============================================================================
# Tilføjer støj til data 
#==============================================================================
def add_noise(data,noise,c = 0.5): #kilde side 229 i DTSP
#    signal = np.convolve(data,noise)
#    return signal
    signal=np.zeros(len(data))    
    for i in range(len(data)):
        signal[i]=data[i]+float(c)*noise[i]
    return signal

#==============================================================================
# Printer signal 
#==============================================================================
## Lav udsnit af signal   
#start = 100000 # start udsnit
#slut = start+1000 # slut udsnit
#data1 = data1[start:slut]
#data2 = data2[start:slut]
#lin = lin[start:slut]

#signal = add_noise(data1,data2)
#
#plt.plot(lin,signal, label='signal')
#plt.plot(lin,data2, label='noise')
#plt.plot(lin,data1, label='data')
#plt.xlabel('Time[Seconds]')
#plt.ylabel('Amplitude [Voltages]')
#plt.legend(loc= "lower right")
#plt.show()

##==============================================================================
## Printer frekvens plot af ren data
##==============================================================================
freq_plot = Freq_plot(lin, data)
print 'Mest betydende frekvens:',freq_plot[1]
print 'Beregningstid for numpy.fft:',freq_plot[2]
print 'Længde af data:',freq_plot[3]
print 'Samplingsfrkvens:',freq_plot[4]

#==============================================================================
# Printer frekvensplot af data med støj
#==============================================================================
#freq_plot2 = Freq_plot(lin, signal)
#print 'Mest betydende frekvens:',freq_plot2[1]
#print 'Beregningstid for numpy.fft:',freq_plot2[2]
#print 'Længde af data:',freq_plot2[3]
#print 'Samplingsfrkvens:',freq_plot2[4]

#==============================================================================
# Printer frekvensplot af ren støj
#==============================================================================
#freq_plot3 = Freq_plot(lin, data2)
#print 'Mest betydende frekvens:',freq_plot3[1]
#print 'Beregningstid for numpy.fft:',freq_plot3[2]
#print 'Længde af data:',freq_plot3[3]
#print 'Samplingsfrkvens:',freq_plot3[4]

#==============================================================================
# Signal to noise ratio
#==============================================================================
#def power(data):
#    N=len(data)    
#    a=0
#    for i in range(N):
#        a = np.abs(data[i])**2
#    return (1/N)*a
#
#P_n = power(data2)
#P_sn = power(signal)*1000
#
#def snr(Pn,Ps):
#    return 10*np.log10((Ps-Pn)/Pn)
#
#SNR = snr(P_n,P_sn)
   