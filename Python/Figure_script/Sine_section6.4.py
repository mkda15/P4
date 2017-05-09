# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:23:37 2017

@author: Martin Kamp Dalgaard
"""

import matplotlib.pyplot as plt
import numpy as np

def sin(x):
    return np.sin(x)

pi = np.pi

# First signal
Fs1 = 5000.0  # sampling rate
Ts1 = 1.0/Fs1 # sampling interval
t1 = np.arange(0,1,Ts1) # time vector

freq1 = 110   # frequency of the signal
y1 = sin(2*pi*freq1*t1) + sin(2*pi*2*freq1*t1)

n1 = len(y1) # length of the signal
k1 = np.arange(n1)
T1 = n1/Fs1
frq1 = k1/T1 # two sides frequency range
frq1 = frq1[range(n1/2)] # one side frequency range

Y1 = np.fft.fft(y1)/n1 # fft computing and normalization
Y1 = Y1[range(n1/2)]

# Second signal
Fs2 = 2000.0  # sampling rate
Ts2 = 1.0/Fs2 # sampling interval
t2 = np.arange(0,1,Ts2) # time vector

freq2 = 440 # frequency of the signal
y2 = sin(2*pi*freq2*t2) + sin(2*pi*2*freq2*t2)

n2 = len(y2) # length of the signal
k2 = np.arange(n2)
T2 = n2/Fs2
frq2 = k2/T2 # two sides frequency range
frq2 = frq2[range(n2/2)] # one side frequency range

Y2 = np.fft.fft(y2)/n2 # fft computing and normalization
Y2 = Y2[range(n2/2)]

# Signals shifting in time
range1 = 2.0*3/(freq1+2*freq1) # 2 periods of the first sine wave
range2 = 7.0*3/(freq2+2*freq2) # 7 periods of the second sine wave
t_sum = np.linspace(0.0, range1 + range2, 10000)
y_sum = np.zeros(len(t_sum))

for i in range(len(t_sum)):
    if t_sum[i] <= range1:
        y_sum[i] = sin(freq1*2*pi*t_sum[i]) + sin(2*freq1*2*pi*t_sum[i])
    else:
        y_sum[i] = sin(freq2*2*pi*t_sum[i]) + sin(2*freq2*2*pi*t_sum[i])

# Graphs

# Figure 6.3(a):
plt.figure(1)
plt.plot(t_sum, y_sum, "g-", label = "Shifting between two tones")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")
plt.axis([0,t_sum[-1],-2.0,2.6])
plt.legend(loc = "upper right")

# Figure 6.3(b):
plt.figure(2)
fig, ax = plt.subplots(2, 1)
ax[0].plot(frq1,abs(Y1), "b-", label = "Frequency spectrum of the first tone")
ax[0].axis([80,250,0,0.8])
ax[0].legend(loc = "upper left")
ax[0].set_ylabel("Amplitude")
ax[1].plot(frq2,abs(Y2),"r-", label = "Frequency spectrum of the second tone")
ax[1].axis([320,1000,0,0.8])
ax[1].legend(loc = "upper left")
ax[1].set_xlabel("Frequency [Hz]")
ax[1].set_ylabel("Amplitude")