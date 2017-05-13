# -*- coding: utf-8 -*-
"""
Created on Sat May 13 12:13:39 2017

@author: Trine
"""
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np




b, a = signal.butter(4, 100, 'low', analog=True)
w, h = signal.freqs(b, a)

x = np.zeros(len(h)) 
labels = ['ehj']
labeld = labels*5

plt.figure()
plt.plot(w, abs(h))
plt.title('Butterworth filter frequency response')
plt.xlabel('Frequency [radians / second]')
plt.ylabel('Amplitude')
plt.margins(0, 0.1)
plt.grid(which='both', axis='both')
plt.axvline(100, color='green') # cutoff frequency
plt.axhline(1/np.sqrt(2), color='blue')
axes = figure().add_subplot(111)
a=axes.get_xticks().tolist()
a[1]='change'
plt.axes.set_xticklabels(a)
plt.show()