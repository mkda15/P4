# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 13:59:42 2017

@author: Martin Kamp Dalgaard
"""

import numpy as np
import matplotlib.pyplot as plt

freq1 = 440
freq2 = 880

periods1 = 4
range1 = 4/440.

x = np.linspace(0.0, range1, 10000)
y_1 = np.zeros(len(x))
y_2 = np.zeros(len(x))
y_sum = np.zeros(len(x))

def sine(x):
    return np.sin(x)

pi = np.pi

for i in range(len(x)):
    y_1[i] = sine(freq1*2*pi*x[i])
    y_2[i] = sine(freq2*2*pi*x[i])
    y_sum[i] = sine(freq1*2*pi*x[i]) + sine(freq2*2*pi*x[i])

# Figure 2.2(a):
plt.figure(1)
plt.plot(x, y_1, "b-", label = r'$y_f(t) = \ \sin(880\pi\cdot t)$')
plt.plot(x, y_2, "r-", label = r'$y_h(t) = \ \sin(1760\pi\cdot t)$')
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")
plt.legend(loc = "upper right")

# Figure 2.2(b):
plt.figure(2)
plt.plot(x,y_sum, "g-", label = r'$y(t) = \ \sin(880\pi\cdot t) + \sin(1760\pi\cdot t)$')
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")
plt.legend(loc = "upper right")