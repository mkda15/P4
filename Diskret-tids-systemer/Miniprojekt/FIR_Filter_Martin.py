# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 21:50:53 2017

@author: Martin Kamp Dalgaard
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

#==============================================================================
# Plot af den ideelle impulsrespons
#==============================================================================

t = np.linspace(0,np.pi,1000)
y = np.zeros(len(t))

delta = np.pi/15.
o1 = np.pi/2 - delta
o2 = np.pi/2 + delta

for i in range(len(t)):
    if o1 >= t[i]:
        y[i] = 1

for i in range(len(t)):
    if t[i] >= o2:
        y[i] = 1

plt.plot(t,y)
plt.axis([0,np.pi,0,2])
plt.axvline(o1, color='yellow')
plt.axvline(o2, color='yellow')
plt.axvline(np.pi/2, color='red')
plt.axvline(np.pi/3, color='green')
plt.axvline(3*np.pi/4, color='green')
plt.title('Den ideelle amplituderespons for filteret')
plt.xlabel('Frekvens [rad / s]')
plt.ylabel('Amplitude')

#==============================================================================
# Definitioner (filterorden, -længde, indeks, knækfrekvenser og impulsrespons)
#==============================================================================

M = 92
l = M+1

n = np.linspace(0,l,l+1)
x = np.linspace(-np.pi,np.pi,len(n))

def hd(n,M,f1,f2): 
    hd = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] == M/2:
            hd[i] = 1 - (o2 - o1)/np.pi
        else:
            hd[i] = (np.sin(o1*(n[i] - M/2.)) / (np.pi*(n[i] - M/2.))) \
            - (np.sin(o2*(n[i] - M/2.)) / (np.pi*(n[i] - M/2.)))
    return hd

#==============================================================================
# Vinduer
#==============================================================================

def rect(n,M): # Det ektangulaere vindue
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
# Beregning af impuls- og amplituderespons
#==============================================================================

w = ha(n,M,0.54)
hd = hd(n,M,o1,o2)

def fft(x,n):
    return np.fft.fft(x)

h = hd * w


H = np.abs(fft(h,n))
plt.figure(2)
plt.plot(20*np.log10(np.abs(np.fft.fft(np.pad(w,(0,1000),'constant',constant_values=0))))[:400])
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())

#==============================================================================
# Graf over amplituderesponsen
#==============================================================================

plt.figure(3)
plt.plot(x, H)
plt.axis([0,np.pi,0,2])
plt.title('Amplituderespons for filteret, Hamming-vinduet, $M = %d$' %(M))
plt.xlabel('Frekvens [rad / s]')
plt.ylabel('Amplitude')
plt.axvline(o1, color='yellow') # Nedre knækfrekvens
plt.axvline(o2, color='yellow') # Øvre knækfrekvens
plt.axvline(np.pi/2, color='red') # Frekvens, der skal elimineres
plt.axvline(np.pi/3, color='green') # Frekvens, der skal beholdes
plt.axvline(3*np.pi/4, color='green') # Frekvens, der skal beholdes

#==============================================================================
# Filteret i Scipy til sammenligning 
#==============================================================================

N = [o1,o2]
plt.figure(4)
b, a = signal.butter(10, N, 'bandstop', analog=True)
w, h = signal.freqs(b, a)
plt.plot(w, abs(h), "b-")
plt.title('Amplitude for Scipys baandstop-filter')
plt.xlabel('Frekvens [rad / s]')
plt.ylabel('Amplitude')
plt.margins(0, 0.1)
plt.axis([0,np.pi,0,2])
plt.axvline(o1, color='yellow') # Nedre knækfrekvens
plt.axvline(o2, color='yellow') # Øvre knækfrekvens
plt.axvline(np.pi/2, color='red') # Frekvens, der skal elimineres
plt.axvline(np.pi/3, color='green') # Frekvens, der skal beholdes
plt.axvline(3*np.pi/4, color='green') # Frekvens, der skal beholdes
plt.show()