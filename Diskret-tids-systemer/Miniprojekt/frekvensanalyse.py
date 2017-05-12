# -*- coding: utf-8 -*-
"""
Created on Thu May 11 22:25:23 2017

@author: Frederik Vardinghus
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import matplotlib.patches as mpatches

N = np.array([2**8,2**12,2**18])
fs = np.array([1,32])
fontsize = 13

#==============================================================================
# Funktioner
#==============================================================================

def sig(x):
    return np.sin(np.pi/3.*x) + np.sin(np.pi/2.*x+2*np.pi/3.) + np.sin(3*np.pi/4.*x+4*np.pi/3.)

def fft(x):
    return np.fft.fft(x)

def F(x):
    return np.abs(fft(x))/np.max(np.abs(fft(x)))

#==============================================================================
# Plot af signal i tid og frekvens ved forskellige
#   samplingsfrekvenser og antal samples
#==============================================================================

fig = 0

for f in range(len(fs)):
    for n in range(len(N)):
        plt.figure(fig)
        t = np.linspace(0,N[n]/float(fs[f]),N[n])
        bins = 2*np.pi*np.linspace(0,fs[f]/2.,N[n]/2.)
        S = F(sig(t))
        if f == 0:
            plt.plot(bins,S[:len(S)/2])
            plt.title('Samples N = %s og samplingsfrekvens fs = %s' %(N[n],fs[f]),fontsize=fontsize)
            plt.xlabel('Frekvens (rad/s)',fontsize=fontsize)
            plt.ylabel('Amplitude',fontsize=fontsize)
            fig += 1
            if N[n]==N[2]:
                plt.figure(fig)
                plt.plot(t[:32],sig(t)[:32])
                plt.title('Signal samplet ved fs = 1 Hz',fontsize=fontsize)
                plt.xlabel('Tid (sek)',fontsize=fontsize)
                plt.ylabel('Amplitude',fontsize=fontsize)
                fig += 1
        else:
            plt.plot(bins[:len(bins)/32.],S[:len(S)/64.])
            plt.title('Samples N = %s og samplingsfrekvens fs = %s' %(N[n],fs[f]),fontsize=fontsize)
            plt.xlabel('Frekvens (rad/s)',fontsize=fontsize)
            plt.ylabel('Amplitude',fontsize=fontsize)
            fig += 1
            if N[n]==N[2]:
                plt.figure(fig)
                plt.plot(t[:1000],sig(t)[:1000])
                plt.title('Signal samplet ved fs = 32 Hz',fontsize=fontsize)
                plt.xlabel('Tid (sek)',fontsize=fontsize)
                plt.ylabel('Amplitude',fontsize=fontsize)
                fig += 1

#==============================================================================
# Filtrering
#==============================================================================

M = 92 # Orden af filter

t = np.arange(M) # Samplingstidspunkter

cut = np.pi/2 # Frekvens som ønskes frafiltreret
delta = np.pi/15 # Bredde af stopbåndet
cut1 = (cut - delta)/(2*np.pi) # Første knækfrekvens
cut2 = (cut + delta)/(2*np.pi) # Anden knækfrekvens

def bs(n,M,ft1,ft2):
    ''' Båndstopfilter af typen FIR I og af orden M
        med knækfrekvenser ft1 og ft2'''
    h = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] == M/2:
            h[i] = 1-2*(ft2 - ft1)
        else:
            h[i] = (1 / (np.pi*(n[i] - M/2.)))*(np.sin(ft1*2*np.pi*(n[i] - M/2.)) \
            - (np.sin(ft2*2*np.pi*(n[i] - M/2.))))
    return h

def ha(n,M,a):
    ''' Hann window, hvis a = 0.5
        Hamming window, hvis a = 0.54 '''
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] >= 0 and n[i] <= M:
            w[i] = a - (1 - a)*np.cos((2*np.pi*n[i])/M)
        else:
            w[i] = 0
    return w

w = ha(t,M,0.54) # Vindue
h = bs(t,M,cut1,cut2) # Båndstopfilter
h = h*w # Windowing af båndstopfilter

bins = np.linspace(0,M/2,M/2)
s = sig(t) # Sampling af signal
s_ideal = s - np.sin(t*np.pi/2+2*np.pi/3)
       
sf = np.convolve(s,h) # Filtrering af signal
SF = F(s)*F(h)
                
plt.figure(fig)
plt.plot(bins,SF[:len(SF)/2])
plt.title('Frekvensspektrum af filtreret signal',fontsize=fontsize)
plt.xlabel('Frekvens (rad/s)',fontsize=fontsize)
plt.ylabel('Amplitude',fontsize=fontsize)
fig +=1 
plt.figure(fig)
plt.plot(t,s_ideal)
plt.plot(t,sf[M/2:len(sf)-M/2+1])
plt.xlabel('Tid (sek)',fontsize=fontsize)
plt.ylabel('Amplitude',fontsize=fontsize)
plt.title('Filtrerede og ideelle signal',fontsize=fontsize)




















        
        



