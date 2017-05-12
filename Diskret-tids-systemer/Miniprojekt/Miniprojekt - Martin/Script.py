# -*- coding: utf-8 -*-
"""
Created on Thu May 11 09:42:06 2017

@author: Martin Kamp Dalgaard
"""

import numpy as np
import matplotlib.pyplot as plt

#==============================================================================
# Kapitel 2: Generering af signal
#==============================================================================

N = 2**15

x = np.linspace(0,11*np.pi,N)

def sig(x):
    return np.sin((np.pi/3.)*x) + np.sin((np.pi/2.)*x+2*np.pi/3.) + np.sin((3*np.pi/4.)*x+4*np.pi/3.)

plt.figure(1)
plt.plot(x,sig(x))
plt.xlabel("Tid [s]")
plt.ylabel("Amplitude")
plt.title("Signal")
plt.savefig("Figures/signal.png")

#==============================================================================
# Kapitel 2: Sampling af signal
#==============================================================================

rang = 11*np.pi
fs = 1

x2 = np.linspace(0,rang,fs*rang)
y2 = sig(x2)

plt.figure(2)
plt.plot(x2,y2)
plt.title(r"Signal, $f_s = \ %d$ Hz" %(fs))
plt.xlabel("Tid [s]")
plt.ylabel("Amplitude")
plt.savefig("Figures/signal_%dhz.png" %(fs))

#==============================================================================
# Kapitel 3: Algoritmer
#==============================================================================

def DFT(x,c):
    X = np.zeros(c,dtype=complex)
    for k in range(len(x)):
        a = 0+0*1j
        for n in range(c):
            a += x[n]*np.exp(-2*np.pi*1j*k*n/float(c))
            X[k] = a
    return X

def is_power2(num): # Checks if a number is a power of 2
	return num != 0 and ((num & (num - 1)) == 0)

def FFT(x):
    N_new = len(x)
    if is_power2(N_new) == False:
        raise ValueError("N skal vaere en potens af 2.")
    elif N_new == 2:
        return DFT(x,N_new) # Returnerer DFT naar data ikke kan deles mere op
    else:
        X_even = FFT(x[::2]) # Deler rekursivt input op - lige dele
        X_odd = FFT(x[1::2]) # Deler rekursivt input op - ulige dele
        factor = np.exp(-2j * np.pi * np.arange(N_new) / N_new) # Twiddlefaktor
        return np.concatenate([X_even + factor[:int(N_new / 2)] * X_odd,
                               X_even + factor[int(N_new / 2):] * X_odd])

#==============================================================================
# Kapitel 3: Frekvensanalyse
#==============================================================================

N = np.array([2**8,2**12,2**18])
fs = np.array([1,32])
fontsize = 13

fig = 3

def fft(x):
    return np.fft.fft(x)

def F(x):
    return np.abs(fft(x))/np.max(np.abs(fft(x)))

for f in range(len(fs)):
    for n in range(len(N)):
        plt.figure(fig)
        t = np.linspace(0,N[n]/float(fs[f]),N[n])
        bins = 2*np.pi*np.linspace(0,fs[f]/2.,N[n]/2.)
        S = F(sig(t))
        if f == 0:
            plt.figure(fig)
            plt.plot(bins,S[:len(S)/2])
            plt.title(r"N = %d, $f_s$ = %d Hz" %(N[n],fs[f]))
            plt.xlabel("Frekvens [rad / s]")
            plt.ylabel("Amplitude")
            plt.axis([0,np.pi,0,1])
            plt.savefig("Figures/freq_%dhz_N%d.png" %(fs[f],N[n]))
            fig += 1
        else:
            f = 1
            plt.figure(fig)
            plt.plot(bins[:len(bins)/16],S[:len(S)/32])
            plt.title(r"N = %d, $f_s$ = %d Hz" %(N[n],fs[f]))
            plt.xlabel("Frekvens [rad / s]")
            plt.ylabel("Amplitude")
            plt.axis([0,np.pi,0,1])
            plt.savefig("Figures/freq_%dhz_N%d.png" %(fs[f],N[n]))
            fig += 1
            if n == 2:
                plt.figure(fig)
                plt.plot(t[:1000],sig(t)[:1000])
                plt.title(r"Signal, $f_s$ = %d Hz" %(fs[f]))
                plt.xlabel("Tid [s]")
                plt.ylabel("Amplitude")
                plt.savefig("Figures/signal_%dhz.png" %(fs[f]))
                fig += 1

#==============================================================================
# Kapitel 4: Den ideelle amplituderespons
#==============================================================================

t = np.linspace(0,np.pi,10000)
Hd = np.zeros(len(t))

cut = np.pi/2 # Frekvens som ønskes frafiltreret
delta = np.pi/15 # Bredde af stopbåndet
cut1_rad = cut - delta # Første knækfrekvens
cut2_rad = cut + delta # Anden knækfrekvens

for i in range(len(t)):
    if cut1_rad >= t[i]:
        Hd[i] = 1

for i in range(len(t)):
    if t[i] >= cut2_rad:
        Hd[i] = 1

plt.figure(10)
plt.plot(t,Hd,label=r"$|H_d(e^{j\omega})|$")
plt.legend(loc = "upper left")
plt.axis([0,np.pi,0,2])
plt.axvline(cut1_rad, color="yellow") # Nedre knækfrekvens
plt.axvline(cut2_rad, color="yellow") # Øvre knækfrekvens
plt.axvline(np.pi/2, color="red") # Frekvens, der skal elimineres
plt.axvline(np.pi/3, color="green") # Frekvens, der skal beholdes
plt.axvline(3*np.pi/4, color="green") # Frekvens, der skal beholdes
plt.xlabel("Frekvens [rad / s]")
plt.ylabel("Amplitude")
plt.savefig("Figures/ideel_amp_respons.png")

#==============================================================================
# Kapitel 4: Design af filter
#==============================================================================

M = 50 # Orden
l = M+1 # Længde af filter

n = np.linspace(0,l,l+1)
x = np.linspace(-np.pi,np.pi,len(n))

def bs_rad(n,M,f1,f2): 
    h = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] == M/2:
            h[i] = 1 - (f2 - f1)/np.pi
        else:
            h[i] = (np.sin(f1*(n[i] - M/2.)) / (np.pi*(n[i] - M/2.))) \
            - (np.sin(f2*(n[i] - M/2.)) / (np.pi*(n[i] - M/2.)))
    return h

def rect(n,M): # Det rektangulaere vindue
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] >= 0 and n[i] <= M:
            w[i] = 1
        else:
            w[i] = 0
    return w, "det rektangulaere vindue", "rekt"

def ha(n,M,a): # Hann-vindue, hvis a = 0.5. Hamming-vindue, hvis a = 0.54.
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] >= 0 and n[i] <= M:
            w[i] = a - (1 - a)*np.cos((2*np.pi*n[i])/M)
        else:
            w[i] = 0
    if a == 0.5:
        return w, "Hann-vinduet", "Hann"
    else:
        return w, "Hamming-vinduet", "Hamming"

def blackman(n,M): # Blackman-vinduet
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] >= 0 and n[i] <= M:
            w[i] = 0.42 - 0.5*np.cos((2*np.pi*n[i])/M) + 0.8*np.cos((4*np.pi*n[i])/M)
        else:
            w[i] = 0
    return w, "Blackman-vinduet", "Blackman"

w, s, q = ha(n,M,0.54)

plt.figure(11)
plt.stem(n,w)
plt.title("%s, M = %d" %(s, M))
plt.axis([0,52,0,1.1])
plt.xlabel(r"$n$")
plt.ylabel(r"$w[n]$")
plt.savefig("Figures/%s_%d.png" %(q, M))

#==============================================================================
# Kapitel 4: Amplituderespons for filteret
#==============================================================================

M2 = 50 # Orden
l2 = M2+1 # Længde af filter

n2 = np.linspace(0,l2,l2+1)
x2 = np.linspace(-np.pi,np.pi,len(n2))

w2, s2, q2 = rect(n2,M2) # Det rektangulære vindue
hd2 = bs_rad(n2,M2,cut1_rad,cut2_rad) # Ideel impulsrespons

h2 = hd2 * w2 # Den faktiske impulsrespons

H2 = np.abs(np.fft.fft(h2))

plt.figure(12)
plt.plot(x2, H2)
plt.axis([0,np.pi,0,2])
plt.title("Amplituderespons, %s, M = %d" %(s2, M2))
plt.xlabel("Frekvens [rad / s]")
plt.ylabel("Amplitude")
plt.axvline(cut1_rad, color="yellow")
plt.axvline(cut2_rad, color="yellow")
plt.axvline(np.pi/2, color="red")
plt.axvline((np.pi/3), color="green")
plt.axvline(3*np.pi/4, color="green")
plt.savefig("Figures/Filter_%s_%d.png" %(q2,M2))

M3 = 100 # Orden
l3 = M3+1 # Længde af filter

n3 = np.linspace(0,l3,l3+1)
x3 = np.linspace(-np.pi,np.pi,len(n3))

w3, s3, q3 = rect(n3,M3) # Det rektangulære vindue
hd3 = bs_rad(n3,M3,cut1_rad,cut2_rad) # Ideel impulsrespons

h3 = hd3 * w3 # Den faktiske impulsrespons

H3 = np.abs(np.fft.fft(h3))

plt.figure(13)
plt.plot(x3, H3)
plt.axis([0,np.pi,0,2])
plt.title("Amplituderespons, %s, M = %d" %(s3, M3))
plt.xlabel("Frekvens [rad / s]")
plt.ylabel("Amplitude")
plt.axvline(cut1_rad, color="yellow")
plt.axvline(cut2_rad, color="yellow")
plt.axvline(np.pi/2, color="red")
plt.axvline((np.pi/3), color="green")
plt.axvline(3*np.pi/4, color="green")
plt.savefig("Figures/Filter_%s_%d.png" %(q3,M3))

M4 = 92 # Orden
l4 = M4+1 # Længde af filter

n4 = np.linspace(0,l4,l4+1)
x4 = np.linspace(-np.pi,np.pi,len(n4))

w4, s4, q4 = ha(n4,M4,0.54) # Hamming-vinduet
hd4 = bs_rad(n4,M4,cut1_rad,cut2_rad) # Ideel impulsrespons

h4 = hd4 * w4 # Den faktiske impulsrespons

H4 = np.abs(np.fft.fft(h4))

plt.figure(14)
plt.plot(x4, H4)
plt.axis([0,np.pi,0,2])
plt.title("Amplituderespons, %s, M = %d" %(s4, M4))
plt.xlabel("Frekvens [rad / s]")
plt.ylabel("Amplitude")
plt.axvline(cut1_rad, color="yellow")
plt.axvline(cut2_rad, color="yellow")
plt.axvline(np.pi/2, color="red")
plt.axvline((np.pi/3), color="green")
plt.axvline(3*np.pi/4, color="green")
plt.savefig("Figures/Filter_%s_%d.png" %(q4,M4))

#==============================================================================
# Filtrering af signal
#==============================================================================

def bs_hz(n,M,f1,f2): # Ideel impulsrespons for båndstopfilteret i hertz
    h = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] == M/2:
            h[i] = 1-2*(f2 - f1)
        else:
            h[i] = (1 / (np.pi*(n[i] - M/2.)))*(np.sin(f1*2*np.pi*(n[i] - M/2.)) \
            - (np.sin(f2*2*np.pi*(n[i] - M/2.))))
    return h

cut = np.pi/2 # Frekvens som ønskes frafiltreret
delta = np.pi/15 # Bredde af stopbåndet
cut1_hz = (cut - delta)/(2*np.pi) # Første knækfrekvens
cut2_hz = (cut + delta)/(2*np.pi) # Anden knækfrekvens

M = 92 # Orden af filter

t = np.arange(M) # Samplingstidspunkter

w, q, s = ha(t,M,0.54) # Vindue
hd = bs_hz(t,M,cut1_hz,cut2_hz) # Båndstopfilter
h = hd*w

bins = np.linspace(0,np.pi,M/2)
s = sig(t) # Sampling af signal
s_ideal = s - np.sin(t*np.pi/2+2*np.pi/3)

sf = np.convolve(s,h) # Filtrering af signal gennem foldning
SF = F(s)*F(h)

plt.figure(15)
plt.plot(bins,SF[:len(SF)/2])
plt.title("Frekvensspektrum af filtreret signal")
plt.xlabel("Frekvens [rad / s]")
plt.ylabel("Amplitude")
plt.savefig("Figures/freq_filt_signal.png")

plt.figure(16)
plt.plot(t,s_ideal,label = "Det ideelle signal")
plt.plot(t,sf[M/2:len(sf)-M/2+1], label = "Det filtrerede signal")
plt.legend(loc = "upper right")
plt.xlabel("Tid [s]")
plt.ylabel("Amplitude")
plt.savefig("Figures/signal_compare.png")