
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as axes
import scipy.fftpack


M = 92
N = 2**15
fs = 2**5
t = N/fs


#n = np.linspace(0,M,N)
#n = np.pad(np.arange(M),(0,N-M),'constant',constant_values=0)
n = np.arange(M)
#tid = np.linspace(0,t,N)

delta = (np.pi/15)

cut1 = ((np.pi/2.) - delta)/(2*np.pi)
cut2 = ((np.pi/2.) + delta)/(2*np.pi)

       
a = 1
       
def sig(x):
    return np.sin(a*np.pi/3.*x) + np.sin(a*np.pi/2.*x+2*np.pi/3.) + np.sin(a*3*np.pi/4.*x+4*np.pi/3.)

def bs(n,M,ft1,ft2): # Bandpass filter
    h = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] == M/2:
            h[i] = 1-2*(ft2 - ft1)
        else:
            h[i] = (1 / (np.pi*(n[i] - M/2.)))*(np.sin(ft1*2*np.pi*(n[i] - M/2.)) \
            - (np.sin(ft2*2*np.pi*(n[i] - M/2.))))
    return h

def ha(n,M,a): # Hann window, if a = 0.5. Hamming window, if a = 0.54.
    w = np.zeros(len(n))
    for i in range(len(n)):
        if n[i] >= 0 and n[i] <= M:
            w[i] = a - (1 - a)*np.cos((2*np.pi*n[i])/M)
        else:
            w[i] = 0
    return w

w = ha(n,M,0.5)
h = bs(n,M,cut1,cut2)

h = h*w
#hw = np.pad(h*w,(0,N-M),'constant',constant_values=0)

def db(x):
    return 20*np.log10(np.abs(x))
def abs(x):
    return np.abs(x)

omega = np.linspace(0,np.pi,len(n))

s = sig(n)
S = np.abs(np.fft.fft(s))

s_ideal = s - np.sin(np.pi/2.*n+2*np.pi/3.)

H = np.abs(np.fft.fft(h))/np.max(np.abs(np.fft.fft(h)))
SF = np.abs(np.fft.fft(h)*np.fft.fft(s))
SF = SF/np.max(SF)
sf = np.convolve(s,h)

omega = np.linspace(0,np.pi,M/2)


#plt.plot(s)
#plt.plot(n,s_ideal)
#plt.plot(n,sf[M/2:len(sf)-M/2+1])

#plt.plot(S[:len(SF)/2.])
#plt.plot(omega,SF[:len(SF)/2.])

#plt.plot(S[:1000]/np.max(S))
#plt.plot(SF[:1000]/np.max(SF))

#plt.plot(S[:100]/np.max(S))
plt.plot(omega,H[:M/2])

plt.xlabel('Frekvens (rad/s)',fontsize=13)
plt.ylabel('Amplitude',fontsize=13)






