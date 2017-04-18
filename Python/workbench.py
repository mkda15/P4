
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as axes
import scipy.fftpack

M = 40
N = 2**15
fs = 2**5
td = 1/fs
t = td*N

x = np.linspace(0,t,N)
xf = np.linspace(0,1/float(2*td),N/2)

cut1 = 2*np.pi*0.2
cut2 = 2*np.pi*0.3

def high(x):
    if x != M/2:
        return -np.sin(cut2*(x-M/2))/(np.pi*(x-M/2))
    elif x == M/2:
        return 1-cut2/np.pi

def low(x):
    if x != M/2:
        return (1/(np.pi*(x-M/2)))*(np.sin(cut1*(x-M/2)))
    elif x == M/2:
        return cut1/np.pi

def window(x):
    return 0.54-0.46*np.cos(2*np.pi*x/M)

def db(x):
    return 20*np.log10(np.abs(x))

y = np.zeros(N)
for i in range(N):
    y[i] = low(x[i])+high(x[i])
y = y/np.max(np.abs(y))
yw = window(x)*y



#Y = np.fft.fft(y)[:N/2]
Y = np.abs(np.fft.fft(yw))
Y = 1-Y/np.max(np.abs(Y))

def f(x):
    return np.sin(np.pi/3*x)

def g(x):
    return np.sin(2*np.pi/3 + np.pi/2*x)

def h(x):
    return np.sin(4*np.pi/3 + 3*np.pi/4*x)

sig = f(x)+g(x)+h(x)
sig_g = f(x)+h(x)
#sig = sig/np.max(np.abs(sig))
SIG = np.abs(np.fft.fft(sig))
#SIG = SIG/np.max(np.abs(SIG))

filt_SIG = Y*SIG
#filt_sig = filt_sig/np.max(np.abs(filt_sig))
filt_sig = np.fft.ifft(filt_SIG)


#plt.plot(x,y)
#plt.plot(xf,Y)
plt.plot(x[0:1500],sig_g[0:1500])
#plt.plot(xf,SIG)
#plt.plot(xf,filt_sig)
plt.plot(x[0:1500],filt_sig[0:1500])

