
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as axes
import scipy.fftpack

M = 50
N = 2**15
fs = 2**5
td = 1/fs
t = td*N

x = np.linspace(0,N*td,N)
xf = np.linspace(0,1/float(2*td),N/2)

cut1 = 2*np.pi*1
cut2 = 2*np.pi*2


def h(x):
    if x != M/2:
        return (1/(np.pi*x))*(np.sin(cut2*(x))-np.sin(cut1*(x)))
    elif x == M/2:
        return 1-(cut2-cut1)/np.pi

def high(x):
    if x != M/2:
        return -(np.sin(cut2*(x-M/2)))/(np.pi*(x-M/2))
    elif x == M/2:
        return 1-cut2/np.pi

def low(x):
    if x != M/2:
        return (1/(np.pi*(x-M/2)))*(np.sin(cut1*(x-M/2)))
    elif x == M/2:
        return cut1/np.pi

y = np.zeros(N)
for i in range(N):
    y[i] = high(x[i])


Y = np.abs(np.fft.fft(y)[:N/2])*2/np.sqrt(N)


#plt.plot(x,y)
plt.plot(xf,Y)


