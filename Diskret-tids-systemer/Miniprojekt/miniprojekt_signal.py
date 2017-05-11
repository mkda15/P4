
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

N = 2**15
inter = 11*np.pi

def f(x):
    return np.sin(np.pi/3*x)

def g(x):
    return np.sin(2*np.pi/3 + np.pi/2*x)

def h(x):
    return np.sin(4*np.pi/3 + 3*np.pi/4*x)

x = np.linspace(0,inter,N)
y = f(x)+g(x)+h(x)

#plt.plot(x,f(x),label='f(x)')
#plt.plot(x,g(x),label='g(x)')
#plt.plot(x,h(x),label='h(x)')
plt.plot(x,y)
plt.legend(loc='upper right')






