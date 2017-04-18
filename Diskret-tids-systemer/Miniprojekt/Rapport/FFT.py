#==============================================================================
# DFT
#==============================================================================
def DFT(x,c):
    X = np.zeros(c,dtype=complex)
    for k in range(len(x)): # Loop til alle samples i frekvens
        a = 0+0*1j
        for n in range(c):
            a += x[n]*np.exp(-2*np.pi*1j*k*n/float(c)) # Loop til samples i tid
            X[k] = a
    return X

#==============================================================================
# FFT
#==============================================================================
def FFT(x):
    N_new = len(x)
    if N % 2 > 0:
        raise ValueError('N skal vaere en potens af 2') # Brug N = potenser af 2
    elif N_new == 2:
        return DFT(x,N_new) # Returnerer DFT naar data ikke kan deles mere op
    else:
        X_even = FFT(x[::2]) # Deler rekursivt input op - lige dele
        X_odd = FFT(x[1::2]) # Deler rekursivt input op - ulige dele
        factor = np.exp(-2j * np.pi * np.arange(N_new) / N_new) # Twiddlefaktor
        return np.concatenate([X_even + factor[:N_new / 2] * X_odd,
                               X_even + factor[N_new / 2:] * X_odd])











