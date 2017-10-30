# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 17:20:53 2017

@author: Srishty Saha
"""
from __future__ import division
import numpy as N
import numpy.fft as F
import numpy as np;
import scipy as sp;
from scipy.stats import signaltonoise
from scipy.signal import argrelextrema
from pylab import *
import operator
import cmath
import time
import numpy
from PIL import Image
import matplotlib.pyplot as plt
import numpy
from scipy import stats
from scipy.fftpack import fft,ifft
from scipy.io import wavfile
fs, data = wavfile.read('C:/Users/Srishty Saha/Desktop/da/fft-python-master/input/audio.wav') # load the data
print (data) 
a = data.T[0] # this is a two channel soundtrack, I get the first track
b=[(ele/2**8.)*2-1 for ele in a]

def czt(x, m=None, w=None, a=1.0):
    """
    Copyright (C) 2000 Paul Kienzle
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  US
    usage y=czt(x, m, w, a)
    Chirp z-transform.  Compute the frequency response starting at a and
    stepping by w for m steps.  a is a point in the complex plane, and
    w is the ratio between points in each step (i.e., radius increases
    exponentially, and angle increases linearly).
    To evaluate the frequency response for the range f1 to f2 in a signal
    with sampling frequency Fs, use the following:
    m = 32;                          ## number of points desired
    w = exp(-2i*pi*(f2-f1)/(m*Fs));  ## freq. step of f2-f1/m
    a = exp(2i*pi*f1/Fs);            ## starting at frequency f1
    y = czt(x, m, w, a);
    If you don't specify them, then the parameters default to a Fourier 
    transform:
      m=length(x), w=exp(2i*pi/m), a=1
    Because it is computed with three FFTs, this will be faster than
    computing the Fourier transform directly for large m (which is
    otherwise the best you can do with fft(x,n) for n prime).
    TODO: More testing---particularly when m+N-1 approaches a power of 2
    TODO: Consider treating w,a as f1,f2 expressed in radians if w is real
    """
    # Convenience declarations
    ifft = F.ifft
    fft = F.fft

    if m is None:
        m = len(x)
    if w is None:
        w = N.exp(2j*N.pi/m)

    n = len(x)

    k = N.arange(m)
    Nk = N.arange(-(n-1), m-1)

    nfft = next2pow(min(m,n) + len(Nk) -1)
    Wk2 = w**(-(Nk**2)/2)
    AWk2 = a**(-k) * w**((k**2)/2)
    w1=(fft(Wk2,nfft) * fft(x*AWk2, nfft))
    plt.plot(w1,'r') 
    plt.show()
    y = ifft(w1);
    f = w**((k**2)/2) * y[n:m+n]
    wavfile.write('C:/Users/Srishty Saha/Desktop/da/fft-python-master/input/audio7.wav', fs, f.astype(data.dtype))
    return y


def SNR(signal):
    max_index, max_value = max(enumerate(signal), key=operator.itemgetter(1))
    leftsignal = signal[0:max_index];
    rightsignal = signal[max_index:];

    leftMin = array(leftsignal);
    rightMin = array(rightsignal);

    findLMin = argrelextrema(leftMin, np.less)[0][-1];
    findRMin = argrelextrema(rightMin, np.less)[0][0]+len(leftsignal);

    x = np.linspace(0, 100,len(signal));



    Anoise = np.mean(list(signal[0:findLMin])+list(signal[findRMin:]))
    #Asignal = 1-(signal[findLMin]+signal[findRMin])/2
    Asignal = 1-Anoise;

    snr_value = 20*np.log10(Asignal/Anoise);
    plot(x[0:findLMin], signal[0:findLMin],'b')
    plot(x[findLMin:findRMin],signal[findLMin:findRMin],'r')
    plot(x[findRMin:],signal[findRMin:],'b');
    plot([x[max_index], x[max_index]],[1, 1- Asignal],'r--');
    plot(x, x*0+Anoise,'b--');
    show();

    print (snr_value)
    return snr_value    
def next2pow(x):
    return 2**int(N.ceil(N.log(float(x))/N.log(2.0)))

s=(czt(b))
print(SNR(s))