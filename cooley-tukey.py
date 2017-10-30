# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 11:59:56 2017

@author: Srishty Saha
"""
from __future__ import division
import numpy as np;
import scipy as sp;
from scipy.stats import signaltonoise
from scipy.signal import argrelextrema
from pylab import *
import operator
import matplotlib.pyplot as plt
import numpy
from scipy import stats
from scipy.fftpack import fft,ifft
from scipy.io import wavfile # get the api
fs, data = wavfile.read('C:/Users/Srishty Saha/Desktop/da/fft-python-master/input/audio.wav') # load the data
print (data)
'''try:
     singleChannel = numpy.sum(data, axis=1)
     print singleChannel
     norm = singleChannel / (max(numpy.amax(singleChannel), -1 * numpy.amin(singleChannel)))
     print stats.signaltonoise(norm)
except:
     # was mono after all
     pass'''

a = data.T[0] # this is a two channel soundtrack, I get the first track
b=[(ele/2**8.)*2-1 for ele in a] # this is 8-bit track, b is now normalized on [-1,1)
c = fft(b) # calculate fourier transform (complex numbers list)
d = len(c)/2
plt.plot(c,'r') 
plt.ylabel('Gain in dB')
plt.xlabel('Frequency in Hz')
plt.show()
f = ifft(c)

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
#def snr(data):
#    singleChannel = data
#    try:
#      singleChannel = numpy.sum(data, axis=1)
#    except:
#      # was mono after all
#      pass
# 
#    norm = singleChannel / (max(numpy.amax(singleChannel), -1 * numpy.amin(singleChannel)))
#    print (norm)
    return (stats.signaltonoise(norm))

     # was mono after all
    
print(SNR(f))
wavfile.write('C:/Users/Srishty Saha/Desktop/da/fft-python-master/input/audio6.wav', fs, f.astype(data.dtype))


#d = len(c)/2  # you only need half of the fft list (real signal symmetry)
#plt.plot(abs(c[:(d-1)]),'r') 
#plt.show()

"""import scipy
import wave
import struct
import numpy
import pylab

from scipy.io import wavfile

rate, data = wavfile.read('C:/Users/neeti/Desktop/UMBC/Sem 3/DAA/Project/fft-python-master/input/audio.wav')

filtereddata = numpy.fft.fft(data)

print (data)

filteredwrite = numpy.fft.ifft(filtereddata)

print (filteredwrite)

wavfile.write('C:/Users/neeti/Desktop/UMBC/Sem 3/DAA/Project/fft-python-master/input/audio.wav', rate, filteredwrite)"""