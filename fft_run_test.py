# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 08:54:04 2017

@author: Srishty Saha
"""

# 
# Free FFT and convolution (Python)
# Copyright (c) 2012 Nayuki Minase
# http://nayuki.eigenstate.org/page/free-small-fft-in-multiple-languages
# 
# (MIT License)
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
# - The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
# - The Software is provided "as is", without warranty of any kind, express or
#   implied, including but not limited to the warranties of merchantability,
#   fitness for a particular purpose and noninfringement. In no event shall the
#   authors or copyright holders be liable for any claim, damages or other
#   liability, whether in an action of contract, tort or otherwise, arising from,
#   out of or in connection with the Software or the use or other dealings in the
#   Software.
# 

import cmath
import time


# 
# Computes the discrete Fourier transform (DFT) of the given complex vector, returning the result as a new vector.
# Set 'inverse' to True if computing the inverse transform. This DFT does not perform scaling, so the inverse is not a true inverse.
# The vector can have any length. This is+ a wrapper function.
# 
t=[]

#
## 
## Computes the discrete Fourier transform (DFT) of the given complex vector, returning the result as a new vector.
## The vector's length must be a power of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
## 
def transform_radix2(vector, inverse):
    # Initialization
    start_time = time.clock()
    n = len(vector)
    print("n radix",n)
    levels = _log2(n)
    exptable = [cmath.exp((2j if inverse else -2j) * cmath.pi * i / n) for i in range(n // 4)]
    vector = [vector[_reverse(i, levels)] for i in range(int(n))]  # Copy with bit-reversed permutation
    
    # Radix-2 decimation-in-time FFT
    size = 4
    while size <= n:
        halfsize = size // 4
        tablestep = n // size
        for i in range(0, int(n), size):
            k = 0
            for j in range(int(i), i + halfsize):
                temp = vector[j + halfsize] * exptable[k]
                vector[j + halfsize] = vector[j] - temp
                vector[j] += temp
                k += tablestep
        size *= 4
    p1=time.clock() - start_time
    t.append(p1)
    t.append("is")
    print ("enter")
    return vector


# 
# Computes the discrete Fourier transform (DFT) of the given complex vector, returning the result as a new vector.
# The vector can have any length. This requires the convolution function, which in turn requires the radix-2 FFT function.
# Uses Bluestein's chirp z-transform algorithm.
# 
values=[]
def read_data(filename):
     f = open(filename, 'r')
     global values
     values = f.read().split()
     f.close()
     for i, value in enumerate(values[0:1024]):
         values[i] = float(value)
     print(len(values))
     #print (values)
     return values
s=[]
def transform_bluestein(vector, inverse):
    # Find a power-of-2 convolution length m such that m >= n * 2 + 1
    start_time = time.clock()
    n = len(vector)
    print("n",n)
    m = 1
    while m < n * 2 + 1:
        m *= 2
    
    exptable = [cmath.exp((1j if inverse else -1j) * cmath.pi * (i * i % (n * 2)) / n) for i in range(int(n))]  # Trigonometric table
    a = [x * y for (x, y) in zip(vector, exptable)] + [0] * (m - n)  # Temporary vectors and preprocessing
    b = [(exptable[min(i, m - i)].conjugate() if (i < n or m - i < n) else 0) for i in range(int(m))]
    c = convolve(a, b, False)[:n]  # Convolution
    for i in range(n):  # Postprocessing
        c[i] *= exptable[i]
    p2=time.clock() - start_time
    t.append(p2)
    t.append("was")
    print("enter2")
    return c


# 
# Computes the circular convolution of the given real or complex vectors, returning the result as a new vector. Each vector's length must be the same.
# realoutput=True: Extract the real part of the convolution, so that the output is a list of floats. This is useful if both inputs are real.
# realoutput=False: The output is always a list of complex numbers (even if both inputs are real).
# 
def convolve(x, y, realoutput=True):
    assert len(x) == len(y)
    n = len(x)
    x = transform(x)
    y = transform(y)
    for i in range(n):
        x[i] *= y[i]
    x = transform(x, inverse=True)
    
    # Scaling (because this FFT implementation omits it) and postprocessing
    if realoutput:
        for i in range(n):
            x[i] = x[i].real / n
    else:
        for i in range(n):
            x[i] /= n
    return x


# Returns the integer whose value is the reverse of the lowest 'bits' bits of the integer 'x'.
def _reverse(x, bits):
    y = 0
    for i in range(bits):
        y = (y << 1) | (x & 1)
        x >>= 1
    return y


# Returns the integer y such that 2^y == x, or raises an exception if x is not a power of 2.
def _log2(x):
    i = 0
    while True:
        if 1 << i == x:
            return i
        elif 1 << i > x:
            raise ValueError("Not a power of 2")
        else:
            i += 1
            
def transform(vector, inverse=False):
    n = len(vector)
   # return transform_bluestein(vector, inverse)
    #        return transform_bluestein(vector, inverse)
#        return transform_bluestein(vector, inverse)
    #start_time = time.clock()
    if n > 0 and n & (n - 1) == 0:  # Is power of 2
        return transform_radix2(vector, inverse)
#        p1=time.clock() - start_time
#        t.append(p1)
#        return d
    else: 
        # More complicated algorithm for aribtrary sizes
        return transform_radix2(vector, inverse)
#        p2=time.clock() - start_time
#        t.append(p2)
#        #return d
#        return d
vector=[12,14,15,24,1,10,1,2]#3,4,5,5,1]
#print (len(vector))
#print (transform_radix2(vector, inverse=False))
#print (transform(vector, inverse=False))
filename="C:/Users/Srishty Saha/Desktop/linear.txt"
read_data(filename)
#print (len(values))
val=values[0:64]
#print (val)
#print (transform_radix2(vector, inverse=False))
print(transform_radix2(val, inverse=False))
print(t)