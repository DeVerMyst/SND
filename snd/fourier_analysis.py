# -*- coding: utf-8 -*-
"""
Fourier Analysis
"""


import numpy as np


def my_specgr(sig_mat, Nfft):

    X = np.fft.fft(sig_mat, n=Nfft, axis=0)/sig_mat.shape[0]

    Nfft2 = Nfft//2

    X = X[:Nfft2,:]

    return X


def peakdet(v, delta):
   
    v = np.array(v)

    nv = len(v)

    mv = np.zeros(nv+2*delta)

    mv[:delta] = np.Inf

    mv[delta:-delta]=v

    mv[-delta:] = np.Inf   

    ind = []
 
    # find all local maxima

    ind = [n-delta for n in range(delta,nv+delta) \
           if mv[n]>mv[n-1] and mv[n]>mv[n+1]]
   

    # take the maximum and closest from original estimation

    indi = np.array(ind,dtype=int) 

    return indi

 

def estimate_frame_prameters(spec_frame, Nfft):

    delta = 100

    power = np.abs(spec_frame)

    phase = np.angle(spec_frame)

    peaks = peakdet(power, delta)

    freq = peaks/Nfft

    amp = 2*power[peaks]

    pha = phase[peaks]

   

    return freq, amp, pha

 

def red_data(f,a,p,N):

    ind = np.argsort(a)
    if N < len(f):
        ind = ind[-N:]
        return f[ind], a[ind], p[ind]        
    else: 
        print('The number of found peak is less than %d' %N)
        return f, a, p


 

def est_ord(a, thresh=.025):

    '''estimate the number of harmonic based on a power ratio, power of all the

    maximum found, cumsum threshold'''

    ind = np.argsort(a)

    cs = np.cumsum(a[ind]**2)

    cs /= cs.max()

    N = np.sum(cs>thresh)

    return N

 

def create_frame(f, a, p, L_frame):
    

    t = np.arange(L_frame) 

    frame = np.zeros(L_frame)

    for n in range(len(f)):

        frame += 2*a[n]*np.cos(2*np.pi*t*f[n] + p[n])

    return frame

 

def create_frames(size_out, X, Nfft, prm):
    
    L_frame, N_frame = size_out

    Frames = np.zeros(size_out)
    a_list = []
    p_list = []
    f_list = []    
    for n in range(N_frame):       

        f, a, p = estimate_frame_prameters(X[:,n], Nfft)

        if prm<1 :
            N_f = est_ord(a, prm)
        elif type(prm)==int and prm>1:
            N_f = prm
        else:
            raise ValueError('Parameter must be <1 or a positive integer')

        
        f, a, p = red_data(f, a, p, N_f)
        
        f_list.append(f)
        a_list.append(a)
        p_list.append(p)        
        
        frame = create_frame(f, a, p, L_frame)

        Frames[:,n] = frame

    return Frames, f_list, a_list, p_list
