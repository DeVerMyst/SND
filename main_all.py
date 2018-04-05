#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""

Created on Thu Mar 29 09:47:31 2018

 

@author: antonyschutz

"""

import numpy as np

import matplotlib.pyplot as pl

 

 

#%%

## tools audio decomposition

def generate_sound(Duration = 5, Fs = 44100, Number_harm = 10):

    ''' Generate random sound (note) of duration Duration second at sampling

    Frequency Fs in Hz with Number_harm harmonic

    Fondamental frequency will be between 80 and Fs/(2*Number_harm)'''

    f0 = 80 + np.random.rand()*int((Fs-80)/(2*Number_harm))

    

    amp = .2 + ( np.random.rand(Number_harm) - .2 )

    a = [f for f in amp]

    f = [g*f0 for g in range(1,1+Number_harm)]

    p = (np.random.rand(Number_harm)-.5)*np.pi

    N = Duration * Fs

    time = np.arange( N ) / Fs   

    signal = np.zeros( N )

        

    for n in range( Number_harm ):    

        signal += a[n] * np.cos( 2 * np.pi * (f[n]*time) + p[n] )

    

    return signal, a, f, time, Fs

 

#%%

def calc_full_Powerfft(x, Fs, Nfft = 4096):

    ''' compute the spectral density of x with Nfft points'''

    Nf = Nfft//2

    y = np.abs( np.fft.fft(x,n=Nfft) / len(x) )**2

    y = y[:Nf]

    freq = np.arange( int( Nf )  ) * Fs / Nfft

    return y, freq

 

#%%

def fen_an(wintype, winlen):

    if wintype == 'rect':

        fen = np.ones(winlen)

    elif wintype == 'hann':

        n = np.arange(winlen)

        fen = .5*(1 - np.cos( 2*np.pi*n/(winlen-1) ))

    elif wintype == 'hamming':

        n = np.arange(winlen)

        fen = .54 - .46*np.cos( 2*np.pi*n/(winlen-1))

    elif wintype == 'blackman':

        n = np.arange(winlen)

        fen = .42 - .5*np.cos( 2*np.pi*n/(winlen-1)) + .08*np.cos( 4*np.pi*n/(winlen-1))      

    else:

        print("PB")

    return fen

 

def windowing(signal, winlen, winovr=50, wintype='hann'):

    ''' windowing of the signal in small part of len winlen,

    with overlapping windows of size winovr.

    The analysis function is of type wintype'''         

    fen = fen_an(wintype, winlen)   

    N = len(signal)

    coef = 1 / ( 1 - winovr/100 )   

    Nwindow = int( np.floor( (N / winlen) * coef ) - np.floor( coef - 1 ) )

    jump = int( (1 - winovr/100)*winlen )   

    sig_mat = np.zeros((winlen,Nwindow))

    nrm_fen = np.zeros(int(N/winlen)*winlen)

    start = 0

    for n in range( Nwindow ):

        sig_mat[:,n] = signal[start:start+winlen] * fen

        nrm_fen[start:start+winlen] += fen

        start+=jump

   

    return sig_mat, nrm_fen

   

def unwindowing(sig_mat, nrm_mat, winlen, sig_len, winovr):

    Nwindow = sig_mat.shape[1]

    out_sig = np.zeros(int(sig_len/winlen)*winlen)   

    jump = int( (1 - winovr/100)*winlen )   

    start = 0   

    for n in range( Nwindow ):

        out_sig[start:start+winlen] += sig_mat[:,n]

        start+=jump   

    out = out_sig/nrm_mat

    out[0]=0

    out[-1]=0   

    return out

 

def red_sig(signal, winlen):

    N = len(signal)   

    return signal[:int(N/winlen)*winlen]

   

 

def my_specgr(sig_mat, Nfft):

#    X = np.fft.fft(sig_mat, n=Nfft, axis=0)/Nfft

    X = np.fft.fft(sig_mat, n=Nfft, axis=0)/sm.shape[0]

    Nfft2 = Nfft//2

    X = X[:Nfft2,:]

    return X

 

#%%

def peakdet(v, delta):

   

    v = np.array(v)

    nv = len(v)

    mv = np.zeros(nv+2*delta)

    mv[:delta] = np.Inf

    mv[delta:-delta]=v

    mv[-delta:] = np.Inf   

    ind = []

 

    # find all local maxima

    ind = [n-delta for n in range(delta,nv+delta) if mv[n]>mv[n-1] and mv[n]>mv[n+1]]

   

    # take the maximum and closest from original estimation

    indi = np.array(ind,dtype=int)

 

    return indi

 

def estimate_frame_prameters(spec_frame, Nfft, percents=0):

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

    ind = ind[-N:]

    return f[ind], a[ind], p[ind]

 

def est_ord(a, thresh=.025):

    '''estimate the number of harmonic based on a power ratio, power of all the

    maximum found, cumsum threshold'''

    ind = np.argsort(a)

    cs = np.cumsum(a[ind]**2)

    cs /= cs.max()

   

    N = np.sum(cs>thresh)

    return N

 

def create_frame(f, a, p, win_len, wintype, Fs):

    t = np.arange(win_len)

    frame = np.zeros(win_len)

    for n in range(len(f)):

        frame += 2*a[n]*np.cos(2*np.pi*t*f[n] + p[n]) * fen_an(wintype, winlen)

    return frame

 

def create_frames(size_out, X, Nfft, winlen, wintype):

    N_frame = X.shape[1]

    Frames = np.zeros(size_out)

    for n in range(N_frame):

       

        f, a, p = estimate_frame_prameters(X[:,n], Nfft)

        N_f = est_ord(a)

        f, a, p = red_data(f, a, p, N_f)

        frame = create_frame(f, a, p, winlen, wintype, Fs)

        Frames[:,n] = frame

    return Frames

 

#%%

winovr = 75

winlen = 256

Nfft = 4096

wintype = 'hann'

 

signal, a, f, time, Fs = generate_sound(Number_harm=10, Fs = 2**15)

 

y,freq = calc_full_Powerfft(signal, Fs)

sm,nm = windowing(signal, winlen, winovr=winovr, wintype=wintype)

 

 

X = my_specgr(sm, Nfft)

 

F = create_frames(sm.shape, X, Nfft, winlen, wintype)

 

s_out = unwindowing(F, nm, winlen, len(signal), winovr)

 

pl.clf()

 

pl.subplot(121)

pl.plot(time,signal)

pl.xlim( (0, 4/ f[0]) )

pl.xlabel('time (sec)')

pl.title('four period of the signal ')

 

pl.subplot(122)

pl.plot(freq,10*np.log10(y))

pl.xlim( (0, freq.max()) )

pl.xlabel('frequency (Hz)')

pl.title('power density (of full signal)')

 

pl.plot(f,10*np.log10( (np.array(a)/2)**2),'or')

 


 

pl.close('all')
pl.subplot(121)
pl.imshow(np.abs(X)**2,aspect='auto')

pl.subplot(122)
pl.plot(np.abs(X)**2)

 
 

pl.subplot(311)

pl.plot(signal)

pl.plot(s_out)

 

pl.subplot(312)

pl.plot(np.abs(np.fft.fft(signal)))

pl.plot(np.abs(np.fft.fft(s_out)))

 

pl.subplot(313)

pl.plot(signal - s_out)

 

 

#%%

#pl.close('all')

#fn = np.arange(len(x))/Nfft

#pl.plot(fn, np.abs(x)**2)

#pl.plot(f,(a/2)**2,'or')

 

 

##%%

#def next_power_of_2(x): 

#    return 1 if x == 0 else 2**(x - 1).bit_length()

#

#import scipy as sp

#def blanchir(s,p):

## blanchiment du signal s par filtrage médian puis estimation AR

## p est la longueur du filtre, le transitoire est supprimé.

#

#    a = filtre_blanchiment[s[:np.minimum(len(s),1000)],p]

#    x =  sp.signal.lfilter(b, a, s)

#    x = s[p:]

#    return x

#   

#def filtre_blanchiment(x,p):

## x: signal

## p: ordre de l'AR

#    e       = np.ones((p,1))

#    N       = len(x)

#    Nfft    = 2^next_power_of_2(N);

#   

#    q       = int(len(x)/5) # longueur du filtre median

#   

#    P       = np.abs(np.fft.fft(x * fen_an('hann', N), n=Nfft))**2 / N

#   

#    if np.mod(q,2):

#        P1 = np.vstack( (P[(-(q-1)/2 +1):] , P[:(q/2)] ))

#        P2 = np.vstack( (P[(q/2):] , P[:(q/2 - 1)] ))       

#        H = sp.linalg.hankel(P1,P2)

#    else:

#        H = hankel([P(end-q/2+1:end);P(1:q/2)],[P(q/2:end);P(1:q/2-1)]);

#    end

#   

#    H       = sort(H,1);

#    P       = H(round(q/3),:).';

#    rx      = real( ifft(P) );

#    a       = toeplitz(rx(1:p)) \ e;

#    a       = a / a(1);
