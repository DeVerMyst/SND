# -*- coding: utf-8 -*-
"""
Windowing and unwindowing of signal
"""
import numpy as np 

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