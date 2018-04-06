#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 15:11:03 2018

@author: antonyschutz
"""
from synth_signal import generate_sound
import matplotlib.pyplot as pl
import numpy as np

from sdn import SDN

#%%
sig, a, f, time, Fs = generate_sound()
sdn = SDN(sig, Fs)


sdn.init(winovr = 50, winlen = 256, wintype = 'hann', length_opt = 'extended')   
sdn.sig2wframe()



#%%
sdn.frame_analysis()

#%%
sdn.frame2sig()

#%%
pl.plot(sdn.sig_modified)
pl.plot(sdn.rec_signal)



#%%
f = np.arange(2048)/4096
pl.clf()
pl.subplot(221) 
pl.plot(f,np.abs(sdn.frames_fft[:,0]))
pl.plot(sdn.frq_est[0], sdn.amp_est[0]/2,'o')


pl.subplot(222) 
pl.plot(sdn.frames[:,0])

pl.subplot(223) 
pl.plot(sdn.sig_mat[:,0])

#%%
pl.close()
ax = pl.subplot(411)
pl.plot(time,sig)

ax = pl.subplot(412)
sdn.plot_signal(ax=ax, color='r')

ax = pl.subplot(413)
pl.plot(np.arange(len(sdn.sig_modified))/Fs,sdn.rec_signal)

ax = pl.subplot(414)
pl.plot(np.arange(len(sdn.sig_modified))/Fs,sdn.sig_modified - sdn.rec_signal)
