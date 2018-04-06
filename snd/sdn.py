#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 15:51:49 2018

@author: antonyschutz
"""

from __future__ import absolute_import, division
import matplotlib.pyplot as pl
import numpy as np
from windowing import windowing_signal, \
    reconstruct_signal, \
    analysis_window

from fourier_analysis import create_frames, my_specgr
    
__version__ = "0.0.0"    
    
# anwin = analysis_window(wintype, winlen, flag)   
    
class SDN(object):    
    
    def __init__(self, signal_input, sampling_frequency):
        
        # signal
        self.signal_input = signal_input    
        self.sampling_frequency = sampling_frequency            
        
        # input parameters
        self.winlen = None        
        self.wintype = None
        self.winovr = None
        self.length_opt = None
        
        # intermediary data
        self.analysis_windows = None
        self.rec_window = None
        self.sig_mat = None
        self.sig_modified = None  
        self.frames = None          
        self.frames_fft = None        
        self.amp_est = None
        self.frq_est = None
        self.pha_est = None          

        # intermediary parameter
        self.Nfft = None        
        # output
        self.rec_signal = None        

 
    def init(self, winlen, wintype='hann', winovr=50, length_opt='extended'):        
        """ Initialisation of SDN objects
        
        Parameters
        ----------
        winlen      :   int
                        length of analysis window
                        
        wintype     :   string
                        The analysis window
                                            rect
                                            hann
                                            hann_nonperiodic
                                            hamming
                                            blackman
                        
                        to add: triangle and other                        

        winovr      :   int
                        size of the overlap in percents      
                        
        length_opt  :   string
                        define if the signal is reduced or extended define also
                        the number of frame
        """                
        self.winlen = winlen 
        self.wintype = wintype
        self.winovr = winovr
        self.length_opt = length_opt

    def sig2wframe(self):
        """ Frame and window """                
        window = analysis_window(self.wintype, self.winlen, 'periodic')
        sig_mat, sig_modified = windowing_signal(self.signal_input, window, 
                                          self.winovr, self.length_opt)
        
        self.sig_mat = sig_mat
        self.sig_modified = sig_modified
        self.analysis_windows = window        
        
    def frame_analysis(self, Nfft = 4096):       
        """ Specgram and spectral content, frames estimation """             
        self.frames_fft = my_specgr(self.sig_mat, Nfft)        
        frames, f, a, p = create_frames(self.sig_mat.shape, self.frames_fft, \
                                        Nfft)
        self.frq_est = f
        self.amp_est = a
        self.pha_est = p        
        self.frames = frames
        self.Nfft = Nfft            
        
        
    def frame2sig(self):
        """ Reconstruct the signal from frame """                        
        window = analysis_window(self.wintype, self.winlen, 'symmetric')
        rec_signal = reconstruct_signal(self.frames, window, 
                                          self.winovr, len(self.sig_modified))
        self.rec_signal = rec_signal
        self.rec_window = window
           

    def plot_signal(self, ax=None, **kwargs):
        """ Plot the signal
        
        Parameters
        ----------
        ax        : matplotlib.Axes
                    The Axes instance in which the image is drawn
        kwargs : matplotlib.artist.Artist
                 Optional extra keyword/value arguments to be passed to
                 the ``ax.imshow()`` function
        """        
        if ax is None:
            ax = pl.gca()        
            
        time_input = np.arange(len(self.signal_input))/self.sampling_frequency            
            
        ax.plot(time_input, self.signal_input, **kwargs)
        ax.set_xlim((time_input.min(), time_input.max()))
        ax.set_xlabel('Time (second)')
        ax.set_ylabel('input signal')        
        