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
        
    def frame_analysis(self, Nfft=4096, opt_N_hrmnc=.05):       
        """ Spectrogram and spectral content, frames estimation 
        
        Parameters
        ----------
        Nfft        :   int
                        Number of point used in the fft computation
                        
        opt_N_hrmnc :   float or int
                        if <1 : All the local maximum of the power spectrum 
                                are estimated and sorted. The peaks which 
                                contribute less than opt_N_hrmnc percents are
                                not taken into account
                        if (int)>0 : The opt_N_hrmnc main peak are estimated
                                

        winovr      :   int
                        size of the overlap in percents      
                        
        length_opt  :   string
                        define if the signal is reduced or extended define also
                        the number of frame
        """       
        
        self.frames_fft = my_specgr(self.sig_mat, Nfft)        

        
        frames, f, a, p = create_frames(self.sig_mat.shape, self.frames_fft, \
                                        Nfft, opt_N_hrmnc)
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
        
        
    def draw_spectrogram(self, prm=True, log10=True, ax=None, **kwargs):
        """ Draw the Spectrogram and estimated parameters
        
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
            
        Nfft, N_frame = self.frames_fft.shape    
        signal_duration = len(self.sig_modified)/self.sampling_frequency            
        frame_time = np.arange(N_frame) * signal_duration/N_frame
            
        power_spectra = np.abs(self.frames_fft)**2
        if log10:
            power_spectra = np.log10( power_spectra  )
        XT = [frame_time.min(), frame_time.max(),0 , .5]
        ax.imshow(power_spectra, extent=XT, \
                  origin='lower', aspect='auto', **kwargs)
        
#        ax.set_xlim((frame_time.min(), frame_time.max()))
        ax.set_xlabel('Time (second)')
        ax.set_ylabel('Frequencies')  
        if prm:
            for n,f in enumerate(self.frq_est):
                a = self.amp_est[n]
                p = self.pha_est[n]           
                t = signal_duration/N_frame*n*np.ones(len(f))
                ax.scatter(t,f,c=a,marker=".",s=1)
         
        
        