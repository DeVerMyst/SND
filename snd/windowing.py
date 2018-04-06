# -*- coding: utf-8 -*-
"""
Windowing and unwindowing of signal
"""
import numpy as np 

def analysis_window(wintype, winlen, flag):
    ''' 
    Common windows used to analyse signal.
    
    wintype     :   string
                    The analysis window
                    rect
                    hann
                    hann_nonperiodic
                    hamming
                    blackman
                    
                    to add: triangle and other
                    
    winlen      :   int
                    length of analysis window
                    
    flag        :   string
                    flag to define if the analysis window is periodic or
                    symmetric. Depends on application, periodic for Fourier
                    analysis and symmetric for filter analysis.
                    
    ''' 
    winname = []                    
    winname.append('rect')    
    winname.append('hann')    
    winname.append('blackman')   
    
    if flag is 'periodic': # fourier analysis
        L = winlen + 1 
    elif flag is 'symmetric': # filter analysis
        L = winlen     
    else:
        print('flag can be: periodic (fourier), symmetric (filter)')
        raise ValueError('bad flag choice for window')    
        
    if wintype == 'rect':

        anwin = np.ones(L)

    elif wintype == 'hann':
        
        n = np.arange(L)
        anwin = .5*(1 - np.cos( 2*np.pi*n/(L-1) ))

    elif wintype == 'hamming':

        n = np.arange(L)
        anwin = .54 - .46*np.cos( 2*np.pi*n/(L-1))
                    
    elif wintype == 'blackman':

        n = np.arange(L)
        anwin = .42 - .50*np.cos( 2*np.pi*n/(L-1)) \
                    + .08*np.cos( 4*np.pi*n/(L-1))      

    else:
        print('available windows are: ', winname)
        raise ValueError('bad window choice')

    return anwin[:winlen]



def windowing_signal(signal_input, window, winovr, length_opt):

    ''' windowing of the signal in small part of len winlen, with overlapping 
    windows of size winovr. The analysis function is of type wintype

    signal_input:   array 1D (to do stereo)
                    The mono audio or 1D signal to cut into window
                    
    windows     :   array (see analysis_window)
                    the analysis window
                    
    winovr      :   int 
                    overlap in percents of the analysis window   
                    
    winlen      :   int  (see analysis_window)
                    length of analysis window                                                    
                    
    length_opt  :   string 
                    define if the signal is reduced (truncated) or extended 
                    (fulfilled with zeros)                   
                    
    '''         

    # CHECK OVERLAP value
    if winovr < 50:
        print('ovr < 50% lead to oscillations in succesive frame window')

    winlen = len(window)
    # CHECK Signal length compared to windows lenght        
    if winlen>len(signal_input):                
        raise ValueError('windows length must be > length of the signal ')    

    # output signal can different size
    signal = signal_input.copy()
    
    N = len(signal)    
        
    coef = 1 / ( 1 - winovr/100 )      
    step = int( ( 1 - winovr/100 ) * winlen )
            
    if length_opt == 'extended':
        
        if int(N/winlen)!=N/winlen:
            
            len_o = N + step
            o = signal.copy() 
            signal = np.zeros(len_o)
            signal[:N] = o
            N+=step
            
    elif (length_opt!='extended') and (length_opt!='reduced'):
        
        print('length_opt option: reduced (default), extended')
        raise ValueError('bad length_opt option')        
        
    Nwindow = int( np.floor( (N / winlen) * coef ) - np.floor( coef - 1 ) )
    
    N_o = (winlen + (Nwindow-1)*step )
    
    signal = signal[: N_o]

    sig_mat = np.zeros((winlen,Nwindow))

    start = 0

    for n in range( Nwindow ):

        sig_mat[:,n] = signal[start:start+winlen] * window

        start+=step

    return sig_mat, signal

   


def reconstruct_signal(sig_mat, window, winovr, sig_len):

    winlen = len(window)
    if len(sig_mat.shape) != 2:
        raise ValueError('sig_mat must be of size winlen*N_frames')         
        
#check sig_mat is on side with winlen        
    if sig_mat.shape[0] != winlen:
        sig_mat = sig_mat.T
        
    Nwindow = sig_mat.shape[1]

    out_sig = np.zeros(sig_len)   
    nrm_win = np.zeros(sig_len)       

    step = int( (1 - winovr/100)*winlen )   

    start = 0   

    for n in range( Nwindow ):

        out_sig[start:start+winlen] += sig_mat[:,n] * window
        nrm_win[start:start+winlen] += window
        start+=step   

    out = out_sig/nrm_win
    out[0]=0
    out[-1]=0   

    return out
