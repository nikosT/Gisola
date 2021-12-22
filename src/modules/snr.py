#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#    Copyright (C) 2021 Fountoulakis Ioannis

#    This file is part of Gisola.

#    Gisola is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, 
#    or any later version.

#    Gisola is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with Gisola.  If not, see <https://www.gnu.org/licenses/>.

import math
import numpy as np
import matplotlib.pyplot as plt
from obspy.signal.trigger import classic_sta_lta
from obspy.signal.filter import envelope
import sys, os, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import config

def check(trace):
    """
    The main function
    """
    try:
        # work on a copy
        # leave original signal untouched
        tr=trace.copy()

        result=getSNR(tr)
        config.logger.info('SNR check for {} passed: {}'.format(trace,result))
        return result
    except:
        config.logger.exception('SNR check error occurred for ' +trace)
        return False  

# the below functions is currently not in used
def smooth1D(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise(ValueError, "smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise(ValueError, "Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise(ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=numpy.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    return y

def trigCalc(trace, STA = 1, LTA = 60, threshold=20):
   """
   Simple STA/LTA trigger Algorithm
   [Inputs]:
   trace: Obspy trace
   STA-LTA in seconds
   threshold: For Pick detection
   [Output]:
   start_spot: trigger position
   """
   #Detrend-demean the trace
   trace.detrend('demean')
   trace.detrend('linear')

   #Sampling rate
   df = trace.stats.sampling_rate

   #Filter the signal (Bandpass 2-15 Hz
   trace.filter('bandpass', freqmin=2.0, freqmax=15.0)

   #Do the STA/LTA (trace!!)
   cft = classic_sta_lta(trace.data, int(STA*df), int(LTA*df))

   #Check where in the cft the values exceed the STA_threshold
   start_spot = np.where(cft>=threshold)

   return [start_spot[0], cft]


def trigCalcEnv(trace, STA = 1, LTA = 60, threshold=20):
   """
   Simple STA/LTA trigger Algorithm, with envelope and smoothing (Earle and Shearer 1994)
   [Inputs]:
   trace: Obspy trace
   STA-LTA in seconds
   threshold: For Pick detection
   [Output]:
   start_spot: trigger position
   """
   
   #Detrend-demean the trace
   trace.detrend('demean')
   trace.detrend('linear')

   #Sampling rate
   df = trace.stats.sampling_rate

   #Filter the signal (Bandpass 2-15 Hz
   trace.filter('bandpass', freqmin=2.0, freqmax=15.0)

   trace = envelope(trace.data)

   #Do the STA/LTA (trace!!)
   cft = classic_sta_lta(trace, int(STA*df), int(LTA*df))

   #Smooth the STA/LTA
   cft = smooth1D(cft, window_len=200, window='hanning')

   start_spot = np.where(cft>=threshold)

   return [start_spot[0], cft]


def SNRampl(trace, start_spot, a_time=10):
   """
   Simple SNR based on amplitude.
   [Inputs]:
   trace: Obpsy trace
   start_spot: trigger position
   a_time (sec): time around the trigger for the SNR calculation
   [Output]:
   SNR value
   """
   df = trace.stats.sampling_rate

   #Translate the seconds in number of points
   slice_pos = int(a_time * df)
   noise_window = trace.data[(start_spot-slice_pos):start_spot]
   signal_window = trace.data[start_spot: (start_spot + slice_pos)]

   #Calculate the SNR(Signal to Ratio) by dividing the RMS amplitude values
   Srms = np.sqrt(1/len(signal_window) * np.sum(np.square(signal_window)))
   Nrms = np.sqrt(1/len(noise_window) * np.sum(np.square(noise_window)))

   SNR = Srms/Nrms

   return SNR


def getSNR(trace, snr_threshold = 5):
  """
  [Inputs]:
  trace: Obspy trace
  [Outputs]:
  True or False
  """
  #Check if the signal has triggers
  trigs = trigCalc(trace, STA = 1, LTA = 20, threshold=8)

  #No triggers
  if len(trigs[0])==0:
     return False
  else:
  #We have trigger
     SNR = SNRampl(trace, trigs[0][0], a_time=5)

     if SNR < snr_threshold:
        return False
     else:
        return True

