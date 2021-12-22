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

import numpy as np

import sys, os, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import config

# Suggested Threshold values:
# broadband=0.9
# strong motion=0.98

def check(trace, threshold=0.98):
    """
    The main function
    """
    try:
        # work on a copy
        # leave original signal untouched
        tr=trace.copy()
        tr.detrend('demean')
        tr.detrend('linear')
        result=not clipAutoDetect(tr, threshold, sigma = 0.0001)[0]
        config.logger.info('Clipping check for {} passed: {}'.format(trace,result))
        return result
    except:
        config.logger.exception('Clipping check error occurred for ' +trace)
        return False    
        
def clipAutoDetect(trace, threshold = 0.9, sigma = 0.0001):
   """
   Auto Detect Clipping Algorithm by Wang and Zhang (2020).
   [Inputs]:
   trace: Obspy trace object
   threshold: Threshold for peak amplitude.Usually 0.9-0.98
   sigma: A small amplitude of random serial to purturb the peak amplitude. Usually 0.0001-0.001

   [Output]:
   Array with containing the Clipped Samples. If Array empty the signal is not Clipped.
   """

   #Create uniformly distributed random numbers array
   r = (np.random.uniform(low = 0.0, high = 1.0, size = trace.data.shape) - 0.5) * sigma * max(trace.data)

   #Add the r signal to the Trace data (avoid the failure for exactly equal or locally rebounding waveform around the peak amplitudes)
   trace.data = trace.data + r
 
   #Get the first order gradient
   g = np.diff(trace.data)

   #Pre-Allocate the array for the Clipped Data
   p = np.zeros(g.shape)

   #1) Two gradients in opposite sign, marked as Clipped Samples
   for j in range(1, len(g)-2, 1):
      if np.sign(g[j] * g[j+1]) <= 0:
         p[j] = trace.data[j]
         p[j+1] = trace.data[j+1] 
  
   #2) Exclude Under the Threshold
   _min=min(trace.data)
   _max=max(trace.data)
   _min_th=threshold * _min
   _max_th=threshold * _max
   _min_abs=abs(_min)
   _max_abs=abs(_max)
   for j in range(0, len(g), 1):
      if trace.data[j] < _max_th and trace.data[j] > _min_th:
         p[j] = 0

   #3) For Back to Zero type - If the gradient is absurdly large (say, about half of the peak amplitude), mark it by assigning it as peak amplitude (t       he sign is from its right neighbour)
   _max_factor=0.01 * _max
   for j in range(1, len(g)-1, 1):
      if (abs(trace.data[j]-trace.data[j-1]) + abs(trace.data[j+1]-trace.data[j])) > _max_th and (abs(trace.data[j+1]) > _max_th or abs(trace.data[j]) < _max_factor):
      
         p[j] = np.sign(trace.data[j+1]) * _max_abs
         p[j+1] = np.sign(trace.data[j+1]) * _max_abs    

   #4) For Back to Zero type - If two adjacent samples are absurdly large (say, about half of the peak amplitude), mark it by assigning it as peak amp       litude (the sign is from its right neighbour).

   for j in range(1, len(trace.data)-2,1):
      if abs(trace.data[j]-trace.data[j-1]) + abs(trace.data[j+2]-trace.data[j+1]) > _max_th and (abs(trace.data[j+2]) > _max_th or (abs(trace.data[j]) < 0.01 * _max and abs(trace.data[j+1]) < 0.01 * _max)):
     
        p[j] = np.sign(trace.data[j+2]) * _max_abs
        p[j+1] = np.sign(trace.data[j+2]) * _max_abs              

   #5) Mark the single unmarked sample that is between two marked samples.
   for j in range(1, len(p)-1, 1):
      if abs(p[j-1])> 0 and abs(p[j+1]) > 0 and abs(p[j]) == 0:
         p[j] = trace.data[j]

   #6) Mark the two unmarked samples that are between two marked samples.
   for j in range(1, len(p)-2, 1):
      if abs(p[j-1])> 0 and abs(p[j+2]) > 0 and abs(p[j+1]) == 0 and abs(p[j]) == 0:
         p[j] = trace.data[j]
         p[j+1] = trace.data[j+1]

   #7) Mark the three unmarked samples that are between two marked samples.
   for j in range(1,  len(p)-3, 1):
      if abs(p[j-1])> 0 and abs(p[j+3]) > 0 and abs(p[j+1]) == 0 and abs(p[j]) == 0 and abs(p[j+2]) == 0:
         p[j] = trace.data[j]
         p[j+1] = trace.data[j+1]
         p[j+2] = trace.data[j+2]

   #8) Mark the four unmarked samples that are between two marked samples.
   for j in range(1,  len(p)-4, 1):
      if abs(p[j-1])> 0 and abs(p[j+4]) > 0 and abs(p[j+1]) == 0 and abs(p[j]) == 0 and abs(p[j+2]) == 0 and abs(p[j+3])==0:
         p[j] = trace.data[j]
         p[j+1] = trace.data[j+1]
         p[j+2] = trace.data[j+2]
         p[j+3] = trace.data[j+3]

   #9) Mark the five unmarked samples that are between two marked samples.
   for j in range(1,  len(p)-5, 1):
      if abs(p[j-1])> 0 and abs(p[j+5]) > 0 and abs(p[j+1]) == 0 and abs(p[j]) == 0 and abs(p[j+2]) == 0 and abs(p[j+3])==0 and abs(p[j+4])==0:        
         p[j] = trace.data[j]
         p[j+1] = trace.data[j+1]
         p[j+2] = trace.data[j+2]
         p[j+3] = trace.data[j+3]
         p[j+4] = trace.data[j+4]

   #10) Back-to-zero type - Exclude the local single marked sample.
   for j in range(1, len(p)-1, 1):
      if p[j-1] == 0 and p[j+1] ==0 and abs(p[j])>0:
         p[j]=0

   #11) Back-to-zero type - Exclude local two marked samples.
   for j in range(1, len(p)-2, 1):
      if p[j-1] == 0 and p[j+2] == 0 and abs(p[j])>0 and abs(p[j+1])>0:
         p[j]=0
         p[j+1]=0

   #12) Back-to-zero type - Exclude local three marked samples.
   for j in range(1, len(p)-3, 1):
      if p[j-1] == 0 and p[j+3] == 0 and abs(p[j])>0 and abs(p[j+1])>0 and abs(p[j+2])>0:
         p[j]=0
         p[j+1]=0
         p[j+2]=0

   #13) Back-to-zero type - Exclude local four marked samples.
   for j in range(1, len(p)-4, 1):
      if p[j-1] == 0 and p[j+4] == 0 and abs(p[j])>0 and abs(p[j+1])>0 and abs(p[j+2])>0 and abs(p[j+3])>0:
         p[j]=0
         p[j+1]=0
         p[j+2]=0
         p[j+3]=0

   #14) Back-to-zero type - Exclude local five marked samples.
   for j in range(1, len(p)-5, 1):
      if p[j-1] == 0 and p[j+5] == 0 and abs(p[j])>0 and abs(p[j+1])>0 and abs(p[j+2])>0 and abs(p[j+3])>0 and abs(p[j+4])>0:
         p[j]=0
         p[j+1]=0
         p[j+2]=0
         p[j+3]=0
         p[j+4]=0
   
   if p[p!=0].size == 0:
      return (False, p)
   else:
      return (True, p)

