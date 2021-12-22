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

# Code initially retrieved from http://geo.mff.cuni.cz/~vackar/mouse/ and modified accordingly

# MouseTrap is an open-source module for Python for detecting a special type of disturbance 
# (so called mouse, ping, or fling step) in seismic records.
# Copyright: Jiří Vackář.
# Version: The most recent version is 0.2.1, released on 2015-04-28.
# License: GNU Lesser General Public License, Version 3 (http://www.gnu.org/copyleft/lesser.html)

import numpy as np
from numpy import pi, sin, cos, arctan
from obspy.core import stream, Trace, UTCDateTime
from obspy.core.util.attribdict import AttribDict
import matplotlib.pyplot as plt
import math
from math import sin, cos, radians
import sys, os, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import config

class mouse:
	def __init__(self, min_step = 2., min_step2 = 0.1, fit_time_before = 115, fit_time_after = 60):
		self.min_step        = min_step
		self.min_step2       = min_step2
		self.fit_time_before = fit_time_before
		self.fit_time_after  = fit_time_after

	def create(self, paz, length, dt, onset, typ=1, demean=True, integrate=True):
		self.nmax = nmax = int(pow(2,math.ceil(math.log(length, 2)))) #Changed
		self.df = df = 1. / (nmax*dt)
		self.dt = dt
		self.onset = onset
		self.typ = typ
		mys = self.mouse = np.zeros(nmax, dtype=complex)
		mys[0] = 0+0j
		for i in range(1, int(nmax/2+1)): # negative frequencies nmax/2+1 .. nmax-1 are complex conjugate of positive frequencies
			f = i*df
			if (typ == 2):
				mys[i] = (1./dt) / (2*np.pi*1j*f) # input velocity= step (i.e. input accel   = delta); mouse II
			elif (typ == 3):
				mys[i] = 1. / dt # mouse type III
			else:
				mys[i] = (1./dt) / (2.*np.pi*1j*f)**2. # input velocity= linear  (i.e. input accel = step); mouse I
				# INPUT ACCEL STEP = 1 m/s^2;
			mys[i] = mys[i] * np.exp(-2*np.pi*1j*f*onset) # time shift
		# Apply effect of the instrument:
		for i in range(nmax):
			fc = i * df * 2 * np.pi * 1j
			H = 1.+0j
			for zero in paz.zeros:
				H *= (fc - zero)
			for pole in paz.poles:
				H /= (fc - pole)
			H *= paz.gain * paz.sensitivity
			mys[i] *= H
		for i in range(1, int(nmax/2)):
			mys[nmax-i] = np.conj(mys[i]) # DFT of the real signal has complex conjugate spectrum
		self.mouse = np.fft.ifft(mys)
		if demean:
			self.demean()
		if integrate:
			self.integrate()

	def demean(self):
		last = int(self.onset / self.dt - 10) # take 10 samples less to avoid possible noc-causal blend of the mouse onset
		sum = 0
		for i in range(last):
			sum += self.mouse[i]
		mean = sum/last
		self.mouse -= mean

	def integrate(self):
		"""
		Integrate the synthetic mouse using cumulative sum (:func:`~numpy.cumsum`).
		"""
		self.mouse = self.mouse.cumsum()
		self.mouse *= self.dt

	def fit_mouse(self, st, t_min=0, t_max=0):
		dt = self.dt; mys = self.mouse

		if type(st) == stream.Stream:
			self.reclen = length = len(st[0])
		else:
			self.reclen = length = len(st)
		if (dt < self.min_step):
			step = int(self.min_step / dt) # grid-search step at least 0.1 s
		else:
			step = 1
		if (dt < self.min_step2):
			step2 = int(self.min_step2 / dt) # grid-search step at least 0.1 s
		else:
			step2 = 1
		if (t_min < 0 or t_min/dt > length):
			t_min = 0
		if (t_max <= 0 or t_max/dt > length):
			t_max = length*dt
		self.a = np.zeros(length); self.fit = 0; self.shift = 0
		if type(st) == stream.Stream:
			self.phi = np.zeros(length); self.theta = np.zeros(length)
			for i in range(3):
				comp = st[i].stats.channel[2] # channels must be ??N, ??E, and ??Z
				if comp == 'N':
					N = st[i].data
				elif comp == 'E':
					E = st[i].data
				elif comp == 'Z':
					Z = st[i].data
				else:
					raise Exception('Unknown component %s' % comp)
		elif type(st) == stream.Trace:
			D = st.data
		else:
			raise TypeError('Parameter st must be of type obspy.core.stream.Stream or obspy.core.stream.Trace.')
		for shift in range(int(length-t_max/dt), int(length-t_min/dt-1), step):
			if type(st) == stream.Stream:
				self.invert(N, E, Z, shift)
			else:
				self.invert_1D(D, shift)
		if step2 < step: # more detailed time grid search around time found in the previous step
			T_min = self.shift-step+step2
			T_max = self.shift+step
			if self.shift == int(length-t_max/dt): # the first inverted time
				T_min = self.shift+step2 # to avoid serching outside the time window
			elif self.shift == int(length-t_min/dt-1): # the last inverted time
				T_max = self.shift
			for shift in range(T_min, T_max, step2):
				if type(st) == stream.Stream:
					self.invert(N, E, Z, shift)
				else:
					self.invert_1D(D, shift)
		self.onset_found = (length - self.shift) * dt
		self.amplitude = self.a[self.shift]
		if type(st) == stream.Stream:
			self.phi = self.phi[self.shift]
			self.theta = self.theta[self.shift]
		
	def invert(self, N, E, Z, shift):
		"""
		Analytically calculate mouse amplitude, phi, and theta from partial derivatives. If find better fit than the previously found, update mouse.fit and mouse.shift values.
		
		It is called from :meth:`~mouse.fit_3D`, no need for using it alone.
		
		:type N:  north-south component
		:param N: :class:`~numpy.ndarray` (trace.data)
		:type E:  east-west component
		:param E: :class:`~numpy.ndarray` (trace.data)
		:type Z:  vertical component
		:param Z: :class:`~numpy.ndarray` (trace.data)
		:type shift: integer
		:param shift: time shift between synthetic mouse and real record; shift=0 means the mouse onset is at the end of the record, shift>0 shifts mouse to earlier time
		"""
		dt = self.dt; mys = self.mouse
		t_mouse = (self.reclen - shift) * dt
		fit_t0 = int((t_mouse - self.fit_time_before) / dt)
		fit_t1 = int((t_mouse + self.fit_time_after) / dt)
		on = self.reclen - int(self.onset/dt)
		ms_N = 0; ms_E = 0; ms_Z = 0; ms_NE = 0; ms_NEZ = 0; mm = 0; fit_nom = 0; fit_denom = 0
		# determination of phi: the azimuth from which the acceleration step comes from
		for i in range(fit_t0, fit_t1):
			ms_N += mys[i+shift-on].real * N[i] # \frac{\sum_i m_i s_i^E}{\sum_i m_i s_i^N} = \frac{\sin \phi}{\cos \phi} = \tg \phi
			ms_E += mys[i+shift-on].real * E[i]
		if (ms_E < ms_N):
			ph = arctan(ms_E / ms_N)
		else:
			ph = pi/2 - arctan(ms_N / ms_E)
		# determination of theta: inclination from which the acceleration step comes from
		for i in range(fit_t0, fit_t1): #for i in range(length):
			ms_Z  += mys[i+shift-on].real * Z[i]
			ms_NE += mys[i+shift-on].real * (N[i]*cos(ph) + E[i]*sin(ph))
		if (ms_Z < ms_NE):
			th = arctan(ms_Z / ms_NE)
		else:
			th = pi/2 - arctan(ms_NE / ms_Z)
			if th>pi/2:
				th -= pi
		# determination of amplitude of acceleration step
		for i in range(fit_t0, fit_t1): #for i in range(length):
			ms_NEZ += mys[i+shift-on].real * (N[i]*cos(ph)*cos(th) + E[i]*sin(ph)*cos(th) + Z[i]*sin(th))
			mm = mm + mys[i+shift-on].real ** 2
		a = ms_NEZ / mm
		# calculate the misfit
		for i in range(fit_t0, fit_t1):
			fit_nom += (N[i]-mys[i+shift-on].real*a*cos(ph)*cos(th))**2 + (E[i]-mys[i+shift-on].real*a*sin(ph)*cos(th))**2 + (Z[i]-mys[i+shift-on].real*a*sin(th))**2
			fit_denom += N[i]**2 + E[i]**2 + Z[i]**2
		fit = 1 - fit_nom / fit_denom
		# save results
		if a < 0:
			a = -a
			th = -th
			ph += np.pi
		self.a[shift] = a
		self.phi[shift] = ph % (2*np.pi)
		self.theta[shift] = th
		# is it the best result yet?
		if (fit > self.fit):
			self.fit = fit
			self.shift = shift
			
	def invert_1D(self, D, shift):
		"""
		Analytically calculate mouse amplitude from partial derivative. If find better fit than the previously found, update mouse.fit and mouse.shift values.
		
		It is called from :meth:`~mouse.fit_3D`, no need for using it alone.
		
		:type D:  a component
		:param D: :class:`~numpy.ndarray` (trace.data)
		:type shift: integer
		:param shift: time shift between synthetic mouse and real record; shift=0 means the mouse onset is at the end of the record, shift>0 shifts mouse to earlier time
		"""
		dt = self.dt; mys = self.mouse
		t_mouse = (self.reclen - shift) * dt
		fit_t0 = int((t_mouse - self.fit_time_before) / dt)
		fit_t1 = int((t_mouse + self.fit_time_after) / dt)
		on = self.reclen - int(self.onset/dt)
		ms_NEZ = 0; mm = 0; fit_nom = 0; fit_denom = 0
		# determination of amplitude of acceleration step
		for i in range(fit_t0, fit_t1): #for i in range(length):
			try:
				ms_NEZ += mys[i+shift-on].real * D[i]
			except:
				print(i, shift, on, i+shift-on, len(mys), len(D))
				exit()
			mm = mm + mys[i+shift-on].real ** 2
		a = ms_NEZ / mm
		# calculate the misfit
		for i in range(fit_t0, fit_t1):
			fit_nom += (D[i]-mys[i+shift-on].real*a)**2
			fit_denom += D[i]**2
		fit = 1 - fit_nom / fit_denom
		# save results
		self.a[shift] = a
		# is it the best result yet?
		if (fit > self.fit):
			self.fit = fit
			self.shift = shift

	def exist(self, Ts=120, likehood=False):
		likehood = self.fit #- abs(self.onset_found - Ts) / 50
		print(likehood)
		if likehood > 0.35:
			return 1
		else:
			return 0
	
	def plot(self, st, outfile='', distance=5e3, ylabel='', xmin=100, xmax=300, legend=False, yaxes=True):
		npts = st[0].stats.npts
		samprate = st[0].stats.sampling_rate
		t = np.arange(0, npts / samprate, 1 / samprate)
		fit_t0 = self.onset_found - self.fit_time_before
		fit_t1 = self.onset_found + self.fit_time_after
		samp_0 = int((self.onset - self.fit_time_before) / self.dt)
		samp_1 = int((self.onset + self.fit_time_after ) / self.dt)
		t2 = np.arange(fit_t0, fit_t1 - 0.5/samprate, 1 / samprate)
		shifts = {st[0].stats.channel[2]:0, st[1].stats.channel[2]:distance, st[2].stats.channel[2]:-distance}
		colors = {'N':'r', 'E':'g', 'Z':'b'}
		plt.rcParams.update({'font.size': 22})
		plt.figure(figsize=(12, 6))
		for tr in st:
			ch = tr.stats.channel[2]
			print(len(t), len(tr.data), len(t2), len(self.mouse[samp_0:samp_1]))
			p = plt.plot(t, tr.data+shifts[ch], colors[ch], label = {'Z':'raw record','N':'','E':''}[ch])
			plt.text(xmin+(xmax-xmin)*0.05, shifts[ch]+distance*0.2, ch, color=colors[ch])
			if ch == 'E':
				m = plt.plot(t2, self.mouse[samp_0:samp_1]*self.amplitude*sin(self.phi)*cos(self.theta)+shifts[ch], colors[ch]+'--')
			elif ch == 'N':
				m = plt.plot(t2, self.mouse[samp_0:samp_1]*self.amplitude*cos(self.phi)*cos(self.theta)+shifts[ch], colors[ch]+'--')
			elif ch == 'Z':
				m = plt.plot(t2, self.mouse[samp_0:samp_1]*self.amplitude*sin(self.theta)+shifts[ch], colors[ch]+'--', label='simulated')
		plt.xlim(xmin, xmax)
		if legend:
			plt.legend(loc='upper right')
		plt.title(st[0].stats.network+':'+st[0].stats.station+', '+str(st[0].stats.starttime.date)+' '+str(st[0].stats.starttime.time)[:8])
		if ylabel:
			plt.ylabel(ylabel)
		plt.xlabel('time [s]')
		if yaxes:
			# plot y2 axis
			plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,3))
			xmin,xmax,ymin,ymax = plt.axis()
			ax2=plt.twinx()
			const = self.mouse[int((self.nmax+self.onset/self.dt) / 2)]
			ax2.set_ylim(ymin/const,ymax/const)
			ax2.set_ylabel('acceleration [m/s$^2$]')
			ax2.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
		else:
			plt.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
		if outfile:
			plt.savefig(outfile, bbox_inches='tight')
		else:
			plt.show()

def prepR(st, t, bitrate=10):
	"""
	This Function prepares the records for the Core of the MouseTrap Program.
	Replaces the the Noise tests functions

	"""
	#Remove the before the event mean of the Record
	#Calculates mean value from part of the records between starttime and the time ``t`` and removes this mean from entire record
	for tr in st:
		idx = int(t * tr.stats.sampling_rate) - 1
		mean = np.mean(tr.data[0:idx])
		tr.data = tr.data - mean

	st.filter('bandpass', freqmin=0.005, freqmax=8) # giannis addon

	#Move to Raw Displacement and decimate
	st.integrate()
	if (st[0].stats.sampling_rate >= 2*bitrate):
		factor = int(st[0].stats.sampling_rate / bitrate)
		while (factor > 16):
			st.decimate(5)
			factor = int(factor / 5)
		st.decimate(factor)


def mouseChaser(st, t_start_origin, t_min, t_max):
   """
   This function tries to find if the input Stream has Mouse or not, using the MouseTrap modules
   [Inputs]:
   st: Obspy Stream or trace
   t_start_origin: Record seconds before the event origin time.
   [Output]:
   True of False
   """
  
   response=st[0].stats.response

   #Find the Poles and Zeros of the Station 
   #The 1st stage has the Disp --> Vel info
   poles = response.get_paz().poles
   zeros = response.get_paz().zeros
   norm_fact = response.get_paz().normalization_factor
   stage_gain = response.get_paz().stage_gain
   sensitivity = response.instrument_sensitivity.value

   #Create the Paz dictionary
   paz = AttribDict()
   paz['zeros'] = zeros
   paz['poles'] = poles
   paz['sensitivity'] = sensitivity
   paz['gain'] = stage_gain


   #Prepare Records (integrate and downsample to 10 Hz)
   prepR(st, t_start_origin, bitrate=10)

   #Create the synthetic mouse
   length = max([tr.count() for tr in st])
   dt = st[0].stats.delta

   Mouse = mouse()
   Mouse.create(paz, length*2, dt, length*dt, 1)

   trs=[]

   #Run per trace (1D)
   for tr in st:
      Mouse.fit_mouse(tr, t_min=t_min, t_max=t_max)
      result=Mouse.exist(t_start_origin)

      config.logger.info('Mouse check for {} passed: {}'.format(tr,not bool(result)))
      if result:
         trs.append(tr.get_id())

   return trs


# inversion must not be less than 120 sec
def check(stream, t_start_origin=10, t_min=0, t_max=120, inventory=None):
    """
    [Inputs]:
    st: obspy stream object (station traces)
    t_start_origin: Record seconds before the event origin time.
    t_min, t_max: times after the startime to fit for the mouse
   
    [Optional]: inventory --> in case the responce is not attached in the traces
    [Outputs]: mouseChaser function
    """
    try:
        # deep copy the trace
        st=stream.copy()

        # if inventory is specified, attach it
        # otherwise, it is assumed that inventory is already attached in stream
        if inventory:
            st.attach_response(inventory)

        return mouseChaser(st, t_start_origin=t_start_origin, t_min=t_min, 
        t_max=t_max)

    except:
        return [tr.get_id() for tr in st]

