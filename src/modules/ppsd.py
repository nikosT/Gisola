#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#    Copyright (C) 2021 Triantafyllis Nikolaos

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

import os, inspect, timeit, sys, numpy as np
from obspy import read
from obspy.signal import PPSD
from obspy.core.inventory.inventory import read_inventory
from obspy.signal.spectral_estimation import PPSD #import load_npz
from obspy.signal.spectral_estimation import get_nhnm
from obspy.signal.spectral_estimation import get_nlnm
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import config


# inversion must not be less than 120 sec
def check(trace):

    try:
        inventory=config.inv.select(network=trace.stats.network, station=trace.stats.station, \
                                location=trace.stats.location, channel=trace.stats.channel)

        result=PPSDcoverage(trace, inventory)
        config.logger.info('PPSD check for {} passed: {}'.format(trace,result))
        return result
    except:
        config.logger.exception('PPSD check error occurred for ' +trace)
        return False  

#low_freq=0.005, high_freq=0.2
def PPSDcoverage(trace, inventory, low_freq=0.01, high_freq=0.2, threshold=70):
    """
    Calculate the percentage of data 
    that are between the noise model (Peterson, 1993)
    low_freq, high_freq: Lower and upper threshold in frequencies to calculate (Hz)
    thershold: percentage inside the (Peterson, 1993)
    """
    # calculate the PPSD
    ppsd = PPSD(trace.stats, metadata=inventory)
    ppsd.add(trace)
    #ppsd.plot(show_mode=True, period_lim=(0.01, 50), show_coverage=True, show_histogram=True)

    # set options for numpy
    precision=5
    np.set_printoptions(precision=precision, suppress=True, 
                        threshold=sys.maxsize)

    # retrieve periods and amplitudes calculated by PPSD
    periods, ampls=ppsd.get_mode()

    # calculate freqs (1/periods)
    freqs=1.0/periods

    # 1. round freqs by precision
    # 2. clip/limit freqs to desired band
    # 3. create list with unique freqs only
    dfreqs=np.unique(np.clip(np.around(freqs, decimals=precision),low_freq,high_freq))

    # retrieve periods, low and high amplitudes based on noise model (Peterson, 1993)
    (mperiods, mlampls), mhampls=get_nlnm(), get_nhnm()[1]

    # calculate mfreqs (1/mperiods)
    mfreqs=1.0/mperiods

    # count how many points are "good"
    count=0
    for frq in dfreqs:
        # find nearest indexes that discrete freq applies to
        # to retrieve the most nearest amplitude
        idx=(np.abs(freqs - frq)).argmin()
        midx=(np.abs(mfreqs - frq)).argmin()

        # if our data (get_mode) are between the noise model
        # then this point is inside the "good" area, so count it
        if ampls[idx] >= mlampls[midx] and ampls[idx] <= mhampls[midx]: 
            count+=1

    # calculate percentage of "good" points
    percentage = count*100/len(dfreqs)

    if percentage < threshold:
        return False
    else:
        return True
    #return percentage

