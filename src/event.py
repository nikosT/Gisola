#!/usr/bin/python3
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

from obspy import UTCDateTime
from obspy.clients.fdsn.client import Client
from obspy.core.event.origin import Origin
from obspy.core.event.magnitude import Magnitude
from obspy.core.event.catalog import Catalog
from obspy.core.event.event import Event
from turfpy.measurement import boolean_point_in_polygon
from geojson import Point, Polygon, Feature
import os, math

import config

def getCatalog(text):

    def getEvent(info):
        """
        Converting text from format:
        2020-06-04T18:12:03,3.0,ML,37.24,20.49,4.1
        to ObsPy Event object 
        """
        evt=Event()
        evt.event_type='earthquake'
        evt.event_type_certainty='suspected'
        org=Origin()
        org.time=UTCDateTime(info[0])
        org.latitude=float(info[3])
        org.longitude=float(info[4])
        org.depth=float(info[5])*1000
        mag=Magnitude()
        mag.mag=float(info[1])
        mag.magnitude_type=info[2]
        mag.origin_id=org.resource_id
        evt.origins.append(org)
        evt.magnitudes.append(mag)
        evt.preferred_origin_id=org.resource_id
        evt.preferred_magnitude_id=mag.resource_id
        return evt

    cat=Catalog()
    cat.append(getEvent(text))
    return cat

def getFDSNWSCatalog(cfg, log, minMag=None, starttime=None, endtime=None, eventid=None):
    """
    # Returns a Catalog of Events based on starttime-endtime or eventid from a FDSNWS-event
    # Cases of Events requests (gisola API):
    # -e eventid (FDSNWS-event required)
    # -e datetime (FDSNWS-event required)
    # -e datetime lat lon mag depth
    # -p starttime endtime (FDSNWS-event required)
    # -r (-p (now-offset) - now)
    # -f eventfile where f: datetime lat lon mag depth (multiple lines)
    """
    log.info('Connecting to FDSNWS-event host:  ' + cfg['Event']['Host'])
    fdsnws = Client(cfg['Event']['Host'])

    if starttime and endtime:
        log.info('Requesting event in timespan: {} - {}'.format(UTCDateTime(starttime), UTCDateTime(endtime)))
        log.info('Minimum Magnitude and type is set to: '+str(minMag)+ \
        ' and '+ (cfg['Watcher']['Magnitudetype'] if cfg['Watcher']['Magnitudetype'] else 'None') +' respectively')

    elif eventid:
        log.info('Requesting event with ID: ' + eventid)

    try:
        return fdsnws.get_events(eventid=eventid or None,
                            starttime=UTCDateTime(starttime) if starttime else None,
                            endtime=UTCDateTime(endtime) if endtime else None,                        
                            minmagnitude=minMag or None if not eventid else None,
                            magnitudetype=cfg['Watcher']['Magnitudetype'] or None if not eventid else None,
                            includeallorigins=True, includeallmagnitudes=True, includearrivals=True,
                            orderby='time-asc' if not eventid else None)

    except Exception as e:
        if e.__class__.__name__=='FDSNNoDataException':
            log.info('No event is found...')
        else:
            log.info(e)

def watch(cfg, log):

    minMag=min([_[0] for _ in cfg['Inversion']['Window']])
    tl=min([_[2] for _ in cfg['Inversion']['Window']])
    tl2=max([_[2] for _ in cfg['Inversion']['Window']])
    endtime=UTCDateTime.now()-math.ceil(tl)-2-cfg['Watcher']['Playback']
    starttime=endtime-math.ceil(tl2)-cfg['Watcher']['Range']

    #starttime='2021-04-18T16:26:13.877220Z' set for debugging
    #endtime='2021-04-18T16:26:33.877220Z' set for debugging

    cat=getFDSNWSCatalog(cfg, log, minMag=minMag, starttime=starttime, endtime=endtime)

    vcat=Catalog()
    for evt in (cat if cat else []):
        try:
            org=getOrigin(cfg,evt,cfg['Watcher']['Historical'])
            if (not cfg['Watcher']['Geobox']) or boolean_point_in_polygon(\
            Feature(geometry=Point((org.longitude, org.latitude))),\
            Polygon([[eval('('+x.replace(')','').replace('(','')+')') for x in cfg['Watcher']['Geobox'].replace(', ', ',').split('),(')]])):

                # check for quality of origin
                if qualityEvent(cfg, evt, org):
                    workdir=os.path.join(cfg['WorkDir'],str(org.time)[:4],os.path.basename(str(evt.resource_id)))
                    if not os.path.exists(workdir): 
                        vcat.append(evt)
                    else:
                        log.info('{}\nalready calculated. It is removed from process'.format(evt))
                else:
                    log.info('{}\ndoes not meet quality thresholds. It is removed from process'.format(evt))

        except Exception as e:
            log.info(e+'\nError in event parsing. Moving to next event, if any')

    log.info(vcat)
    return vcat        

def getOrigin(cfg, evt, historical):

    if not historical:
        # retrieve "best" origin info (if any), else last found
        try:
            org=evt.preferred_origin()
            if not org: raise
        except:
            try:
                org=sorted(evt.origins, key=lambda o: o.creation_info.creation_time)[-1]
            except:
                org=evt.origins[-1]

    else:
        try:
            for i, org in enumerate(sorted(evt.origins, key=lambda o: o.creation_info.creation_time)):
                # first found (mimics real-time status)
                if qualityEvent(cfg,evt,org):
                    return org
        except:
            for i, org in enumerate(evt.origins):
                print(org)
                # first found (mimics real-time status)
                if qualityEvent(cfg,evt,org):
                    return org

    return org

def getMagnitude(evt, org):
    # retrieve associated magnitude
    return [m for m in evt.magnitudes if m.origin_id==org.resource_id][0]

def getFocalMechanism(evt):
    # retrieve "best" focal info (if any), else last found
    try:
        fm=evt.preferred_focal_mechanism()
        if not fm: raise
    except:
        fm=evt.focal_mechanisms[-1]

    return fm

def qualityEvent(cfg, evt, org):
    try:
        inittime=min([o.creation_info.creation_time for o in evt.origins])
        magnitude=getMagnitude(evt, org)

        # if data exit for realtime scenario
        windowRules=config.rules(magnitude.mag, cfg['Inversion']['Window'])
        maxtl=max([tl for tl in windowRules])

        if org.time+math.ceil(maxtl)+2>UTCDateTime.now()-cfg['Watcher']['Playback']:
            return False

        # if quality check pass or timeout occurred 
#        if int(round(org.creation_info.creation_time-inittime,0))>=cfg['Watcher']['Quality']['Timeout'] or \
# here could be modified according to needs
        if org.time+math.ceil(maxtl)+2+cfg['Watcher']['Quality']['Timeout']>=UTCDateTime.now()-cfg['Watcher']['Playback'] or \
           org.evaluation_mode=='manual' or \
           (round(org.time_errors.uncertainty,1)<=cfg['Watcher']['Quality']['Time'] and \
            round(org.depth_errors.uncertainty/1000,1)<=cfg['Watcher']['Quality']['Depth'] and \
            round(org.depth_errors.uncertainty/1000,1)!=0 and \
            round(org.latitude_errors.uncertainty,1)<=cfg['Watcher']['Quality']['Latitude'] and \
            round(org.longitude_errors.uncertainty,1)<=cfg['Watcher']['Quality']['Longitude'] and \
            round(magnitude.mag_errors.uncertainty,1)<=cfg['Watcher']['Quality']['Magnitude']):

            return True

    except Exception as e:
        print(e)
        return False

