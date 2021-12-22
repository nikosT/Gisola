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

import numpy as np, os.path
from obspy import read_inventory, Inventory, Stream
from obspy.clients.fdsn.client import Client
from obspy.geodetics.base import kilometers2degrees
from obspy.geodetics.base import gps2dist_azimuth # for fast run: sudo -H pip3 install geographiclib !!!!
import operator, logging, math

# import local configuration 
import config

@config.time
def getInventory(save=True):
    """
    Downloading inventory from stationXML file or FDSNWS
    """
    # find min and max distance of accepted rules
    # for defining inventory geobox selection
    mindist=min([rule[0] for rule in config.distRules])
    maxdist=max([rule[1] for rule in config.distRules])

    # get the maximum time window to retrieve data
    maxtl=max([tl for tl in config.windowRules])

    if 'ppsd' in config.cfg['Stream']['Modules']:
        starttime=config.org.time-10-3600 # half hour
    else:
        starttime=config.org.time-10

    endtime=config.org.time+math.ceil(maxtl)+1

    config.logger.info(('Retrieving inventory with working epochs from {}' + \
    ' to {} and stations between {} and {} km from event\'s ' + \
    'location').format(starttime, endtime, mindist, maxdist))

    inv=Inventory()

    # use queue priority as priority host service
    for service in config.cfg['Inventory']['Service']:
        # loops will finish when no service is left for requesting data

        if service[0]=='FDSNWS':
            try:
                # connect to host
                fdsnws = Client(service[1], eida_token=service[2])
                config.logger.info('Connecting to FDSNWS host ' +service[1])

                # download inventory from host based on config rules
                _inv=fdsnws.get_stations(starttime=starttime, 
                     endtime=endtime, latitude=config.org.latitude, \
                     longitude=config.org.longitude, \
                     minradius=kilometers2degrees(mindist), \
                     maxradius=kilometers2degrees(maxdist), \
                     includerestricted=True, level='response')

            except Exception as e:
                if e.__class__.__name__=='FDSNNoDataException':
                    config.logger.info('Could not retrieve inventory from FDSNWS host '+service[1])
                else:
                    config.logger.info(e)
                _inv=Inventory()
                pass

        elif service[0]=='StationXML':
            try:
                config.logger.info('Reading Inventory from ' +service[1])
                _inv=read_inventory(service[1], format='STATIONXML').select(starttime=starttime, \
                     endtime=endtime, latitude=config.org.latitude, longitude=config.org.longitude, \
                     minradius=kilometers2degrees(mindist), \
                     maxradius=kilometers2degrees(maxdist))

            except:
                config.logger.info('Could not load Inventory file '+service[1])
                _inv=Inventory()
                pass

        # remove duplicates (keep only the first occurance)
        for name in list(set([sta.code for net in inv for sta in net])):
            _inv=_inv.remove(station=name)
        # append rest
        inv+=_inv

    _inv=inv
    # get only selected streams
    chas=[]
    # for each accepted channel get all possible component combinations
    for cha in [case for rule in config.distRules for case in rule[2]]:
        chas+=[cha+orient for comp in config.cfg['Inventory']['Components'] \
        for orient in comp]
    chas=sorted(list(set(chas)))

    config.logger.info('Accepted channels: '+ str(chas))

    inv=Inventory()
    for cha in chas:
        inv+=_inv.select(channel=cha)

    # keep unique station names
    _inv=Inventory()
    for name in list(set([sta.code for net in inv for sta in net])):
        _inv+=inv.select(station=name)

    # if magnitude priority is greater than default keep only the whiteList
    if config.cfg['Inventory']['WhiteList']['Filepath']:
        if config.magPriority > int(config.cfg['Inventory']['WhiteList']['Priority']):
            _inv2=Inventory()
            for sta in config.whiteList:
                _inv2+=_inv.select(network=sta[0],station=sta[1])
            _inv=_inv2

        else:
            # remove blacklisted (priority=0) stations
            for sta in config.blackList:
                _inv=_inv.remove(network=sta[0],station=sta[1])

    config.inv=_inv

    # save inventory
    if save:
        config.logger.info(('Saving inventory in one StationXML file in ' + \
        '{}').format(os.path.join(config.workdir,'inventory.xml')))
        config.inv.write(os.path.join(config.workdir,'inventory.xml'), \
        format='StationXML')

@config.time
def selectDist():
    """
    Selecting subset of inventory based on Distance Rules
    """
    _inv=Inventory()
    for rule in config.distRules:
        config.logger.info(('Selecting stations with channel codes {} and that' + \
        ' they are between {} and {} km from event\'s ' + \
        'location').format(rule[2], rule[0], rule[1]))
        for cha in rule[2]:
            _inv+=config.inv.select(latitude=config.org.latitude, \
                longitude=config.org.longitude, \
                minradius=kilometers2degrees(rule[0]), \
                maxradius=kilometers2degrees(rule[1]), \
                channel=cha+'*')
 
    config.inv=_inv

    # save inventory
    config.logger.info('Saving inventory in XML file {}'.format(os.path.join(
    config.workdir,'inventory.xml')))
    config.inv.write(os.path.join(config.workdir,'inventory.xml'), format='STATIONXML')  

@config.time
@config.log
def selectAzm():
    """
    Selecting subset of streams based on Azimuthal Coverage and Rules
    """
    config.logger.info('Selecting subset of streams based on Azimuthal '+ \
    'Coverage and Rules')
    count=0

    # disable them in order to choose the preferred
    for tr in config.st:
        tr.enabled=False

    stations=[]
    for _station in list(set([tr.stats.station for tr in config.st])):
        _inv=config.inv.select(station=_station)[0][0]
 
        # get priority
        priority=1 # default value in case of null stations file
        if config.cfg['Inventory']['WhiteList']['Filepath']:
            priority=int(config.cfg['Inventory']['WhiteList']['Priority'])
            for elem in config.whiteList:
                if _station==elem[0]:
                    priority=elem[1]
                    break

        stations.append((_station, priority, _inv.latitude, _inv.longitude, \
        gps2dist_azimuth(config.org.latitude, config.org.longitude, \
        _inv.latitude, _inv.longitude)))
   
    for deg in range(0,360,45):
        #sector=[_ for _ in stations if _[4][1]>=deg and _[4][1]<deg+45]

        sector=[(_[0],_[1],'VEL' if config.st.select(station=_[0])[0].stats.response.instrument_sensitivity.input_units in config.vel_types else 'ACC', config.st.select(station=_[0]).count(),-float(_[4][0])) for _ in stations if _[4][1]>=deg and _[4][1]<deg+45]

        # sort them by priority value
        sector=sorted(sector, key = operator.itemgetter(1,2,3,4), reverse=True)

        # cut to the max stations per sector
        sector=sector[:config.cfg['Inventory']['Azimuth'][1]]
        config.logger.info('Sector [{},{}): {}'.format(deg, deg+45, sector))

        if sector:
            count+=1

        # mark traces to be used
        for quad in sector:
            for tr in config.st.select(station=quad[0]):
                tr.enabled=True

    # write station, priority, lat, lon, distance, azimuth and backazimuth
    with open(os.path.join(config.workdir,'locations'), 'w') as f:
        f.writelines('{}\n'.format(x) for x in sorted(stations, key=operator.itemgetter(4)))

    # check number of sectors
    if count < config.cfg['Inventory']['Azimuth'][0]:
        config.logger.info('Azimuthal coverage less than configuration set. ' +\
        'Terminating process')
        raise

    # keep only the passed traces
    config.st=Stream([tr for tr in config.st if tr.enabled]).sort()

    config.locations=stations
