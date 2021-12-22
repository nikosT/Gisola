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

import os, inspect, sys, datetime, argparse, textwrap, os.path
from obspy import UTCDateTime, Stream
from obspy.core.event.origin import Origin
from obspy.geodetics.base import degrees2kilometers,  gps2dist_azimuth
from geojson import Point, Polygon, Feature
from turfpy.measurement import boolean_point_in_polygon, centroid
from turfpy.transformation import transform_scale
import numpy as np
import multiprocessing

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
parentdir = os.path.dirname(parentdir)
sys.path.insert(0, parentdir)

# local libs
import src.modules.ppsd as ppsd
import src.inventory as inventory
import src.stream as stream
import src.config as config

# app meta data
__copyright__ = "🄯 "+str(datetime.datetime.utcnow().year)+", Institute of Geodynamics - National Observatory of Athens"
__credits__ = ["Nikolaos Triantafyllis (triantafyl@noa.gr), Ioannis Venetis (venetis@unipi.gr), Ioannis Fountoulakis (ifountoul@noa.gr), Erion-Vasilis Pikoulis (pikoulis@ceid.upatras.gr), Efthimios Sokos (esokos@upatras.gr), Christos Evangelidis (cevan@noa.gr)"]
__author__= 'Nikolaos Triantafyllis (triantafyl@noa.gr)'
__license__ = "GPLv3"
__version__ = "1.0"
__maintainer__ = "Nikolaos Triantafyllis"
__email__ = "triantafyl@noa.gr"
__status__ = "Production"


# functions
def apply(trace, dirname='PPSD'):
    """
    Check for PPSD and return True/False
    """
    try:
        os.makedirs(os.path.join(parentdir,'src','modules', dirname))
    except:
        pass
    with open(os.path.join(parentdir,'src','modules', dirname, trace.get_id()), 'w') as f:
        f.write(str(ppsd.check(trace)))

def rescale(geobox, dist):
    """
    Rescale function
    """
    #Move to Polygon
    f = Feature(geometry=Polygon(geobox))
    #Calculate the cetroid of the polygon
    c = centroid(f)
    #Calculate the maximum distance of the polygon from the centroid in earth.
    dists = [gps2dist_azimuth(c['geometry']['coordinates'][0], c['geometry']['coordinates'][1], coord[0], coord[1])[0]/1000 \
            for coord in f['geometry']['coordinates'][0]]
    scale_ = 1 + dist/np.max(dists) 
   
    #Do the rescale
    # return Polygon
    return transform_scale(f, scale_, origin='centroid') 

class LineWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(' ', text).strip()
        return textwrap.wrap(text, width)


# programs parameters
parser = argparse.ArgumentParser(prog='PPSD_APP', description='Automatic PPSD inspection tool\nFind more info at: https://github.com/nikosT/Gisola',     formatter_class=LineWrapRawTextHelpFormatter)
parser.add_argument('-c', '--config', metavar=('FILEPATH'), help='override default configuration file (./config.yaml)', nargs='?', default='./config.yaml', type=str)
parser.add_argument('-l', '--log', metavar=('FILEPATH'), help='override default main log file (./log)', nargs='?', default='./log', type=str)
parser.add_argument("-th", "--threshold", metavar=('FLOAT'), help='assign threshold for setting PPSD coverage to False', action="store", nargs=1, type=float, default=70)
parser.add_argument("-r", "--range", metavar=('FLOAT'), help='set PPSD coverage check time range in minutes', action="store", nargs=1, type=float, default=24*60)
parser.add_argument("--verbose", help="increase output verbosity", action="store_true")
parser.add_argument('-v', '--version', action='version', version='%(prog)s v'+__version__)

args=parser.parse_args(args=None if sys.argv[1:] else ['--help'])


# reads configuration file
cfg=config.read(args.config)

# if no stationfile in white list found, exit
if not cfg['Inventory']['WhiteList']['Filepath']:
    print('Station file cannot be null in WhiteList-Filepath')
    sys.exit(0)

# setup log
logger = config.setup_logger('ppsd_app_log', args.log)
inventory.config.cfg=cfg
inventory.config.logger=logger
stream.config.logger=logger
ppsd.config.logger=logger


logger.info('Starting PPSD calculation')

# setup inventory request
# modify at your own needs
channels=list(set([x for _ in cfg['Inventory']['Distance'] for x in _[2][2]]))
inventory.config.distRules=[[0,degrees2kilometers(180),channels]]
inventory.config.windowRules=[3600*24-3600-1-10]
inventory.config.org=Origin()
inventory.config.org.latitude=0
inventory.config.org.longitude=0
inventory.config.org.time=UTCDateTime.now()-3600*24+10
inventory.config.magPriority=0
inventory.config.cfg['Stream']['Modules']=[]
inventory.config.blackList=[]
inventory.config.whiteList=[]
inventory.config.workdir='./'

if cfg['Inventory']['WhiteList']['Filepath']:
    with open(os.path.join(os.path.dirname(args.config),cfg['Inventory']['WhiteList']['Filepath']), 'r') as f:
        lines=f.readlines()
    for line in lines:
        if not(line.startswith('#') or line.startswith('\n')):
            net, sta, *_=line.split()
            priority=int(_[0]) if _ else int(cfg['Inventory']['WhiteList']['Priority'])

            if priority==0:
                inventory.config.blackList.append([net,sta])

inventory.getInventory(save=False)


# remove stations outside the Geobox
logger.info('Upscale Geobox based on longest distance and remove stations outside of it')

geobox=[[eval('('+x.replace(')','').replace('(','')+')') for x in cfg['Watcher']['Geobox'].replace(', ', ',').split('),(')]]
dist=max([rule[2][1] for rule in cfg['Inventory']['Distance']])

newgeobox=rescale(geobox,dist)

for sta in inventory.config.inv:
    if not boolean_point_in_polygon(Feature(geometry=Point((sta[0][0].longitude, sta[0][0].latitude))),newgeobox):
        inventory.config.inv = inventory.config.inv.remove(station=sta[0].code)


# setup waveform request
# remove seedlink from service is very slow
stream.config.cfg['Stream']['Service']=[s for s in stream.config.cfg['Stream']['Service'] if not s[0]=='SeedLink']

st=Stream()

_inv=inventory.config.inv.copy()

# one request per station
for station in sorted(list(set([sta.code for net in inventory.config.inv for sta in net]))):
    stream.config.st=Stream()
    stream.config.inv=_inv.select(station=station)
    try:
        stream.getWaveforms(save=False)

        # write True/False in trace file based on PPSD check
        with multiprocessing.Pool() as p:
            res=p.map(apply, stream.config.st)

    except:
        logger.exception('Error occurred for ' +station)
        pass

logger.info('End of PPSD calculation')

