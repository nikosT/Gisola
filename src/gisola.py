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

import datetime, sys, argparse, yaml, logging
from obspy.clients.fdsn.client import Client
from obspy import Catalog, Stream, read_events, read
from argparse import RawTextHelpFormatter
from obspy.core import UTCDateTime
import textwrap, os, glob

# local libs
import event, inventory, stream, isola, plot, web, config

__copyright__ = "🄯 "+str(datetime.datetime.utcnow().year)+", Institute of Geodynamics - National Observatory of Athens"
__credits__ = ["Nikolaos Triantafyllis (triantafyl@noa.gr), Ioannis Venetis (venetis@unipi.gr), Ioannis Fountoulakis (ifountoul@noa.gr), Erion-Vasilis Pikoulis (pikoulis@ceid.upatras.gr), Efthimios Sokos (esokos@upatras.gr), Christos Evangelidis (cevan@noa.gr)"]
__author__= 'Nikolaos Triantafyllis (triantafyl@noa.gr)'
__license__ = "GPLv3"
__version__ = "1.0"
__maintainer__ = "Nikolaos Triantafyllis"
__email__ = "triantafyl@noa.gr"
__status__ = "Production"

os.environ['COLUMNS'] = "90"

def print_help(func):
    def inner():
        print("""  ____ _           _       
 / ___(_)___  ___ | | __ _ 
| |  _| / __|/ _ \| |/ _` |
| |_| | \__ \ (_) | | (_| |
 \____|_|___/\___/|_|\__,_|
                           """) 
        print('R e a l – T i m e   M o m e n t   T e n s o r   C o m p u t a t i o n\n')
        print('Version: ' + __version__)
        print('License: ' + __license__)
        print('Author: ' + __author__)
        print ('Credits: '+ textwrap.TextWrapper().fill(text=''.join(__credits__)))
        print( __copyright__+'\n')
        func()
    return inner

class LineWrapRawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(' ', text).strip()
        return textwrap.wrap(text, width)


parser = argparse.ArgumentParser(prog='Gisola', description='Real-Time Moment Tensor Computation\nFind more info at: https://github.com/nikosT/Gisola',     formatter_class=LineWrapRawTextHelpFormatter)

# decorate (overload) print help message
parser.print_help=print_help(parser.print_help)

help_event='EVENT can be in any of the following formats: \n(i) DATETIME MAGNITUDE TYPE LATITUDE LONGITUDE DEPTH\ne.g.: ' + __file__ + ' -e ' + datetime.datetime.utcnow().isoformat() + ' 3.0 ML 37.24 20.49 4.1\n(ii) DATETIME\ne.g.: ' + __file__ + ' -e ' + datetime.datetime.utcnow().isoformat()+'\n(iii) EVENTID (event identifiers are data center specific)\ne.g.: ' + __file__ + ' -e noa2020owyrp\nIn cases (ii) and (iii) the rest of the information is retrieved by the FDSNWS-event \nIn more than one results, only the first event is returned \nPassing milliseconds is optional \nPassing more events is also posible using -e or --event prefix for each event'

parser.add_argument('-c', '--config', metavar=('FILEPATH'), help='override default configuration file (./config.yaml)', nargs='?', default='./config.yaml', type=str)
parser.add_argument('-e', '--event', metavar=('EVENT'), help=help_event, action='append', type=str, nargs='+', default=None)
parser.add_argument('--event-file', metavar=('FILEPATH'), help='parse and run a file with EVENT lines', type=argparse.FileType('r'), default=None)
parser.add_argument('--event-xml', metavar=('FILEPATH'), help='parse and run a file in QuakeML', type=argparse.FileType('r'), default=None)
parser.add_argument('-d', '--datetime-range', metavar=('TIME'), help='invoke MT computation for all events found in specific datetime range (FDSN bounded)', action='append', nargs=2, default=None)
parser.add_argument("--real-time", action="store_true", help='invoke --datetime-range for real-time use (FDSN bounded)', default=False)
parser.add_argument("--home", metavar=('YEAR'), help='Render Home page for particular year (based on config.yaml results dir)', action='store', nargs=1, default=None)

parser.add_argument("-s", "--station", metavar=('STA[.NEZ]'), help='override default stations selection. Optionally, components could be also specified.\nIt can be combined with --remove for the reverse result', action="extend", nargs="+", type=str, default=[])
parser.add_argument("--remove", help='it can only be used with --station. Invokes reverse result', action="store_true", default=False)
parser.add_argument("--disable-filters", help='disable all filter checks in process', action="store_true", default=False)
parser.add_argument("-cru", "--crustal", metavar=('FILEPATH'), help='override default crustal models', action="extend", nargs="+", type=str, default=None)
parser.add_argument("-dist", "--distance", metavar=('STRING'), help='override default distance rule', action="append", nargs=4, type=str, default=None)
parser.add_argument("-tl", "--time-window", metavar=('FLOAT'), choices=[245.76, 327.68, 409.6], help='override default inversion time-window', action="extend", nargs=1, type=float, default=None)
parser.add_argument("-ts", "--time-shift", metavar=('INT'), help='override default time search grid', action="append", nargs=3, type=int, default=None)
parser.add_argument("-f", "--frequency", metavar=('FLOAT'), help='override default inversion frequency band', action="append", nargs=4, type=float, default=None)
parser.add_argument("-r", "--revise", action="store_true", default=False, help='Trigger revised MT procedure. Computation assumes all available input exists\n(e.g. Greens\' Functions) and runs only the inversion and plot procedures.\nOnly --frequency, --time-shift, --station (--remove), --config' )
parser.add_argument('-l', '--log', metavar=('FILEPATH'), help='override default main log file (./log)', nargs='?', default='./log', type=str)
parser.add_argument('-v', '--version', action='version', version='%(prog)s v'+__version__)

args=parser.parse_args(args=None if sys.argv[1:] else ['--help'])

###### Procedure  ######

_cfg=config.read(args.config)
cat=Catalog()

logger = config.setup_logger('gisola_log', args.log)

if args.home:
    web.commandRenderHome(cfg=_cfg,year=args.home[0])

if args.event_xml:
    cat=read_events(args.event_xml)

# if --event-file is given override all other event input
if args.event_file:
    with args.event_file as _f:
        args.event=[_e.rstrip().split(' ') for _e in _f.readlines() if not _e[0]=='#']

# if --event-file or --event argument is passed
if args.event and not args.revise:
    for text in args.event:
        try:
            # if event is FDSN bound
            if len(text)==1:
                try:
                    cat.append(event.getFDSNWSCatalog(_cfg, logger, eventid=text[0])[0])
                except:
                    cat.append(event.getFDSNWSCatalog(_cfg, logger, starttime=UTCDateTime(text)-1, endtime=UTCDateTime(text)+1)[0])
            # create a manual object
            elif len(text)==6:
                cat.append(event.getCatalog(text)[0])
        except:
            pass

# by pass configuration options
if args.disable_filters:
    _cfg['Stream']['Modules']=[]

if args.distance:
    _cfg['Inventory']['Distance']=[[0,0, [_[0], _[1], _[2:]]] for _ in args.distance]

if args.crustal:
    _cfg['Green']['Crustal']=[]
    for text in args.crustal:
        _cfg['Green']['Crustal'].append({'Filepath': text, 'Geobox': None})

if args.time_window:
    _cfg['Inversion']['Window']=[[0,0,_] for _ in args.time_window]

if args.time_shift:
    _cfg['Inversion']['TimeShift']=[[0,0,[*_]] for _ in args.time_shift]

if args.frequency:
    _cfg['Inversion']['Frequency']=[[0,0,[*_]] for _ in args.frequency]

# if --datetime-range is given append event found by FDSN
if args.datetime_range:
    start=args.datetime_range[0][0]
    end=args.datetime_range[0][1]

    minMag=min([_[0] for _ in _cfg['Inversion']['Window']])
    cat+=event.getFDSNWSCatalog(_cfg, logger, minMag=minMag, starttime=UTCDateTime(start), endtime=UTCDateTime(end))

if args.real_time:
    cat+=event.watch(_cfg, logger)

# by-pass azimuthal coverage rule
if args.station and not args.remove:
    _cfg['Inventory']['Azimuth']=[1,_cfg['Green']['MaxStations']]

if not args.revise:
    # for each event returned from catalog object
    for evt in cat:
        try:
            config.init(evt, _cfg)
        except Exception as e:
            print('Fatal Error in init: '+str(e))
            continue

        try:
            #config.logger.info(('Loading initial configuration\n '+ \
            #'{}').format(config.dump(config.cfg)))
            org=event.getOrigin(config.cfg,evt,config.cfg['Watcher']['Historical'])
            config.logger.info('Starting MT calculation for Event {}\n{}\n{}\n'.format(evt, org,event.getMagnitude(evt,org)))
            config.logger.info('Setting working dir to: {}'.format(config.workdir))

            #if args.station and args.remove:
            #    for sta in args.station:
            #        config.blackList.append(sta.split('.')[0])

            #elif args.station and not args.remove:
                # make gisola keep whitelist only (however distance rule 
                # must return the stations here)
            #    config.magPriority=1
            #    config.cfg['Inventory']['WhiteList']['Priority']=0
            #    config.whiteList=[]
            #    for sta in args.station:
            #        config.whiteList.append([sta.split('.')[0],1])
    
            inventory.getInventory()
            config.logger.info(config.inv)

            if not args.station:
                inventory.selectDist()
                config.logger.info(config.inv)
          
            stream.getWaveforms()

            stream.clean()

            inventory.selectAzm()

            if args.station and args.remove:
                for sta in args.station:
                    res=sta.split('.')

                    if len(res)==3:
                        for comp in list(res[2]):
                            for tr in config.st.select(station=res[0], channel=res[1], component=comp):
                                config.st.remove(tr)
                    else:
                        for tr in config.st.select(station=res[0], channel='*' if len(res)==1 else res[1]):
                            config.st.remove(tr)

            elif args.station and not args.remove:
                _st=Stream()
                for sta in args.station:
                    res=sta.split('.')
                    if len(res)==3:
                        for comp in list(res[2]):
                            for tr in config.st.select(station=res[0], channel=res[1], component=comp):
                                _st+=Stream(tr)
                    else:
                        for tr in config.st.select(station=res[0], channel='*' if len(res)==1 else res[1]):
                            _st+=Stream(tr)
                config.st=_st

            stream.prioritize()

            isola.getGreens()
            isola.getInversions()
            isola.gatherResults()
            plot.allplots()
 
        except:
            config.logger.exception('Fatal Error occurred')
            config.logger.info('Moving to next event, if any')

# if revise flag is enabled
else:
    # for each event returned from catalog object
    for eventid in args.event:
        try:
            eventid=eventid[0]
            _year=''.join([s for s in eventid if s.isdigit()][:4])
            _evt_dir=os.path.join(_cfg['WorkDir'],_year,eventid)
            dirs=[x for x in glob.glob(os.path.join(_evt_dir, '*')) if os.path.isfile(os.path.join(_evt_dir,os.path.basename(x),'output','event.xml'))]
            # if event id with no year in name, try to find it
            if not dirs:
                for root, dir, _ in os.walk(_cfg['WorkDir']):
                    if eventid in dir:
                        _evt_dir=os.path.join(root, eventid)
                        dirs=[x for x in glob.glob(os.path.join(_evt_dir, '*')) if os.path.isfile(os.path.join(_evt_dir,os.path.basename(x),'output','event.xml'))]
                        break
            latest=max(dirs, key=os.path.getctime)
            workdir=os.path.join(_evt_dir,os.path.basename(latest))
            evt=read_events(os.path.join(workdir,'output','event.xml'))[0]
            config.logger = config.setup_logger(str(UTCDateTime.now()), os.path.join(workdir,'output','log.revised'))

            # gather all inv1.dat and find the best (based on correlation)
            corr=[]
            inversions=sorted(os.listdir(os.path.join(workdir,'inversions')))
            for inversion in inversions:
                with open(os.path.join(workdir,'inversions', inversion,'inv1.dat')) as _:
                    content=_.readlines()
                for i,line in enumerate(content):
                    if line.startswith(' Selected source position for subevent'):
                        corr.append(float(content[int(content[i+1].split()[1])+\
                        2].split()[2]))
                        break

            # open best corr inv1 file
            bestinvdir=inversions[corr.index(max(corr))]

            with open(os.path.join(workdir,'allstat','allstat'+bestinvdir.split('.')[0]+'.dat'),'r') as _f:
                lines=_f.readlines()
            lines=list(map(str.split,lines))

            for line in lines:
                line[5]=line[5] if not args.frequency else str(args.frequency[0][0])
                line[6]=line[6] if not args.frequency else str(args.frequency[0][1])
                line[7]=line[7] if not args.frequency else str(args.frequency[0][2])
                line[8]=line[8] if not args.frequency else str(args.frequency[0][3])

            for sta in args.station:
                res=sta.split('.')

                if len(res)==1:
                    res.extend([None,'NEZ'])

                if len(res)==3:
                    for line in lines:
                        if line[0]==res[0]:
                            if bool(int(line[2])):
                                line[2]=int('N' in list(res[2])) if not args.remove else int('N' not in list(res[2]))
                            if bool(int(line[3])):
                                line[3]=int('E' in list(res[2])) if not args.remove else int('E' not in list(res[2]))
                            if bool(int(line[4])):
                                line[4]=int('Z' in list(res[2])) if not args.remove else int('Z' not in list(res[2]))
                            line[1]=(int(line[2]) or int(line[3]) or int(line[4]))
                            break

            if not args.remove and args.station:
                sta_list=[sta.split('.')[0] for sta in args.station]

                for line in lines:
                    if line[0] not in sta_list:
                        line[1]=line[2]=line[3]=line[4]=0
             
            config.st=read(os.path.join(workdir,'streams_corrected.mseed'))
       
            for line in lines:
                if line[2]==0:
                    if config.st.select(station=line[0], component='N'):
                        config.st.remove(config.st.select(station=line[0], component='N')[0])
                if line[3]==0:
                    if config.st.select(station=line[0], component='E'):
                        config.st.remove(config.st.select(station=line[0], component='E')[0])
                if line[4]==0:
                    if config.st.select(station=line[0], component='Z'):
                        config.st.remove(config.st.select(station=line[0], component='Z')[0])

            lines=[list(map(str,line)) for line in lines]
            with open(os.path.join(workdir,'allstat','allstat'+bestinvdir.split('.')[0]+'.dat.revise'),'w') as _f:
                _f.writelines(' '.join(line)+'\n' for line in lines)

            isola.calculateRevisedInversions(_cfg,logger,workdir,bestinvdir)
            isola.gatherResults(cfg=_cfg, evt=evt, workdir=workdir, bestinvdir=bestinvdir, revise=True)

            plot.allplots(revise=True)
 
        except:
            logger.exception('Fatal Error occurred')
            logger.info('Moving to next event, if any')

