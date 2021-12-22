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

import math, os, os.path, operator
from statistics import median
import numpy as np, multiprocessing
from collections import OrderedDict
import shutil, subprocess, logging, glob
from obspy.geodetics.base import gps2dist_azimuth # for fast run: sudo -H pip3 install geographiclib
from obspy.core.event.source import FocalMechanism,  MomentTensor, NodalPlanes, \
                                    NodalPlane, PrincipalAxes, Axis, Tensor
from obspy.core.event.base import WaveformStreamID, DataUsed, CreationInfo
from obspy import read, UTCDateTime
from obspy.core.event.event import Event
from obspy.core.event.origin import Origin, OriginQuality
from obspy.core.event.magnitude import Magnitude
from obspy.geodetics import kilometers2degrees
from obspy.core import AttribDict
#import geo.sphere
# local
import config, event, modules.kagan

def kaganCalc(x):
    return modules.kagan.get_kagan_angle(*x[0:3], *x[3:6])

def cart2earth(lat,lon,x,y):
    """
    BSD 3-Clause License
    Copyright (c) 2017, Juno Inc.
    All rights reserved.
    https://github.com/gojuno/geo-py/blob/master/geo/sphere.py
    Calculate lat and lon from origin lat lon and cart position
    """
    def destination(point, distance, bearing):
        '''
            Given a start point, initial bearing, and distance, this will
            calculate the destina?tion point and final bearing travelling
            along a (shortest distance) great circle arc.

            (see http://www.movable-type.co.uk/scripts/latlong.htm)
        '''
        EARTH_MEAN_RADIUS = 6371008.8

        lon1, lat1 = (math.radians(coord) for coord in point)
        radians_bearing = math.radians(bearing)

        delta = distance / EARTH_MEAN_RADIUS

        lat2 = math.asin(
            math.sin(lat1)*math.cos(delta) +
            math.cos(lat1)*math.sin(delta)*math.cos(radians_bearing)
        )
        numerator = math.sin(radians_bearing) * math.sin(delta) * math.cos(lat1)
        denominator = math.cos(delta) - math.sin(lat1) * math.sin(lat2)

        lon2 = lon1 + math.atan2(numerator, denominator)

        lon2_deg = (math.degrees(lon2) + 540) % 360 - 180
        lat2_deg = math.degrees(lat2)

        return (lon2_deg, lat2_deg)

    dist=math.sqrt(x**2+y**2)*1000
    try:
        if x==0 and y==0:
            azm=0
        elif y>=0 and x>=0:
            azm=math.atan(x/y)
        elif y<0 and x>=0:
            azm=180-math.atan(x/abs(y))
        elif y<=0 and x<=0:
            azm=180+math.atan(abs(x)/abs(y))
        elif y>0 and x<0:
            azm=270+math.atan(abs(y)/abs(x))
    except ZeroDivisionError:
        azm=90 if x>0 else -90

    lon2,lat2=destination((lon,lat), dist, azm)
    return round(lat2, 4), round(lon2, 4)


def calculateVariance(observed, synthetic, tl):
    """
    Calculates variance reduction of stations' streams
    """
    with np.errstate(divide='raise'):
        try:
            dt = round((tl/8192.0),2)
            obs = np.array(observed)
            syn = np.array(synthetic)

            ds = obs-syn
            dsn = (np.linalg.norm(ds)**2)*dt
            d = (np.linalg.norm(obs)**2)*dt

            return round(1-(dsn/float(d)),2)
        except FloatingPointError:
            return None

def createSources():
    """
    Creating the Source files (source.dat)
    """
    # for each accepted grid rule, perform the following actions
    for i,rule in enumerate(config.gridRules):

        # create respective source dir
        os.makedirs(os.path.join(config.sourcedir,'grid'+str(i)))

        # calculate source points
        # for x,y (y=x) distance
        x=[]
        for drange in rule['Distance']:
            x=np.append(x,np.unique(np.append(-np.arange(*drange), \
            np.arange(*drange))))

        # for z depth
        z=[]
        for zrange in rule['Depth']:    
            z=np.append(z,np.unique(np.append(-np.arange(*zrange), \
            np.arange(*zrange))))
       
        # use event's depth as offset
        z+=round(config.org.depth/1000,1) # km
        # keep only >= 1 km 
        z=z[z>=1.0]

        # all possible source points
        mesh=np.array(np.meshgrid(x,x,z)).T.reshape(-1, 3)

        # write one file per chunk
        for k in range(math.ceil(len(mesh)/config.cfg['Green']['MaxSources'])):

            # write chunck size in text buffer and then in file
            # with the step of chuck size (=MaxSources)
            text=""
            for j,elem in enumerate(mesh[k*config.cfg['Green']['MaxSources']: \
            (k+1)*config.cfg['Green']['MaxSources']]):
                text+='{}\t{:.4f}\t{:.4f}\t{:.4f}\n'.format(\
                k*config.cfg['Green']['MaxSources']+j+1, *elem)
            with open(os.path.join(config.sourcedir,'grid'+str(i),\
            'source'+str(k)+'.dat'), 'w') as _:
                _.write(text)

        config.logger.info(('Grid index: {}\nCalculated source points: {}\n' + \
        'Dispatched in {} source file(s)').format(i, \
        len(mesh), k+1))

def createCrustals():
    """
    Copying the Crustal Model files (crustal.dat)
    """
    # create respective source dir
    os.makedirs(config.crustaldir)

    for i,rule in enumerate(config.crustalRules):
        try:
            config.logger.info(('Crustal model Index: {}\nCopying ' + \
            'crustal file: {} to {} directory').format(i, \
            rule['Filepath'], config.crustaldir))
            shutil.copy(rule['Filepath'], os.path.join(config.crustaldir, \
            'crustal'+str(i)+'.dat'))
        except:
            config.logger.info('Crustal model file: {} not found'.format( \
            rule['Filepath']))
            config.logger.info('Continue with next crustal model, if any')
            continue

def createGrdat():
    """
    Creating the Greens Function configuration (grdat.hed)
    """
    maxdist=max([rule[1] for rule in config.distRules])
    xl=2000000 if maxdist<=100 else int(np.ceil(20*1000*maxdist))
    # max freq used for calculation
    maxfreq=2.5*max([rule[3] for rule in config.freqRules])

    # create respective source dir
    os.makedirs(config.grdatdir)

    for i,tl in enumerate(config.windowRules):
        nfreq=np.ceil(tl*maxfreq)
        text=('&input\nnfreq={:n}\ntl={}\naw=1.0\nxl={}\nikmax=100000\n' + \
        'uconv=0.1E-03\nfref=1.\n/end\n').format(nfreq,tl,xl)
 
        with open(os.path.join(config.grdatdir,\
        'grdat'+str(i)+'.hed'), 'w') as _:
            _.write(text)

        config.logger.info('Greens\' Functions configuration index: {}'.format(i))

def createStations():
    """
    Creating the necessary Greens Function configuration (station.dat)
    """
    text=''

    lstations=[]
    for i,sta in enumerate(list(set([tr.stats.station for tr in config.st]))):
        distance, azimuth, _ = gps2dist_azimuth(config.org.latitude, \
        config.org.longitude, config.inv.select(station=sta)[0][0].latitude, \
        config.inv.select(station=sta)[0][0].longitude)

        lstations.append([distance,distance*math.cos(math.radians(azimuth))/1000.0, \
        distance*math.sin(math.radians(azimuth))/1000.0,0,sta])

    # sort by distance and write the necessary info
    for i,sta in enumerate(sorted(lstations, key=operator.itemgetter(0))):
        text+='{:d}\t{:.4f}\t{:.4f}\t{:.4f}\t{}\n'.format(i+1,*sta[1:])
    
    with open(os.path.join(config.workdir,'station.dat'), 'w') as _:
        _.write(text)

@config.time
def calculateGreens():
    """
    Calculating Greens Function (gr.hes)
    """
    def getIdx(text):
        return ''.join(list(filter(str.isdigit, text)))

    # create respective source dir
    os.makedirs(config.greendir)

    for crustal in sorted(os.listdir(config.crustaldir)):

        for grdat in sorted(os.listdir(config.grdatdir)):

            for grid in sorted(os.listdir(config.sourcedir)):
                for source in sorted(os.listdir(os.path.join(\
                config.sourcedir,grid))):
                    try:
                        grhes='gr.{}.{}.{}.{}.hes'.format(getIdx(crustal),
                                                      getIdx(grdat),
                                                      getIdx(grid),
                                                      getIdx(source))

                        command='{} {} {} {} {} {}\n'.format(\
                                        os.path.join(os.getcwd(),config.cfg['Green']['ExePath']),
                                        os.path.join('station.dat'),
                                        os.path.join('crustals', crustal),
                                        os.path.join('grdat', grdat),
                                        os.path.join('sources', grid, source),
                                        os.path.join('greens', grhes))

                        proc=subprocess.Popen(command, cwd=config.workdir, \
                        shell=True, universal_newlines=True,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                        out,err=proc.communicate()
                        config.logger.info(err)
                        config.logger.info(out)
                    except:
                        config.logger.info('It could not create {}'.format(grhes))
                        continue
 
def createInpinv():
    """
    Creating Inversions configuration files (inpinv)
    """
    os.makedirs(config.inpinvdir)

    for i,srule in enumerate(config.shiftRules):
        text='{} {} {}'.format(*srule)
        with open(os.path.join(config.inpinvdir, \
        'inpinv'+str(i)+'.dat'), 'w') as _:
            _.write(text)

def createAllstat():
    """
    Creating Inversions configuration files (inpinv)
    """
    # create respective source dir
    os.makedirs(config.allstatdir)

    # read all stations from station.dat file
    with open(os.path.join(config.workdir,'station.dat'), 'r') as _:
        stations=_.readlines()

    # create allstat dir
    for i,rule in enumerate(config.freqRules):
        text=''
        for stat in stations:
            sta=stat.split()
            _st=config.st.select(station=sta[4])

            if _st:
                text+='{} {} {} {} {} {}\n'.format(sta[4],
                      1 if _st.count() else 0,
                      1 if _st.select(component='N') else 0,
                      1 if _st.select(component='E') else 0,
                      1 if _st.select(component='Z') else 0,
                      ' '.join(map(str, rule)))

            else:
                text+='{} {} {} {} {} {}\n'.format(sta[4], 0, \
                      0, 0, 0, ' '.join(map(str, rule)))


            with open(os.path.join(config.allstatdir, \
            'allstat'+str(i)+'.dat'), 'w') as _:
                _.write(text[:-1])

def createRaw():
    """
    Creates raw file of data; contains 4 columns: time, N, E, Z data
    if there's less than 3 components data, it fills with ajacent data
    (dummy data) in order ISOLA to work
    """
    # create respective source dir
    os.makedirs(config.rawdir)

    # read all stations from station.dat file
    with open(os.path.join(config.workdir,'station.dat'), 'r') as _:
        stations=_.readlines()

    # break sts in unique triplets
    sts=list(set([(tr.stats.station, tr.stats.location, \
    tr.stats.channel[:-1]) for tr in config.st])) 

    # create rawfiles for each accepted tl
    for i,tl in enumerate(config.windowRules):

        # one dir for each tl found
        os.makedirs(os.path.join(config.rawdir,str(i)))

        _st=config.st.copy()
        # downsampling at known frequency at 8192 elements
        _st.resample(8192/tl)
        # be sure than no point exceeds number
        for tr in _st:
            tr.data=tr.data[:8192]

        for stat in stations:
            _st2=_st.select(station=stat.split()[4])

            text=""
            if _st2.count():
                for k,time in enumerate(_st2[0].times()):
                    text+='{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\n'.format(time, 
                          _st2.select(component='N')[0].data[k] if _st2.select(component='N') else 0,
                          _st2.select(component='E')[0].data[k] if _st2.select(component='E') else 0,
                          _st2.select(component='Z')[0].data[k] if _st2.select(component='Z') else 0
                          )

            else: # dummy data
                for _ in range(8192):
                    text+='{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\n'.format(0, 0, 0, 0)

            with open(os.path.join(config.rawdir,str(i),stat.split()[4].upper()+'raw.dat'), 'w') as _:
                _.write(text)

@config.time
def calculateInversions():

    # create respective source dir
    os.makedirs(os.path.join(config.inversiondir))

    for allstat in sorted(os.listdir(config.allstatdir)):

        for inpinv in sorted(os.listdir(config.inpinvdir)):

            for grhes in sorted(os.listdir(config.greendir)):

                _, icrustal, igrdat, igrid, isource, _ = grhes.split('.')

                invdir='{}.{}.{}.{}.{}.{}'.format(allstat.split('.')[0][7:], 
                       inpinv.split('.')[0][6:], icrustal, igrdat, igrid, isource)

                os.makedirs(os.path.join(config.inversiondir,invdir))

                command='{} {} {} {} {} {} {} {} {}\n'.format(os.path.join(os.getcwd(),config.cfg['Inversion']['ExePath']),
                            os.path.join('..', '..','allstat', allstat),
                            os.path.join('..', '..','inpinv', inpinv),
                            os.path.join('..', '..','grdat', 'grdat'+igrdat+'.hed'),
                            os.path.join('..', '..','greens', grhes),
                            os.path.join('..', '..','crustals', 'crustal'+icrustal+'.dat'),
                            os.path.join('..', '..','sources', 'grid'+igrid, 'source'+isource+'.dat'),
                            os.path.join('..', '..','station.dat'),
                            os.path.join('..', '..','raw', igrdat))
                #print(command)

                proc=subprocess.Popen(command, cwd=os.path.join(config.inversiondir, invdir), shell=True, universal_newlines=True,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                out,err=proc.communicate()
                config.logger.info(err)
                config.logger.info(out)


def calculateRevisedInversions(_cfg,logger,workdir,bestinvdir, restore=False):

    allstat='allstat'+bestinvdir.split('.')[0]+'.dat'+('.revise' if not restore else '')
    inpinv='inpinv'+bestinvdir.split('.')[1]+'.dat'
    igrdat=bestinvdir.split('.')[3]
    icrustal=bestinvdir.split('.')[2]
    igrid=bestinvdir.split('.')[4]
    grhes='gr.{}.{}.{}.*.hes'.format(icrustal,igrdat,igrid)

    for grhes in sorted(glob.glob(os.path.join(workdir,'greens',grhes))):
        isource= grhes.split('.')[-2]

        command='{} {} {} {} {} {} {} {} {}\n'.format(os.path.join(os.getcwd(),_cfg['Inversion']['ExePath']),
                os.path.join('..', '..','allstat', allstat),
                os.path.join('..', '..','inpinv', inpinv),
                os.path.join('..', '..','grdat', 'grdat'+igrdat+'.hed'),
                os.path.join('..', '..','greens', os.path.basename(grhes)),
                os.path.join('..', '..','crustals', 'crustal'+icrustal+'.dat'),
                os.path.join('..', '..','sources', 'grid'+igrid, 'source'+isource+'.dat'),
                os.path.join('..', '..','station.dat'),
                os.path.join('..', '..','raw', igrdat))
        #print(command)
        invdir='{}.{}.{}.{}.{}.{}'.format(allstat.split('.')[0][7:], 
               inpinv.split('.')[0][6:], icrustal, igrdat, igrid, isource)

        proc=subprocess.Popen(command, cwd=os.path.join(os.path.join(workdir,'inversions'), invdir), shell=True, universal_newlines=True,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out,err=proc.communicate()
        logger.info(err)
        logger.info(out)



@config.time
def gatherResults(cfg=None, evt=None, workdir=None, bestinvdir=None, revise=False):
    """
    Retrieve all the necessary info derived from Moment Tensor calculation
    """
    if revise:
        config.evt=evt
        config.cfg=cfg
        config.workdir=workdir
        config.bestinvdir=bestinvdir
        config.revise=True
    else:
        workdir=config.workdir
        evt=config.evt

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
        config.bestinvdir=inversions[corr.index(max(corr))]

    with open(os.path.join(workdir,'inversions',config.bestinvdir,'inv1.dat')) as _:
        content=_.readlines()

    # get tl value that was used for the best inversion
    with open(os.path.join(workdir,'grdat','grdat'+config.bestinvdir.split('.')[3]+\
    '.hed'),'r') as f:
        tl=float(f.readlines()[2].split('=')[1])

    config.besttl=tl

    # number of max sources per source file
    sourcefiles=sorted(glob.glob(os.path.join(workdir,'sources','grid'+\
                config.bestinvdir.split('.')[-2]+'/source*.dat')))

    with open(sourcefiles[0],'r') as f:
        maxsources=len(f.readlines())

    # read all source files from all inversions
    srctext=[]
    for i,src in enumerate(sourcefiles):
        with open(src,'r') as f:
            srctext+=f.readlines()

    # read all inv1.dat files from all inversions of the best inversion
    invtext=[]
    for i,inv1 in enumerate(sorted(glob.glob(os.path.join(workdir,\
    'inversions','.'.join(config.bestinvdir.split('.')[:-1])+'.*/inv1.dat')))):
        with open(inv1,'r') as f:
            temp=f.readlines()[3:]

        for i,line in enumerate(temp):
            if line.startswith(' Selected source position for subevent'):
                invtext+=temp[:i-1]
                break

    # merge source and inv1.dat files from all inversions
    # of the best inversion to one list
    # convert all elements to floats
    solutions=[]
    for i, source in enumerate(srctext):
        #print(source.split()+invtext[i].split()[1:])
        solutions+=[list(map(float,source.split()+invtext[i].split()[1:]))]

    # find best
    # ['3', '-2.0000', '0.0000', '8.6000', '89', '0.550889', '0.1180E+15', 
    # '84.198', '336', '68', '-162', '239', '74', '-22']
    best=max(solutions, key=lambda x: x[5])

    # read all time inversions from corr00.dat
    # 1   -4.0400    0.5474        115    46   -81        283    43   -99   
    # 81.22   -0.00 0.128037E-07 0.857937E+16
    timetext=[]
    for corr in sorted(glob.glob(os.path.join(workdir,\
    'inversions','.'.join(config.bestinvdir.split('.')[:-1])+'.*/corr00.dat'))):
        with open(corr,'r') as f:
            text=f.readlines()[2:]
        pos=int(os.path.basename(os.path.dirname(corr)).split('.')[-1])*maxsources
        text=list(map(lambda x: list(map(float,x.split())),text))
        timetext+=list(map(lambda x: [x[0]+pos]+x[1:],text))

    # calculate quality metrics
    # stvar
    thres=0.9*best[5]
    threstimetext=list(filter(lambda x: x[2]>=thres, timetext))
    stvar=len(threstimetext)/len(timetext)

    # fmvar
    threstimetext=list(map(lambda x: [*best[8:11]]+ [*x[6:9]], threstimetext))
    with multiprocessing.Pool() as p:
       res=p.map(kaganCalc, threstimetext)
    fmvar=np.mean(res)

    # merge all time inversions from best inversions to one dir
    srcs=[sol[0] for sol in solutions]
 
    # keep all correlations with the best's x,y
    correlations=list(filter(lambda x: x if \
    solutions[srcs.index(x[0])][1:3]==best[1:3] else None, timetext))
    # attach depth value instead of src id
    correlations=list(map(lambda x: solutions[srcs.index(x[0])][0:4]+x[1:], \
    correlations))

    # save all best time solutions to one file
    with open(os.path.join(workdir,'output',('solutions' if not revise else 'solutions.revise')),'w') as f:
        f.writelines('{}\n'.format(x)[1:-2]+'\n' for x in sorted(solutions))

    config.solutions=solutions

    # save all time solutions to one file
    with open(os.path.join(workdir,'output',('correlations' if not revise else 'correlations.revise')),'w') as f:
        f.writelines('{}\n'.format(x)[1:-2]+'\n' for x in sorted(correlations))

    config.correlations=correlations

    # save only the best inversion
    # create ObsPy objects based on QuakeML standard
    fm=FocalMechanism()
    mt=MomentTensor()
    org=Origin()
    mag=Magnitude()

    org.time=event.getOrigin(config.cfg,evt,config.cfg['Watcher']['Historical']).time+float(best[4])*(tl/8192)
    _lat, _lon=cart2earth(event.getOrigin(config.cfg,evt,config.cfg['Watcher']['Historical']).latitude,
                          event.getOrigin(config.cfg,evt,config.cfg['Watcher']['Historical']).longitude, best[1], best[2])
    org.latitude=_lat
    org.longitude=_lon
    org.depth= best[3]*1000 # meters
    org.depth_type='from moment tensor inversion'
    org.time_fixed=False
    org.epicenter_fixed=False
    org.origin_type='centroid'
    org.creation_info=CreationInfo(agency_id=config.cfg['Citation']['Agency'], \
                     author=config.cfg['Citation']['Author'], \
                     version=config.cfg['Citation']['Version'], \
                     creation_time=UTCDateTime.now())
    mag.origin_id=org.resource_id
    mag.magnitude_type='Mw'
    mag.creation_info=org.creation_info
    mt.creation_info=org.creation_info
    mt.derived_origin_id=org.resource_id
    mt.moment_magnitude_id=mag.resource_id
    mt.category='regional'
    mt.inversion_type='zero trace' # aka deviatoric
    mt.iso=0
    org.evaluation_mode='automatic'
    org.evaluation_status=('preliminary' if not revise else 'reviewed')

    # fill Focal Mechanism and Moment Tensor values from best inv1.dat file
    for i,line in enumerate(content):
            if line.startswith(' SINGULAR values, incl. vardat'):
                minsn, maxsn, conum = content[i+1].split()

            elif line.startswith(' moment (Nm)'):
                mag.mag=float(content[i+1].split()[2])
                mt.scalar_moment=float(line.split()[2])
                mt.double_couple=float(content[i+3].split()[3])/100.0
                mt.clvd=float(content[i+4].split()[3])/100.0
                nd=content[i+5].split() + content[i+6].split()
                fm.nodal_planes=NodalPlanes(nodal_plane_1=NodalPlane(\
                                strike=float(nd[1]), dip=float(nd[2]), 
                                rake=float(nd[3])), \
                                nodal_plane_2=NodalPlane(strike=float(nd[5]), \
                                dip=float(nd[6]), rake=float(nd[7])))
                ax=content[i+7].split()+content[i+8].split()+content[i+9].split()
                fm.principal_axes=PrincipalAxes(p_axis=Axis(azimuth=float(ax[4]),\
                                  plunge=float(ax[5])), \
                                  t_axis=Axis(azimuth=float(ax[10]), \
                                  plunge=float(ax[11])), \
                                  n_axis=Axis(azimuth=float(ax[16]), \
                                  plunge=float(ax[17])))
             
            elif line.startswith(' varred='):
                mt.variance=float(line.split()[1])
                mt.variance_reduction=float(line.split()[1])*100

    # open best corr inv3.dat file and get Tensor 6 values
    with open(os.path.join(workdir,'inversions', config.bestinvdir, 'inv3.dat')) as _:
        line=list(map(float,_.readlines()[0].split()))
    mt.tensor=Tensor(m_rr=float(line[2]), m_tt=float(line[3]), m_pp=float(line[4]),
                     m_rt=float(line[5]), m_rp=float(line[6]), m_tp=float(line[7]))

    # read best allstat.dat
    with open(os.path.join(workdir,'allstat','allstat'+ \
    config.bestinvdir.split('.')[0]+'.dat'+('.revise' if revise else '')), 'r') as _f:
        text=_f.readlines()
    mt.data_used=[DataUsed(wave_type='combined', station_count=len(text), \
                          component_count=sum([int(comp) for sta in text \
                          for comp in sta.split()[2:5]]), \
                          shortest_period=1/float(text[0].split()[-1]), \
                          longest_period=1/float(text[0].split()[-4]))]

    # add waveform info from Stream object
    st=read(os.path.join(workdir,'streams_corrected.mseed'), headonly=True)

    # read station, priority, distance, azimuth and back azimuth
    with open(os.path.join(workdir,'locations'), 'r') as f:
        stationinfo=f.readlines()
   
    # filter station names to those that are being used in the inversion
    text2=list(filter(lambda x: int(x.split()[1]) and (int(x.split()[2]) \
    or int(x.split()[3]) or int(x.split()[4])), text))
    # get only the station name
    text2=[_[0] for _ in text2]

    stationinfo=list(map(eval,stationinfo))

    # filter stationsinfo to those that are being used in the inversion
    stationinfo2=list(filter(lambda x: x[0] not in text2, stationinfo))

    # get all distances together only from used stations
    dist=[float(_[4][0])/1000 for _ in stationinfo2]
    azm=np.array([float(_[4][1]) for _ in stationinfo2])

    # find station azimuthal gap
    azm.sort()
    # circular shift and substract
    azm=azm-np.roll(azm,1)
    mag.azimuthal_gap=max(np.where(azm<=0, 360+azm, azm))

    components = AttribDict()
    components.namespace = 'custom'
    components.value = AttribDict()
    i=1
    # for all stations (used and unused)
    for line in text:
        w=[0,0,0]
        tr=st.select(station=line.split()[0])[0]
        stainfo=[_ for _ in stationinfo if _[0]==tr.stats.station][0]
        sta, sw, w[0], w[1], w[2], *freqs=line.split()
        tr.stats['frequencies']=freqs
        tr.stats['distance']=round(float(stainfo[4][0])/1000,1)
        tr.stats['latitude']=float(stainfo[2])
        tr.stats['longitude']=float(stainfo[3])
        file=os.path.join(workdir,'inversions',config.bestinvdir,sta+'fil.dat')
        obs={'Time':None, 'N':None, 'E':None, 'Z':None}
        syn={'Time':None, 'N':None, 'E':None, 'Z':None}
        # load observed data
        obs['Time'], obs['N'], obs['E'], obs['Z'] = np.loadtxt(file,unpack=True, usecols=[0,1,2,3])
        # load synthetic data
        syn['Time'], syn['N'], syn['E'], syn['Z'] = np.loadtxt(file[:-7]+'syn.dat',unpack=True, usecols=[0,1,2,3])

        for j, orient in enumerate(['N', 'E', 'Z']):
            tr.stats['channel']=tr.stats.channel[:2]+orient
            tr.stats['weight']=int(int(sw) and int(w[j]))
            tr.stats['variance']=calculateVariance(obs[orient], syn[orient], tl) if tr.stats['weight'] else 'None'

            if tr.stats['weight']:
                fm.waveform_id.append(WaveformStreamID(network_code=tr.stats.network, \
                station_code=tr.stats.station, location_code=tr.stats.location, \
                channel_code=tr.stats.channel))

            for element in ['network', 'station', 'location', 'channel', 'latitude', 'longitude', 'variance', 'weight', 'distance', 'frequencies']:
                if element=='network':
                    components.value['component_'+str(i)]=AttribDict()
                    components.value['component_'+str(i)].namespace='custom'
                    components.value['component_'+str(i)].value=AttribDict()
                components.value['component_'+str(i)].value[element]=AttribDict()
                components.value['component_'+str(i)].value[element].namespace='custom'
                components.value['component_'+str(i)].value[element].type='attribute'
                components.value['component_'+str(i)].value[element].value=tr.stats[element]
            i+=1

    fm.extra=AttribDict()
    fm.extra.components=components

    mag.station_count=mt.data_used[0].station_count
    mag.evaluation_mode=org.evaluation_mode
    mag.evaluation_status=org.evaluation_status
    fm.azimuthal_gap=mag.azimuthal_gap
    fm.triggering_origin_id= event.getOrigin(config.cfg,evt,config.cfg['Watcher']['Historical']).resource_id
    fm.evaluation_mode=org.evaluation_mode
    fm.evaluation_status=org.evaluation_status
    fm.creation_info=org.creation_info
    fm.moment_tensor=mt
    org.quality=OriginQuality(used_station_count=mt.data_used[0].station_count, \
                              azimuthal_gap=mag.azimuthal_gap, \
                              minimum_distance=kilometers2degrees(min(dist)), \
                              maximum_distance=kilometers2degrees(max(dist)), \
                              median_distance=kilometers2degrees(median(dist)))

    # define quality value
    if mt.variance >= 0.6 and mt.data_used[0].station_count > 4:
        quality='A'
    elif (mt.variance >= 0.4 and mt.variance < 0.6 and \
    mt.data_used[0].station_count >= 4) or (mt.variance >= 0.7 and \
    (mt.data_used[0].station_count==2 or mt.data_used[0].station_count==3)):
        quality='B'
    elif (mt.variance >= 0.15 and mt.variance < 0.4 and \
    mt.data_used[0].station_count > 4) or (mt.variance >= 0.2 and \
    mt.variance < 0.4 and mt.data_used[0].station_count == 4) or \
    (mt.variance >= 0.2 and mt.variance < 0.7 and \
    mt.data_used[0].station_count == 3) or (mt.variance >= 0.3 and \
    mt.variance < 0.7 and mt.data_used[0].station_count == 2):
        quality='C'
    else:
        quality='D'

    if mt.clvd <= 0.2:
        quality+='1'
    elif mt.clvd > 0.2 and mt.clvd <= 0.5:
        quality+='2'
    elif mt.clvd > 0.5 and mt.clvd <= 0.8:
        quality+='3'
    elif mt.clvd > 0.8:
        quality+='4'

    mt.extra = OrderedDict()
    mt.extra['correlation'] = {'namespace': 'custom', 'value': best[5]}
    mt.extra['quality'] = {'namespace': 'custom', 'value': quality}
    mt.extra['min_singular'] = {'namespace': 'custom', 'value': minsn}
    mt.extra['max_singular'] = {'namespace': 'custom', 'value': maxsn}
    mt.extra['condition_number'] = {'namespace': 'custom', 'value': conum}
    mt.extra['stvar'] = {'namespace': 'custom', 'value': stvar}
    mt.extra['fmvar'] = {'namespace': 'custom', 'value': fmvar}

    if revise:
            # check and keep only ONE (this) revision
        _fm=list(filter(lambda x: x.evaluation_status=='reviewed',evt.focal_mechanisms))
        if _fm:
            evt.focal_mechanisms=list(filter(lambda x: not x.resource_id==_fm[0].resource_id, evt.focal_mechanisms))
            evt.origins=list(filter(lambda x: not x.resource_id==_fm[0].moment_tensor.derived_origin_id,evt.origins))
            evt.magnitudes=list(filter(lambda x: not x.origin_id==_fm[0].moment_tensor.derived_origin_id,evt.magnitudes))

    # this double for loop is needed in order to maintain the components attribdict
    for _ in evt.focal_mechanisms:
        for _i in range(1,len(_.extra.components.value)+1):
            _.extra['components']['value']['component_'+str(_i)]['value']=AttribDict()

    evt.focal_mechanisms.append(fm)
    evt.origins.append(org)
    evt.magnitudes.append(mag)

    #if revise:
    #    if evt.preferred_focal_mechanism().moment_tensor.variance <= mt.variance:
    #        evt.preferred_focal_mechanism_id=fm.resource_id

    #else:
    evt.preferred_focal_mechanism_id=fm.resource_id

    evt.write(os.path.join(workdir,'output','event.xml'), format="QUAKEML")
    evt.write(os.path.join(workdir,'output','event_sc.xml'), format="SC3ML") 

def getGreens():
    """
    All steps required for Greens' Functions calculation
    """
    config.logger.info('Creating sources files based on Grid rules')
    config.logger.info(config.dump(config.gridRules))
    createSources()

    config.logger.info('Selecting crustal model files based on Crustal rules')
    config.logger.info(config.dump(config.crustalRules))
    createCrustals()

    config.logger.info('Creating Greens\' Functions configuration files (grdat.hed)')
    createGrdat()

    config.logger.info('Creating Stations file (station.dat)')
    createStations()

    config.logger.info('Performing Greens\' Functions Computation')
    calculateGreens()

def getInversions():
    """
    All steps required for Inversions calculation
    """
    config.logger.info('Creating Inversions configuration files (inpinv) ' + \
    'based on Time-shift rules')
    config.logger.info(config.dump(config.shiftRules))
    createInpinv()

    config.logger.info('Select Stations for inversion (allstat)')
    createAllstat()

    config.logger.info('Creating raw data files for ISOLA parsing')
    createRaw()

    config.logger.info('Calculating inversions')
    calculateInversions()

