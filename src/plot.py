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

import logging, glob, os, os.path
import numpy as np
import math, subprocess
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
from obspy import read
from obspy.geodetics.base import degrees2kilometers, kilometers2degrees, gps2dist_azimuth
# for fast run: sudo -H pip3 install geographiclib !
from mpl_toolkits.basemap import Basemap
mpl.use('Agg')
from obspy.imaging.beachball import beach

# import local configuration
import config, event, web
from obspy import Inventory, UTCDateTime
import matplotlib.tri as tri
import matplotlib.transforms as transforms
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable
from obspy.core.event import read_events
import multiprocessing

def emsc(filepath='emsc.txt'):

    filepath=os.path.join(config.outputdir, ('emsc.txt' if not config.revise else 'emsc.revise.txt'))

    ball=os.path.join(config.inversiondir, config.bestinvdir,'dsretc.lst')

    # preferred origin
    fm=event.getFocalMechanism(config.evt)
    #fm=evt.focal_mechanisms[-1] # temp

    org=[org for org in config.evt.origins if org.resource_id==fm.triggering_origin_id][0]
    mag=[mag for mag in config.evt.magnitudes if mag.origin_id==org.resource_id][0]

    org2=[org for org in config.evt.origins if org.resource_id==fm.moment_tensor.derived_origin_id][0]
    mag2=[mag for mag in config.evt.magnitudes if mag.origin_id==org2.resource_id][0]

    f1,f2,f3,f4=fm.extra.components.value['component_1']['attrib']['{custom}frequencies'][1:-1].split('\', \'')

    mm=[fm.moment_tensor.tensor.m_rr, fm.moment_tensor.tensor.m_tt, \
        fm.moment_tensor.tensor.m_pp, fm.moment_tensor.tensor.m_rt, \
        fm.moment_tensor.tensor.m_rp, fm.moment_tensor.tensor.m_tp]

    base=int('{:.3e}'.format(max(list(map(abs,mm)))).split('e')[1])
    mm=list(map(lambda m: float('{:.3e}'.format(m).split('e')[0])*\
       (10**(int('{:.3e}'.format(m).split('e')[1])-base)),mm))

    with open(ball,'r') as _:
        _ball=''.join([line for line in _.readlines()])

    _text=('______________________________________\n\n' \
           '    {}\n' \
           '      {}\n' \
           '______________________________________\n\n' \
           '======================================\n' \
           '======= Moment Tensor Solution =======\n' \
           '======================================\n\n' \
           'Hypocenter Solution (NOA)\n' \
           'Origin Time :  {}\n' \
           'Lat:  {:.4f}     Lon:  {:.4f}\n' \
           'Depth (km) : {:d}\n' \
           'Mw :  {:.1f}\n\n' \
           '======================================\n' \
           'Centroid Solution\n' \
           'Centroid Time :  {:+.2f} (sec) relative to origin time\n' \
           'Centroid Lat:  {:.4f}       Lon:  {:.4f}\n' \
           'Centroid Depth : {:d}\n\n' \
           '======================================\n' \
           'No of Stations:  {:d}    ({})\n' \
           'Freq band (Hz)\n'
           '{}-{}  tapered {}-{} and {}-{}\n' \
           'Variance Reduction (%): {}\n\n' \
           'Moment Tensor (Nm):  Exponent 10**{:d}\n' \
           '  Mrr      Mtt     Mpp\n' \
           '{:.3f}   {:.3f}   {:.3f}\n' \
           '  Mrt      Mrp     Mtp\n' \
           '{:.3f}   {:.3f}   {:.3f}\n' \
           'VOL (%)      : {}\n' \
           'DC (%)       : {:.1f}\n' \
           'CLVD (%)     : {:.1f}\n' \
           'Quality      : {}\n\n' \
           'Computed by: {}\n\n' \
           'Best Double Couple: Mo= {:.3e}  Nm\n' \
           'NP1:   Strike   Dip   Rake\n' \
           '        {:d}      {:d}     {:d}\n' \
           'NP2:   Strike   Dip   Rake\n' \
           '        {:d}      {:d}     {:d}\n' \
           '{}\n\n'
           'Moment Tensor Solution computed using Gisola\n' \
           'https://github.com/nikosT/Gisola' \
           ).format(config.cfg['Citation']['Agency'].split(' - ')[0], \
                    config.cfg['Citation']['Agency'].split(' - ')[1], 
                    org.time.strftime("%Y%m%d %H:%M:%S.%f")[:-4], \
                    org.latitude, org.longitude, int(round(org.depth/1000,0)), \
                    mag2.mag, org2.time-org.time, org2.latitude, org2.longitude,\
                    int(round(org2.depth/1000,0)), \
                    len(list(set([_.station_code for _ in fm.waveform_id]))), \
                    '-'.join(list(set([_.station_code for _ in fm.waveform_id]))),\
                    f2,f3,f1[1:],f2,f3,f4[:-1], int(round(fm.moment_tensor.variance_reduction,0)), \
                    base, *mm[:3], *mm[3:], fm.moment_tensor.iso*100, \
                    fm.moment_tensor.double_couple*100, fm.moment_tensor.clvd*100, \
                    fm.moment_tensor.extra['quality'].value, config.cfg['Citation']['Author'], \
                    fm.moment_tensor.scalar_moment, int(fm.nodal_planes.nodal_plane_1.strike), \
                    int(fm.nodal_planes.nodal_plane_1.dip), \
                    int(fm.nodal_planes.nodal_plane_1.rake), int(fm.nodal_planes.nodal_plane_2.strike), \
                    int(fm.nodal_planes.nodal_plane_2.dip), int(fm.nodal_planes.nodal_plane_2.rake), _ball)

    with open(filepath, 'w') as _:
        _.write(_text)

def text(filepath='mt.txt'):

    filepath=os.path.join(config.outputdir, ('mt.txt' if not config.revise else 'mt.revise.txt'))

    ball=os.path.join(config.inversiondir, config.bestinvdir,'dsretc.lst')

    # preferred origin
    fm=event.getFocalMechanism(config.evt)

    org=[org for org in config.evt.origins if org.resource_id==fm.triggering_origin_id][0]
    mag=[mag for mag in config.evt.magnitudes if mag.origin_id==org.resource_id][0]

    org2=[org for org in config.evt.origins if org.resource_id==fm.moment_tensor.derived_origin_id][0]
    mag2=[mag for mag in config.evt.magnitudes if mag.origin_id==org2.resource_id][0]

    # Hypocenter info
    _text=('{}\n{} Automatic {}\n{} Moment Tensor Solution {}\n\n{}\n' \
          '{:^68}\n{}\nOrigin Time:\t{}\t\t\tMagnitude:\t{:.1f} {}\n\n' \
          'Latitude: {:.4f} N\t\tLongitude: {:.4f} E\t\t' \
          'Depth: {:.1f} km\n\n').format('-'*68,'-'*29,'-'*28,'-'*22,'-'*22,\
          '-'*68, 'Hypocenter','-'*68, org.time.strftime("%Y-%m-%d %H:%M:%S.%f")[:-4],\
          mag.mag,mag.magnitude_type,org.latitude,org.longitude,round(org.depth/1000,0))

    # Centroid info
    _text+=('{}\n{:^68}\n{}\nCentroid Time:\t {}\t\t\tMagnitude:\t{:.1f} {}\n\n' \
          'Latitude: {:.4f} N\t\tLongitude: {:.4f} E\t\tDepth: {:.1f} km\n\n' \
          'Moment: {:.3e} Nm\n\n').format('-'*68, 'Centroid','-'*68, \
          org2.time.strftime("%Y-%m-%d %H:%M:%S.%f")[:-4],\
          mag2.mag,mag2.magnitude_type,org2.latitude, org2.longitude,\
          int(round(org2.depth/1000,0)), fm.moment_tensor.scalar_moment)

    _text+=('{}\nNP1 --\tStrike:\t\t{:d}\t\tDip:\t\t{:d}\t\tRake:\t\t{:d}\n\n' \
    'NP2 --\tStrike:\t\t{:d}\t\tDip:\t\t{:d}\t\tRake:\t\t{:d}\n\n').format('-'*68, \
    int(fm.nodal_planes.nodal_plane_1.strike), int(fm.nodal_planes.nodal_plane_1.dip), \
    int(fm.nodal_planes.nodal_plane_1.rake), int(fm.nodal_planes.nodal_plane_2.strike), \
    int(fm.nodal_planes.nodal_plane_2.dip), int(fm.nodal_planes.nodal_plane_2.rake))

    _text+=('{}\nP-Axis --\tAzimuth:\t\t{:d}\t\tPlunge:\t\t{:d}\n\nT-Axis --\t' \
    'Azimuth:\t\t{:d}\t\tPlunge:\t\t{:d}\n\nB-Axis --\tAzimuth:\t\t{:d}\t\t' \
    'Plunge:\t\t{:d}\n\n').format('-'*68,int(fm.principal_axes.p_axis.azimuth),\
    int(fm.principal_axes.p_axis.plunge),int(fm.principal_axes.t_axis.azimuth),\
    int(fm.principal_axes.t_axis.plunge),int(fm.principal_axes.n_axis.azimuth),\
    int(fm.principal_axes.n_axis.plunge))

    _text+=('{}\nMrr:\t{:.3e}\t\tMtt:\t{:.3e}\t\tMpp:\t{:.3e}\n\nMrt:\t{:.3e}\t\tMrp:\t' \
    '{:.3e}\t\tMtp:\t{:.3e}\n\n').format('-'*68, fm.moment_tensor.tensor.m_rr, \
    fm.moment_tensor.tensor.m_tt, fm.moment_tensor.tensor.m_pp, \
    fm.moment_tensor.tensor.m_rt, fm.moment_tensor.tensor.m_rp, fm.moment_tensor.tensor.m_tp)

    # Bechball
    with open(ball,'r') as _:
        _text+=''.join([' '*10+line for line in _.readlines()[2:]])

    # Quality Metrics
    _text+=('{}\n{:^68}\n{}\nCorrelation:\t{:.3f}\tVariance Reduction:\t{:.3f}\t' \
    'Quality:\t{}\n\nDC(%):\t{:.2f}\t\tCLVD(%):\t{:.2f}\t\tISO(%):\t{:.2f}\n\n' \
    'Min Singular:\t{:.3e}\t\tMax Singular:\t{:.3e}\n\nCondition Number:\t{:.2f}\t' \
    'STVAR:\t{:.2f}\tFMVAR:\t{:.2f}\n\n').format('-'*68, 'Quality Metrics', '-'*68, \
    float(fm.moment_tensor.extra['correlation'].value), fm.moment_tensor.variance, \
    fm.moment_tensor.extra['quality'].value, fm.moment_tensor.double_couple*100, \
    fm.moment_tensor.clvd*100, fm.moment_tensor.iso*100, \
    float(fm.moment_tensor.extra['min_singular'].value), \
    float(fm.moment_tensor.extra['max_singular'].value), \
    float(fm.moment_tensor.extra['condition_number'].value), \
    float(fm.moment_tensor.extra['stvar'].value), \
    float(fm.moment_tensor.extra['fmvar'].value))

    # Used Stations
    org_quality=[org for org in config.evt.origins if org.resource_id==\
                fm.moment_tensor.derived_origin_id][0].quality

    _text+=('{}\n{:^71}\n{}\nNo of Stations:\t{}\t\tAzimuthal Gap (deg):\t{:.3f}\n\n'\
    'Min Distance(km):\t{:.3f}\t\tMax Distance(km):\t{:.3f}\n\n{}\n').format('-'*68, \
    'Used Stations', '-'*68,org_quality.used_station_count, org_quality.azimuthal_gap,\
    degrees2kilometers(org_quality.minimum_distance), degrees2kilometers(org_quality.maximum_distance),'-'*68)

    components=fm.extra.components.value
    _text+='Stat.\tN\t\tE\t\tZ\tDist. (km)\t\tFrequencies\n\n'
    for i in range(1,len(components)+1,3):
        _text+='{}    {}    {}    {}    {}    {}-{} -- {}-{}\n\n'.format(
        components['component_'+str(i)]['attrib']['{custom}station'], 
        components['component_'+str(i)]['attrib']['{custom}variance'], 
        components['component_'+str(i+1)]['attrib']['{custom}variance'], 
        components['component_'+str(i+2)]['attrib']['{custom}variance'], 
        components['component_'+str(i)]['attrib']['{custom}distance'],
        *components['component_'+str(i)]['attrib']['{custom}frequencies'][1:-1]\
        .replace('\'', '').replace(' ', '').split(',') )

    # evaluation info and credits
    _text+=('{}\n{:^68}\n{}\nEvalutation Mode:\t{}\tEvaluation Status:\t' \
    '{}\n\nAuthor: {}\tVersion: {}\n\nAgency: {}\n').format('-'*68, \
    'Evaluation Info', '-'*68, fm.evaluation_mode, fm.evaluation_status, \
    fm.creation_info.author, fm.creation_info.version, fm.creation_info.agency_id)

    with open(filepath, 'w') as _:
        _.write(_text)

def beachball(filepath='beachball.png'):
    """
    Plotting the focal mechanism
    """
    filepath=os.path.join(config.outputdir, ('beachball.png' if not config.revise else 'beachball.revise.png'))

    fm=event.getFocalMechanism(config.evt)

    # radius of the ball
    radius = 100

    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)

    # Moment Tensor
    try:
        nofill=True
        focal1 = beach([fm.moment_tensor.tensor.m_rr,fm.moment_tensor.tensor.m_tt,\
        fm.moment_tensor.tensor.m_pp,fm.moment_tensor.tensor.m_rt,\
        fm.moment_tensor.tensor.m_rp,fm.moment_tensor.tensor.m_tp], xy=(0.0,0.0), \
        width=2*radius,axes=None,alpha=1, facecolor='r', zorder=1)

        ax.add_collection(focal1)

    except:
        # fall back to dc only with color
        nofill=False
    # Planes
    focal2 = beach([fm.nodal_planes.nodal_plane_1.strike,\
    fm.nodal_planes.nodal_plane_1.dip,fm.nodal_planes.nodal_plane_1.rake], \
    nofill=nofill, facecolor='r', xy=(0.0,0.0),axes=None,width=2*radius,zorder=2)

    ax.add_collection(focal2)

    ax.autoscale_view(tight=False, scalex=True, scaley=True)

    # Convert from spherical to cartesian coordinates the poles 
    r = (90 - fm.principal_axes.t_axis.plunge)*radius/90
    yt = np.cos(np.radians(fm.principal_axes.t_axis.azimuth)) * r  
    xt = np.sin(np.radians(fm.principal_axes.t_axis.azimuth)) * r 

    r = (90 - fm.principal_axes.p_axis.plunge)*radius/90
    yp = np.cos(np.radians(fm.principal_axes.p_axis.azimuth)) * r 
    xp = np.sin(np.radians(fm.principal_axes.p_axis.azimuth)) * r + 0.0
  
    # plot the axis
    ax.plot(xt,yt, color='white', marker='o', markersize=9, zorder=3)
    ax.plot(xp,yp, color='black', marker='o', markersize=9, zorder=4)
    ax.axison=False
    plt.axis('scaled')
    ax.set_aspect(1)

    fig.savefig(filepath, transparent=True)

def atlas(filepath='map.png'):
    """
    Plotting a map with epicenter and possible stations according to distance
    and a map with with the beachball
    """
    filepath=os.path.join(config.outputdir, ('map.png' if not config.revise else 'map.revise.png'))

    fm=event.getFocalMechanism(config.evt)

    org=[org for org in config.evt.origins if \
    org.resource_id==fm.moment_tensor.derived_origin_id][0]

    width=degrees2kilometers(org.quality.maximum_distance)*2*1000+10000 # ofset of 10km

    parallels = [round(org.latitude,2), round((org.latitude-org.quality.maximum_distance/2),2),
                 round((org.latitude+org.quality.maximum_distance/2),2)]
    meridians = [round(org.longitude,2), round((org.longitude-org.quality.maximum_distance/2),2),
                 round((org.longitude+org.quality.maximum_distance/2),2)]

    fig = plt.figure(figsize=(15,15))

    # add map
    m = Basemap(projection='laea', lat_0 = org.latitude,
                lon_0 = org.longitude, lat_ts=org.latitude,
                resolution = 'h', width = width, height = width)

    #labels are [left,right,top,bottom]
    pars=m.drawparallels(parallels,labels=[1,1,0,0], color='grey', fontsize=15)
    mers=m.drawmeridians(meridians,labels=[0,0,1,1], color='grey', fontsize=15)
    for _p in pars:
        pars[_p][1][0].set_rotation(-90)
        pars[_p][1][1].set_rotation(90)

    m.drawrivers(color='#7777ff')
    m.drawcoastlines(color='0.2')
    m.drawcountries(color='0.4')
    m.drawmapboundary(fill_color='#7777ff')
    m.fillcontinents(color='#ddaa66', lake_color='#7777ff')
    x,y = m(org.longitude, org.latitude)

    # add beachball
    ax = plt.gca()

    # mt
    try:
        nofill=True
        focal1 = beach([fm.moment_tensor.tensor.m_rr,fm.moment_tensor.tensor.m_tt,\
        fm.moment_tensor.tensor.m_pp,fm.moment_tensor.tensor.m_rt,\
        fm.moment_tensor.tensor.m_rp,fm.moment_tensor.tensor.m_tp], xy=(x,y), \
        width=width/13, axes=None, alpha=1, facecolor='r', zorder=9, linewidth=0.8)

        ax.add_collection(focal1)

    except:
        nofill=False
    # dc
    focal2 = beach([fm.nodal_planes.nodal_plane_1.strike,\
    fm.nodal_planes.nodal_plane_1.dip,fm.nodal_planes.nodal_plane_1.rake],\
    xy=(x,y),nofill=nofill, facecolor='r', axes=None,width=width/13,zorder=10, linewidth=0.8)

    ax.add_collection(focal2)

    # axes
    r = (90 - fm.principal_axes.t_axis.plunge)*width/13/2/90
    yt = np.cos(np.radians(fm.principal_axes.t_axis.azimuth)) * r  
    xt = np.sin(np.radians(fm.principal_axes.t_axis.azimuth)) * r 
    r = (90 - fm.principal_axes.p_axis.plunge)*width/13/2/90
    yp = np.cos(np.radians(fm.principal_axes.p_axis.azimuth)) * r 
    xp = np.sin(np.radians(fm.principal_axes.p_axis.azimuth)) * r + 0.0
    ax.plot(xt+x,yt+y, color='white', marker='o', markersize=4, zorder=11)
    ax.plot(xp+x,yp+y, color='black', marker='o', markersize=4, zorder=12)

    # add stations
    components=fm.extra.components.value
    for i in range(1,len(components)+1,3):
        _sta=components['component_'+str(i)]['attrib']['{custom}station']
        _lat=float(components['component_'+str(i)]['attrib']['{custom}latitude']) 
        _lon=float(components['component_'+str(i)]['attrib']['{custom}longitude'])
        _w1=bool(float(components['component_'+str(i)]['attrib']['{custom}weight']))
        _w2=bool(float(components['component_'+str(i+1)]['attrib']['{custom}weight']))
        _w3=bool(float(components['component_'+str(i+2)]['attrib']['{custom}weight']))
        _weight=_w1 or _w2 or _w3
        _x,_y = m(_lon, _lat)
        # set color to green if station used else red and zorder respectively
        m.scatter(_x, _y, 1500, color='#52D017' if _weight else '#686868', marker='^',\
        zorder=7 if _weight else 6, linewidths=1.5, edgecolor='k' if _weight else '#686868')
        t=plt.text(_x+1800, _y+3000, _sta, family='monospace', fontsize=30, \
        weight='bold', color='k' if _weight else 'r', zorder=8 if _weight else 6)

    # adding sectors
    circle2 = plt.Circle((x,y), radius=width/2,linestyle='--', linewidth=0.8, fill=False, alpha=0.6)
    ax.add_patch(circle2)
    x2=(width/2/2)*math.cos(math.radians(45))
    y2=(width/2/2)*math.sin(math.radians(45))
    plt.axline((x, y-width/2/2), (x, y+width/2/2), linestyle='--', linewidth=0.8, color='k', alpha=0.6)
    plt.axline((x-width/2/2, y), (x+width/2/2, y), linestyle='--', linewidth=0.8, color='k', alpha=0.6)
    plt.axline((-x2-x, -y2-y), (x2+x, y2+y), linestyle='--', linewidth=0.8, color='k', alpha=0.6)
    plt.axline((x2+x, y-y2), (x-x2, y2+y), linestyle='--', linewidth=0.8, color='k', alpha=0.6)
    plt.axline((-x2, -y2), (x2, y2), linestyle='--', linewidth=0.8, color='k', alpha=0.6)

    fig.tight_layout()
    fig.savefig(filepath)

def streams(st=None, filepath='misfit.png'):
    """
    Plotting the misfits of observed and synthetic timeseries
    of each enabled stream for each station
    """
    filepath=os.path.join(config.outputdir, ('streams.png' if not config.revise else 'streams.revise.png'))

    components=event.getFocalMechanism(config.evt).extra.components.value

    orients=['N', 'E', 'Z']
    formatter = FormatStrFormatter('%.1e')
    plt.rcParams.update({'font.size': 10})

    fig, axs = plt.subplots(int(len(components)/3), 3, constrained_layout=False, figsize=(15, int(len(components)/3)))

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.6)

    # if only one station found (corner case)
    if int(len(components)/3) == 1:
        axs=np.array([axs], dtype=object)

    # set title (first row, second column)
    axs[0,1].annotate('Displacement data in meters', xy=(0.5, 2.1),
                          xycoords='axes fraction',
                          fontsize=13,
                          horizontalalignment='center',
                          verticalalignment='top')

    for _i in range(len(components)):
        i=int(_i/3)
        j=_i%3
        net=components['component_'+str(_i+1)]['attrib']['{custom}network']
        sta=components['component_'+str(_i+1)]['attrib']['{custom}station']
        cha=components['component_'+str(_i+1)]['attrib']['{custom}channel']
        var=components['component_'+str(_i+1)]['attrib']['{custom}variance']
        var= None if var=='None' else var
        w=bool(float(components['component_'+str(_i+1)]['attrib']['{custom}weight']))

        # plot data
        if w:
            tr=config.st.select(station=sta, channel=cha)[0]

            axs[i,j].plot(tr.times(), tr.data, color='black')
 
        else:
            axs[i,j].text(0.5, 0.5, 'unused',
                          horizontalalignment='center',
                          verticalalignment='center',
                          fontsize=11, color='grey',
                          transform=axs[i,j].transAxes)

        # for the first column
        if j==0:
            axs[i,j].set_ylabel(net+'\n'+sta+' '+cha[:2], fontsize=11)

        # for the first row
        if i==0:
            axs[i,j].set_title(orients[j], fontsize=11)

        # set at the last row
        if i==int(len(components)/3)-1:
            axs[i,j].set_xlabel('sec', fontsize=9)

        axs[i,j].yaxis.set_major_formatter(formatter)
        axs[i,j].grid()
        axs[i,j].autoscale(enable=True, axis='both', tight=True)
        axs[i,j].tick_params(labelsize=8)
        axs[i,j].locator_params(axis='y', nbins=3)

    fig.savefig(filepath, bbox_inches='tight', pad_inches=0.3)


def misfit(filepath='misfit.png'):
    """
    Plotting the misfits of observed and synthetic timeseries
    of each enabled stream for each station
    """
    filepath=os.path.join(config.outputdir, ('misfit.png' if not config.revise else 'misfit.revise.png'))

    components=event.getFocalMechanism(config.evt).extra.components.value

    orients=['N', 'E', 'Z']
    formatter = FormatStrFormatter('%.1e')
    plt.rcParams.update({'font.size': 10})

    fig, axs = plt.subplots(int(len(components)/3), 3, constrained_layout=False, figsize=(15, int(len(components)/3)))

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.6)

    # if only one station found (corner case)
    if int(len(components)/3) == 1:
        axs=np.array([axs], dtype=object)

    # set title (first row, second column)
    axs[0,1].annotate('Observed vs Synthetic waveforms in meters', xy=(0.5, 2.1),
                          xycoords='axes fraction',
                          fontsize=13,
                          horizontalalignment='center',
                          verticalalignment='top')

    for _i in range(len(components)):
        i=int(_i/3)
        j=_i%3
        net=components['component_'+str(_i+1)]['attrib']['{custom}network']
        sta=components['component_'+str(_i+1)]['attrib']['{custom}station']
        cha=components['component_'+str(_i+1)]['attrib']['{custom}channel']
        var=components['component_'+str(_i+1)]['attrib']['{custom}variance']
        var= None if var=='None' else var
        w=bool(float(components['component_'+str(_i+1)]['attrib']['{custom}weight']))

        # for the first column
        if j==0:
            obs=[[],[],[],[]]
            syn=[[],[],[],[]]
            # open once per station
            # load observed data (Time,N,E,Z)
            obs[0],obs[1],obs[2],obs[3] = np.loadtxt(os.path.join(config.inversiondir, config.bestinvdir,sta+'fil.dat'),unpack=True, usecols=[0,1,2,3])
            # load synthetic data (Time,N,E,Z)
            syn[0],syn[1],syn[2],syn[3] = np.loadtxt(os.path.join(config.inversiondir, config.bestinvdir,sta+'syn.dat'),unpack=True, usecols=[0,1,2,3])

            axs[i,j].set_ylabel(net+'\n'+sta+' '+cha[:2], fontsize=11)

        # for the first row
        if i==0:
            axs[i,j].set_title(orients[j], fontsize=11)

        # set at the last row
        if i==int(len(components)/3)-1:
            axs[i,j].set_xlabel('sec', fontsize=9)

        axs[i,j].plot(obs[0], obs[j+1], color='black' if var else 'grey', linewidth=2, label="observed")
        axs[i,j].plot(syn[0], syn[j+1], color='red' if var else 'grey', linewidth=2, label="synthetic")

        axs[i,j].yaxis.set_major_formatter(formatter)
        axs[i,j].grid()
        axs[i,j].autoscale(enable=True, axis='both', tight=True)
        axs[i,j].tick_params(labelsize=8)
        axs[i,j].locator_params(axis='y', nbins=3)
   
        axs[i,j].annotate(var if var else '', xy=(0.99, 0.95),
                          xycoords='axes fraction',
                          fontsize=10, color='b',
                          horizontalalignment='right',
                          verticalalignment='top')

        # set legend
        if i==0 and j==2:
            axs[i,j].plot([],[], label='observed', color='black')
            axs[i,j].plot([],[], label='synthetic', color='red')
            axs[i,j].plot([],[], label='unused', color='grey')
            axs[i,j].plot([],[], ' ', label='variance reduction')

            handles, labels = axs[i,j].get_legend_handles_labels()
            pos = list(axs[i,j].get_position().bounds)

            l=axs[i,j].legend(handles[2:], labels[2:], bbox_to_anchor=(pos[0]-pos[0]/6, pos[1]+3*pos[3]), fontsize=8)
            l.get_texts()[3].set_color('blue')

    fig.savefig(filepath, bbox_inches='tight', pad_inches=0.3)


def top(solutionspath=None, tl=None, filepath='top.png'):
    """
    Illustrates an example for plotting beachballs and data points on line
    with specific color based on values with a labelled colorbar.
    """
    solutionspath=os.path.join(config.outputdir,('solutions' if not config.revise else 'solutions.revise'))

    tl=config.besttl

    filepath=os.path.join(config.outputdir, ('top.png' if not config.revise else 'top.revise.png'))

    if not tl:
        # get tl value that was used for the best inversion
        with open(os.path.join(config.workdir,'grdat','grdat'+os.path.basename(config.bestinvdir).split('.')[3]+\
        '.hed'),'r') as f:
            tl=float(f.readlines()[2].split('=')[1])

    dt=tl/8192
    if solutionspath:
        with open(solutionspath, 'r') as f:
            solutions=f.readlines()
        # convert to list of lists of floats
        solutions=list(map(lambda x: list(map(float,x.split(','))),solutions))

    else:
        solutions=solutions

    # find best
    # ['3', '-2.0000', '0.0000', '8.6000', '89', '0.550889', '0.1180E+15', 
    # '84.198', '336', '68', '-162', '239', '74', '-22']
    best=max(solutions, key=lambda x: x[5])

    # get all solutions with best's depth
    solutions=list(filter(lambda x: x if x[3]==best[3] else None,solutions))

    x = [sol[1] for sol in solutions]
    y = [sol[2] for sol in solutions]
    time = [sol[4]*dt for sol in solutions]
    corr = [sol[5] for sol in solutions]
    dc = [sol[7] for sol in solutions]
    beachballs=[[*sol[8:11]] for sol in solutions]

    # creates figure
    fig, ax = plt.subplots(figsize=(15,15))

    fig.suptitle('Top cross section correlation plot for best times\nCross section depth: {:.1f} km'.format(best[3]), fontsize=13)
    # creates a grid on the (main) plot
    ax.grid()

    # sets labels
    plt.ylabel('South (km) <-'+' - '*10+'> North (km)')
    plt.xlabel('West (km) <-'+' - '*10+'> East (km)')

    # sets font size
    plt.rcParams.update({'font.size': 10})

    levels=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    if len(list(set([_ for _ in x])))>=3 and len(list(set([_ for _ in y])))>=3:
        cont = plt.tricontour(x, y, corr, len(levels), levels=levels, \
                          linewidths=0.5, colors='k')

        plt.tricontourf(x, y, corr, len(levels), levels=levels, \
                    cmap=plt.cm.cool)

        plt.clabel(cont)

    cm=plt.cm.copper_r # for dc

    # plotting beachballs on specific x-axis and y-axis
    # with a color based on the data_dc values (normalized to 0-1)
    _index = corr.index(max(corr))

    for i in range(len(beachballs)):
        if i == _index:
            # sets best beachball's color value to red
            _color = 'red'
            _width=50
        else:
            # sets beachball's color value
            _color = cm(dc[i]/100.0)
            _width=40

	    # draws beachball
        b = beach([beachballs[i][0], beachballs[i][1],
                       beachballs[i][2]], xy=(x[i], y[i]),
                       width=_width, linewidth=0.5, facecolor=_color)

        # holds the aspect but fixes positioning:
        b.set_transform(transforms.IdentityTransform())

        # brings the all patches to the origin (0, 0).
        for p in b._paths:
            p.vertices -= [x[i], y[i]]

        # uses the offset property of the collection to position the patches
        b.set_offsets((x[i], y[i]))
        b._transOffset = ax.transData

        # annotate the time value
        ax.annotate('{:.2f}'.format(time[i]), xy=(x[i],y[i]+0.45), color='k', #weight='bold', 
                fontsize=12, ha='center', va='center', zorder=10)

  	    # adds beachball to plot
        ax.add_collection(b)

    ax.set_yticks(np.arange(min(y)-1, max(y)+2, 2))
    ax.set_xticks(np.arange(min(x)-1, max(x)+2, 2))

    # set values to colorbar
    norm = colors.Normalize(vmin=0, vmax=100)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="1%", pad=0.05)
    bar=colorbar.ColorbarBase(cax, cmap=plt.cm.copper_r, norm=norm, orientation='horizontal')
    bar.set_label("DC (%)", labelpad=-19)
    cax.xaxis.set_ticks_position("top")

    cax2 = divider.append_axes("right", size="1%", pad=0.05)
    bar2=colorbar.ColorbarBase(cax2, cmap=plt.cm.cool, norm=norm, orientation='vertical')
    bar2.set_label("Correlation (%)", labelpad=-19)

    fig.tight_layout(pad=1.5)

    fig.savefig(filepath)

def northeast(solutionspath=None, tl=None, filepath='northeast.png'):
    """
    Illustrates an example for plotting beachballs and data points on line
    with specific color based on values with a labelled colorbar.
    """
    @ticker.FuncFormatter
    def major_formatter(x, pos):
        return str(-x) if x < 0 else str(x)

    solutionspath=os.path.join(config.outputdir,('solutions' if not config.revise else 'solutions.revise'))
    tl=config.besttl

    filepath=os.path.join(config.outputdir, ('northeast.png' if not config.revise else 'northeast.revise.png'))

    if not tl:
        # get tl value that was used for the best inversion
        with open(os.path.join(workdir ,'grdat','grdat'+os.path.basename(bestinvdir).split('.')[3]+\
        '.hed'),'r') as f:
            tl=float(f.readlines()[2].split('=')[1])

    dt=tl/8192

    if solutionspath:
        with open(solutionspath, 'r') as f:
            solutions=f.readlines()
        # convert to list of lists of floats
        solutions=list(map(lambda x: list(map(float,x.split(','))),solutions))

    else:
        solutions=solutions

    # find best
    # ['3', '-2.0000', '0.0000', '8.6000', '89', '0.550889', '0.1180E+15', 
    # '84.198', '336', '68', '-162', '239', '74', '-22']
    best=max(solutions, key=lambda x: x[5])

    x=[[],[]]
    depth=[[],[]]
    time=[[],[]]
    corr=[[],[]]
    dc=[[],[]]
    beachballs=[[],[]]

    for i in range(2):
        # get all solutions with best's depth
        _solutions=list(filter(lambda x: x if x[i+1]==best[i+1] else None,solutions))
        x[i] = [sol[1-i+1] for sol in _solutions]
        depth[i] = [sol[3] for sol in _solutions]
        time[i] = [sol[4]*dt for sol in _solutions]
        corr[i] = [sol[5] for sol in _solutions]
        dc[i] = [sol[7] for sol in _solutions]
        beachballs[i]=[[*sol[8:11]] for sol in _solutions]

    # creates figure
    fig, axs = plt.subplots(1,2, figsize=(30,15), sharey=False)

    fig.suptitle('South-North and West-East cross section correlation plot for best times', fontsize=13)
    plt.rcParams.update({'font.size': 10})

    label=[['South', 'North'],['West', 'East']]
    levels=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    for i,ax in enumerate(axs):
        # creates a grid on the (main) plot
        ax.grid()

        if len(list(set([_ for _ in x[i]])))>=3 and len(list(set([_ for _ in depth[i]])))>=3:
            cont = ax.tricontour(x[i], -np.array(depth[i]), corr[i], len(levels), \
                                 levels=levels, linewidths=0.5, colors='k')
            ax.tricontourf(x[i], -np.array(depth[i]), corr[i], len(levels),
                           levels=levels, cmap=plt.cm.cool)
            ax.clabel(cont)

        cm=plt.cm.copper_r # for dc

        # plotting beachballs on specific x-axis and y-axis
        # with a color based on the data_dc values (normalized to 0-1)
        _index = corr[i].index(max(corr[i]))
 
        for j in range(len(beachballs[i])):
            if j == _index:
                # sets best beachball's color value to red
                _color = 'red'
                _width=50
            else:
                # sets beachball's color value
                _color = cm(dc[i][j]/100.0)
                _width=40

	        # draws beachball
            b = beach([beachballs[i][j][0], beachballs[i][j][1],
                       beachballs[i][j][2]], xy=(x[i][j],-depth[i][j]),
                       width=_width, linewidth=0.5, facecolor=_color)


            # holds the aspect but fixes positioning:
            b.set_transform(transforms.IdentityTransform())

            # brings the all patches to the origin (0, 0).
            for p in b._paths:
                p.vertices -= [x[i][j], -depth[i][j]]

            # uses the offset property of the collection to position the patches
            b.set_offsets((x[i][j], -depth[i][j]))
            b._transOffset = ax.transData

            # annotate the time value
            ax.annotate('{:.2f}'.format(time[i][j]), xy=(x[i][j]+0.6,-depth[i][j]), color='k', #weight='bold', 
                fontsize=12, ha='center', va='center', zorder=10)

  	        # adds beachball to plot
            ax.add_collection(b)

        # sets labels
        ax.set_xlabel(label[i][0]+' (km) <-'+' - '*10+'> '+label[i][1]+' (km)', fontsize=10)

        # set axis and labels on the axis
        ax.set_yticks(np.arange(-max(depth[i])-1, -min(depth[i])+2, 1))
        ax.set_xticks(np.arange(min(x[i])-1, max(x[i])+2, 2))

        ax.yaxis.set_major_formatter(major_formatter)

        # set values to colorbar
        norm = colors.Normalize(vmin=0, vmax=100)

        # both plots dc colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("top", size="1%", pad=0.05)
        bar=colorbar.ColorbarBase(cax, cmap=plt.cm.copper_r, norm=norm, orientation='horizontal')
        bar.set_label("DC (%)", labelpad=-19)
        cax.xaxis.set_ticks_position("top")

        # only right plot correlation colorbar
        if i==1:
            cax2 = divider.append_axes("right", size="1%", pad=0.05)
            bar2=colorbar.ColorbarBase(cax2, cmap=plt.cm.cool, norm=norm, orientation='vertical')
            bar2.set_label("Correlation (%)", labelpad=-19)

    # only left plot y axis label
    axs[0].set_ylabel('Depth (km)', fontsize=10)

    fig.tight_layout(pad=1.5)

    fig.savefig(filepath)

def time(correlationspath=None, filepath='time.png'):
    """
    Illustrates an example for plotting beachballs and data points on line
    with specific color based on values with a labelled colorbar.
    """
    @ticker.FuncFormatter
    def major_formatter(x, pos):
        return str(-x) if x < 0 else str(x)
    
    correlationspath=os.path.join(config.outputdir,('correlations' if not config.revise else 'correlations.revise'))
    tl=config.besttl

    filepath=os.path.join(config.outputdir, ('time.png' if not config.revise else 'time.revise.png'))

    if correlationspath:
        with open(os.path.join(correlationspath), 'r') as f:
            correlations=f.readlines()
        # convert to list of lists of floats
        correlations=list(map(lambda x: list(map(float,x.split(','))),correlations))

    else:
        correlations=correlations
    # find best
    # 4.0, -4.0, 2.0, 0.5, -4.04, 0.4485, 153.0, 31.0, 48.0, 19.0, 67.0, 112.0,
    # 64.82, -0.0, 1.46047e-08, 2961630000000000.0
    best=max(correlations, key=lambda x: x[5])

    time = [sol[4] for sol in correlations]
    depth = [sol[3] for sol in correlations]
    corr = [sol[5] for sol in correlations]
    dc = [sol[12] for sol in correlations]
    beachballs=[[*sol[6:9]] for sol in correlations]

    # creates figure
    fig, ax = plt.subplots(figsize=(15,15))

    fig.suptitle('Time correlation plot for best position x,y: {:.1f},{:.1f} km'.format(best[1],best[2]), fontsize=13)
    # creates a grid on the (main) plot
    ax.grid()

    # sets labels
    plt.ylabel('Depth (km)')
    plt.xlabel('Time relative to origin time (sec)')

    # sets font size
    plt.rcParams.update({'font.size': 10})

    levels=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    if len(list(set([_ for _ in time])))>=3 and len(list(set([_ for _ in depth])))>=3:
        cont = plt.tricontour(time, -np.array(depth), corr, len(levels), levels=levels, \
                              linewidths=0.5, colors='k')

        plt.tricontourf(time, -np.array(depth), corr, len(levels), levels=levels, \
                        cmap=plt.cm.cool)

        plt.clabel(cont)

    cm=plt.cm.copper_r # for dc

    # plotting beachballs on specific x-axis and y-axis
    # with a color based on the data_dc values (normalized to 0-1)
    _index = corr.index(max(corr))

    for i in range(len(beachballs)):
        if i == _index:
            # sets best beachball's color value to red
            _color = 'red'
            _width=50
        else:
            # sets beachball's color value
            _color = cm(dc[i]/100.0)
            _width=40

	    # draws beachball
        b = beach([beachballs[i][0], beachballs[i][1],
                       beachballs[i][2]], xy=(time[i], -depth[i]),
                       width=_width, linewidth=0.5, facecolor=_color)

        # holds the aspect but fixes positioning:
        b.set_transform(transforms.IdentityTransform())

        # brings the all patches to the origin (0, 0).
        for p in b._paths:
            p.vertices -= [time[i], -depth[i]]

        # uses the offset property of the collection to position the patches
        b.set_offsets((time[i], -depth[i]))
        b._transOffset = ax.transData

  	    # adds beachball to plot
        ax.add_collection(b)

    # set axis and labels on the axis
    stime=sorted(list(set(time)))
    if len(stime)==1:
        diff_time=1
    else:
        diff_time=abs(stime[1]-stime[0])

    ax.set_yticks(np.arange(-max(depth)-1, -min(depth)+2, 1))
    ax.set_xticks(np.arange(min(stime)-diff_time, max(stime)+diff_time, diff_time))

    ax.yaxis.set_major_formatter(major_formatter)

    # set values to colorbar
    norm = colors.Normalize(vmin=0, vmax=100)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="1%", pad=0.05)
    bar=colorbar.ColorbarBase(cax, cmap=plt.cm.copper_r, norm=norm, orientation='horizontal')
    bar.set_label("DC (%)", labelpad=-19)
    cax.xaxis.set_ticks_position("top")

    cax2 = divider.append_axes("right", size="1%", pad=0.05)
    bar2=colorbar.ColorbarBase(cax2, cmap=plt.cm.cool, norm=norm, orientation='vertical')
    bar2.set_label("Correlation (%)", labelpad=-19)

    fig.tight_layout(pad=1.5)

    fig.savefig(filepath)

def caller(func):
    func()

@config.time
def allplots(evt=None, revise=False):

    functions=[atlas, top, northeast, time, misfit, streams, beachball, text, emsc, web.renderEvent, web.email, web.renderHome]

    if revise:
        config.outputdir=os.path.join(config.workdir,'output')
        config.inversiondir=os.path.join(config.workdir,'inversions')
        config.revise=revise

    config.evt=read_events(os.path.join(config.outputdir, 'event.xml'))[0]

    with multiprocessing.Pool() as p:
       res=p.map(caller, functions)

    if config.cfg['Notification']['Command']:

        # set the full path of the best inversion directory
        command=config.cfg['Notification']['Command'].replace('$bestinvdir',config.bestinvdir)
        command=config.cfg['Notification']['Command'].replace('$outputdir',config.outputdir)

        proc=subprocess.Popen(command, \
             shell=True, universal_newlines=True,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out,err=proc.communicate()
        config.logger.info(err)
        config.logger.info(out)

