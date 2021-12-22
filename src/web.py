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

import os, shutil, glob
from obspy import UTCDateTime, read_events
from obspy.geodetics.base import degrees2kilometers
from jinja2 import Template
import smtplib, ssl
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

# local lib
import config, event

def getYears(cfg):

    l=sorted([str(x) for x in os.listdir(cfg['WorkDir']) if not os.path.isfile(os.path.join(cfg['WorkDir'],x))], reverse=True)
    try:
        l.remove('all')
    except:
        pass

    return l + ['all']

def getList(cfg, year, all=False):
    mtlist=[]

    for folder in sorted(os.listdir(os.path.join(cfg['WorkDir'],str(year))), reverse=True):
        # if not dir, ignore it
        if os.path.isfile(os.path.join(cfg['WorkDir'],str(year),folder)):
            continue

        try:
            dirs=[x for x in glob.glob(os.path.join(cfg['WorkDir'],str(year),folder, '*')) if os.path.isfile(os.path.join(x,'output','event.xml'))]
            rundir=max(dirs, key=os.path.getctime)

            evt=read_events(os.path.join(rundir,'output','event.xml'))[0]

            fm=evt.preferred_focal_mechanism()

            # hypo
            horg=[_ for _ in evt.origins if _.resource_id==fm.triggering_origin_id][0]

            org=[_ for _ in evt.origins if _.resource_id==fm.moment_tensor.derived_origin_id][0]
            mag=[_ for _ in evt.magnitudes if _.origin_id==org.resource_id][0]


            beachball=os.path.join(folder,os.path.basename(rundir), 'output', ('beachball.png' if org.evaluation_status=='preliminary' else 'beachball.revise.png'))
            atlas=os.path.join(folder,os.path.basename(rundir), 'output', ('map.png' if org.evaluation_status=='preliminary' else 'map.revise.png'))
            index=os.path.join(folder,os.path.basename(rundir), 'output', 'index.html')

            # for all case
            if all:
                beachball=os.path.join('..',str(year),beachball)
                atlas=os.path.join('..',str(year),atlas)
                index=os.path.join('..',str(year),index)

            mtlist.append('{},{:.4f},{:.4f},{:d},{:.1f},{:.3e},{},{:d},{:d},{:d},{},{},{},{},{:.4f},{:.4f}\n'.format(str(org.time).split('.')[0].replace('T', ' '), org.latitude, org.longitude, int(round(org.depth/1000)), mag.mag, fm.moment_tensor.scalar_moment, beachball, int(fm.nodal_planes.nodal_plane_1.strike), int(fm.nodal_planes.nodal_plane_1.dip), int(fm.nodal_planes.nodal_plane_1.rake), atlas, fm.moment_tensor.extra['quality'].value, index,org.evaluation_status,horg.latitude,horg.longitude))
        except:
            continue

    # write mtlist to file
    #with open(str(year), 'w') as f:
    #    f.write(''.join([mt for mt in mtlist]))
    mtlist.sort(reverse=True)
    return mtlist

def renderHome(year=None):

    # Read file content and create template object
    template = Template(open(os.path.join('web','templates','home.html'), 'r').read())

    years=getYears(config.cfg)

    if not year:
        org=event.getOrigin(config.cfg,config.evt,config.cfg['Watcher']['Historical'])
        year=str(org.time)[:4]

    mtlist=getList(config.cfg, year)

    # Render HTML Template String
    out = template.render(years=years, current=str(year), lines=mtlist, agency=config.cfg['Citation']['Agency'], website=config.cfg['Citation']['Website'], threshold=config.cfg['Citation']['Quality'], list='.?show=1')

    # output html year
    with open(os.path.join(config.cfg['WorkDir'],str(year), 'list_all.html'), 'w') as f:
        f.write(out)

    if config.cfg['Citation']['Logo']:
        shutil.copy(config.cfg['Citation']['Logo'], os.path.join(config.cfg['WorkDir'],str(year),'logo.png'))
 
    # Read file content and create template object
    template = Template(open(os.path.join('web','templates','map.html'), 'r').read())

    # Render HTML Template String
    out = template.render(years=years, current=str(year), lines=mtlist[::-1], agency=config.cfg['Citation']['Agency'], website=config.cfg['Citation']['Website'], threshold=config.cfg['Citation']['Quality'])

    # output html year
    with open(os.path.join(config.cfg['WorkDir'],str(year), 'index.html'), 'w') as f:
        f.write(out)

    # filter based on Quality Threshold
    qualities=['A1','A2','A3','A4','B1','B2','B3','B4','C1','C2','C3','C4','D1','D2','D3','D4']
    qualities=qualities[:qualities.index(config.cfg['Citation']['Quality'])+1]
    mtlist=list(filter(lambda x: (x.split(',')[11] in qualities), mtlist))

    # Read file content and create template object
    template = Template(open(os.path.join('web','templates','home.html'), 'r').read())

    # Render HTML Template String
    out = template.render(years=years, current=str(year), lines=mtlist, agency=config.cfg['Citation']['Agency'], website=config.cfg['Citation']['Website'], threshold=config.cfg['Citation']['Quality'], list='.')

    # output html year
    with open(os.path.join(config.cfg['WorkDir'],str(year), 'list.html'), 'w') as f:
        f.write(out)

    # create years file
    #with open(os.path.join(config.cfg['WorkDir'],'years'), 'w') as f:
    #    f.write(' '.join(years))


    # if given year is the current/latest year
    try:
        if year==years[0]:
            with open(os.path.join(config.cfg['WorkDir'], 'index.html'), 'w') as f:
                f.write("<html><head></head><script>window.location.href = '"+str(year)+"';</script><body></body></html>")
    except:
        pass

def commandRenderHome(cfg=None,year=None):

    # Read file content and create template object
    template = Template(open(os.path.join('web','templates','home.html'), 'r').read())

    years=getYears(cfg)

    if not year:
        org=event.getOrigin(config.cfg,config.evt,cfg['Watcher']['Historical'])
        year=str(org.time)[:4]

    if year=='all':
        mtlist=[]
        try:
            os.makedirs(os.path.join(cfg['WorkDir'],str(year)))
        except:
            pass
        for y in sorted(years, reverse=True)[1:]:
            mtlist+=getList(cfg, y, all=True)
    else:
        mtlist=getList(cfg, year)

    # Render HTML Template String
    out = template.render(years=years, current=str(year), lines=mtlist, agency=cfg['Citation']['Agency'], website=cfg['Citation']['Website'], threshold=cfg['Citation']['Quality'], list='.?show=1')

    # output html year
    with open(os.path.join(cfg['WorkDir'],str(year), 'list_all.html'), 'w') as f:
        f.write(out)

    if cfg['Citation']['Logo']:
        shutil.copy(cfg['Citation']['Logo'], os.path.join(cfg['WorkDir'],str(year),'logo.png'))

    # Read file content and create template object
    template = Template(open(os.path.join('web','templates','map.html'), 'r').read())

    # Render HTML Template String
    out = template.render(years=years, current=str(year), lines=mtlist[::-1], agency=cfg['Citation']['Agency'], website=cfg['Citation']['Website'], threshold=cfg['Citation']['Quality'])

    # output html year
    with open(os.path.join(cfg['WorkDir'],str(year),'index.html'), 'w') as f:
        f.write(out)

    # filter based on Quality Threshold
    qualities=['A1','A2','A3','A4','B1','B2','B3','B4','C1','C2','C3','C4','D1','D2','D3','D4']
    qualities=qualities[:qualities.index(cfg['Citation']['Quality'])+1]
    mtlist=list(filter(lambda x: (x.split(',')[11] in qualities), mtlist))

    # Read file content and create template object
    template = Template(open(os.path.join('web','templates','home.html'), 'r').read())

    # Render HTML Template String
    out = template.render(years=years, current=str(year), lines=mtlist, agency=cfg['Citation']['Agency'], website=cfg['Citation']['Website'], threshold=cfg['Citation']['Quality'], list='.')

    # output html year
    with open(os.path.join(cfg['WorkDir'], str(year), 'list.html'), 'w') as f:
        f.write(out)

    # create years file
    #with open(os.path.join(cfg['WorkDir'],'years'), 'w') as f:
    #    f.write(' '.join(years))

    # if given year is the current/latest year
    try:
        if year==years[0]:
            with open(os.path.join(cfg['WorkDir'], 'index.html'), 'w') as f:
                f.write("<html><head></head><script>window.location.href = '"+str(year)+"';</script><body></body></html>")
    except:
        pass

def renderEvent():

    fm=event.getFocalMechanism(config.evt)
    # hypo
    org=[_ for _ in config.evt.origins if _.resource_id==fm.triggering_origin_id][0]
    mag=[_ for _ in config.evt.magnitudes if _.origin_id==org.resource_id][0]
    # mt
    org2=[_ for _ in config.evt.origins if _.resource_id==fm.moment_tensor.derived_origin_id][0]
    mag2=[_ for _ in config.evt.magnitudes if _.origin_id==org2.resource_id][0]

    # Read file content and create template object
    template = Template(open(os.path.join('web','templates','index.html'), 'r').read())

    components=fm.extra.components.value

    stations=[]
    for i in range(1,len(components)+1,3):
        stations.append('{} {} {} {} {:d}'.format(
        components['component_'+str(i)]['attrib']['{custom}station'], 
        components['component_'+str(i)]['attrib']['{custom}variance'] if not components['component_'+str(i)]['attrib']['{custom}variance']=='None' else '-', 
        components['component_'+str(i+1)]['attrib']['{custom}variance'] if not components['component_'+str(i+1)]['attrib']['{custom}variance']=='None' else '-', 
        components['component_'+str(i+2)]['attrib']['{custom}variance'] if not components['component_'+str(i+2)]['attrib']['{custom}variance']=='None' else '-', 
        int(round(float(components['component_'+str(i)]['attrib']['{custom}distance']),0))))

    # Render HTML Template String
    out = template.render(home=os.path.join('..','..','..'),
                          agency=config.cfg['Citation']['Agency'], 
                          logo=os.path.join('..','..','..','logo.png'),
                          website=config.cfg['Citation']['Website'],
                          quakeml='event.xml',
                          txt='mt.txt' if org2.evaluation_status=='preliminary' else 'mt.revise.txt',
                          emsc='emsc.txt' if org2.evaluation_status=='preliminary' else 'emsc.revise.txt',
                          beachball='beachball.png' if org2.evaluation_status=='preliminary' else 'beachball.revise.png',
                          atlas='map.png' if org2.evaluation_status=='preliminary' else 'map.revise.png',
                          misfit='misfit.png' if org2.evaluation_status=='preliminary' else 'misfit.revise.png',
                          streams='streams.png' if org2.evaluation_status=='preliminary' else 'streams.revise.png',
                          top='top.png' if org2.evaluation_status=='preliminary' else 'top.revise.png',
                          time='time.png' if org2.evaluation_status=='preliminary' else 'time.revise.png',
                          northeast='northeast.png' if org2.evaluation_status=='preliminary' else 'northeast.revise.png',
                          datetime=org2.time.strftime("%Y-%m-%d %H:%M:%S.%f")[:-4],
                          latitude='{:.4f}'.format(org2.latitude), 
                          longitude='{:.4f}'.format(org2.longitude),
                          depth='{:d}'.format(int(round(org2.depth/1000,0))),
                          cdepth='{:.1f}'.format(round(org2.depth/1000,1)),
                          mw='{:.1f}'.format(round(mag2.mag,1)),
                          mag_type='{}'.format(mag2.magnitude_type),
                          evt_id=os.path.basename(str(config.evt.resource_id)),
                          evt_mode=fm.evaluation_mode,
                          evt_status=fm.evaluation_status,
                          corr='{:.2f}'.format(float(fm.moment_tensor.extra['correlation'].value)),
                          dc='{:.2f}'.format(fm.moment_tensor.double_couple*100),
                          min_sin='{:.3e}'.format(float(fm.moment_tensor.extra['min_singular'].value)),
                          max_sin='{:.3e}'.format(float(fm.moment_tensor.extra['max_singular'].value)),
                          var='{:.2f}'.format(fm.moment_tensor.variance),
                          clvd='{:.2f}'.format(fm.moment_tensor.clvd*100),
                          con_num='{:.2f}'.format(float(fm.moment_tensor.extra['condition_number'].value)),
                          quality=fm.moment_tensor.extra['quality'].value,
                          iso='{:.2f}'.format(fm.moment_tensor.iso*100),
                          stvar='{:.2f}'.format(float(fm.moment_tensor.extra['stvar'].value)),
                          fmvar='{:.2f}'.format(float(fm.moment_tensor.extra['fmvar'].value)),
                          hdatetime=org.time.strftime("%Y-%m-%d %H:%M:%S.%f")[:-4],
                          hlatitude='{:.4f}'.format(org.latitude), 
                          hlongitude='{:.4f}'.format(org.longitude),
                          hdepth='{:.1f}'.format(round(org.depth/1000,1)),
                          hmag='{:.1f}'.format(mag.mag),
                          hmag_type='{}'.format(mag.magnitude_type),
                          sta_count=org2.quality.used_station_count,
                          azm_gap='{:.2f}'.format(org2.quality.azimuthal_gap),
                          min_dist='{:.1f}'.format(degrees2kilometers(org2.quality.minimum_distance)),
                          max_dist='{:.1f}'.format(degrees2kilometers(org2.quality.maximum_distance)),
                          stations=stations,
                          frqs=components['component_'+str(1)]['attrib']['{custom}frequencies'][1:-1],
                          nodals=fm.nodal_planes,
                          axes=fm.principal_axes,
                          mrr='{:.3e}'.format(fm.moment_tensor.tensor.m_rr),
                          mtt='{:.3e}'.format(fm.moment_tensor.tensor.m_tt),
                          mpp='{:.3e}'.format(fm.moment_tensor.tensor.m_pp),
                          mrt='{:.3e}'.format(fm.moment_tensor.tensor.m_rt),
                          mrp='{:.3e}'.format(fm.moment_tensor.tensor.m_rp),
                          mtp='{:.3e}'.format(fm.moment_tensor.tensor.m_tp),
                          mo='{:.3e}'.format(fm.moment_tensor.scalar_moment),)

    # output html year
    with open(os.path.join(config.outputdir, 'index.html'), 'w') as f:
        f.write(out)

def email():

    if not config.cfg['Notification']['Email']['Smtp']:
        return 

    fm=config.evt.preferred_focal_mechanism()
    # mt
    org2=[_ for _ in config.evt.origins if _.resource_id==fm.moment_tensor.derived_origin_id][0]
    mag2=[_ for _ in config.evt.magnitudes if _.origin_id==org2.resource_id][0]

    msg=MIMEMultipart("alternative")
    msg["From"] = config.cfg['Notification']['Email']['Sender']

    msg["Subject"] = '{} || {} {} || {}, {} || {} km || {}'.format( \
                     org2.time.strftime("%Y-%m-%d %H:%M:%S.%f")[:-4], \
                     '{:.1f}'.format(round(mag2.mag,1)), \
                     '{}'.format(mag2.magnitude_type), \
                     '{:.4f}'.format(org2.latitude), \
                     '{:.4f}'.format(org2.longitude), \
                     '{:d}'.format(int(round(org2.depth/1000,0))),\
                      fm.moment_tensor.extra['quality'].value)


    html = """<html><body><h4 style="margin: 0px; color: #000000; line-height: 140%; text-align: center; word-wrap: break-word; font-weight: normal; font-family: 'Montserrat',sans-serif; font-size: 22px;">{}<br><br><a href="{}" target="_self"><span style="display:block;padding:10px 20px;line-height:120%;"><span style="font-size: 14px; line-height: 16.8px;"><img src="https://raw.githubusercontent.com/nikosT/Gisola/main/material/icons/map.png" alt="Map" title="Map" width="32"></span></span></a><a href="{}" target="_self" style="box-sizing: border-box;display: inline-block;font-family:'Montserrat',sans-serif;text-decoration: none;-webkit-text-size-adjust: none;text-align: center;color: #FFFFFF; background-color: #3AAEE0; border-radius: 4px; -webkit-border-radius: 4px; -moz-border-radius: 4px; width:auto; max-width:100%; overflow-wrap: break-word; word-break: break-word; word-wrap:break-word; mso-border-alt: none;"><span style="display:block;padding:10px 20px;line-height:120%;"><span style="font-size: 14px; line-height: 16.8px;">More Info</span></span></a></h4><br><br><br><br><br><br><table width="100%" cellpadding="0" cellspacing="0" border="0"> <tr> <td style="padding-right: 0px;padding-left: 0px;" align="center"> <a href="https://github.com/nikosT/Gisola" target="_self"> <img align="center" border="0" src="https://raw.githubusercontent.com/nikosT/Gisola/main/gisola.png" alt="Image" title="Gisola" style="outline: none;text-decoration: none;-ms-interpolation-mode: bicubic;clear: both;display: inline-block !important;border: none;height: auto;float: none;width: 28%;max-width: 246.4px;" width="246.4"/> </a> </td> </tr></table> </td> </tr> </tbody></table><table style="font-family:'Montserrat',sans-serif;" role="presentation" cellpadding="0" cellspacing="0" width="100%" border="0"> <tbody> <tr> <td style="overflow-wrap:break-word;word-break:break-word;padding:20px;font-family:'Montserrat',sans-serif;" align="left"> <div style="color: #7d7d7d; line-height: 140%; text-align: center; word-wrap: break-word;"> <p style="font-size: 14px; line-height: 140%;"><span style="font-size: 12px; line-height: 16.8px;">🄯 {} Gisola - More info <a href="https://github.com/nikosT/Gisola" target="_self">https://github.com/nikosT/Gisola</a><br /></span></p> </div> </td> </tr> </tbody></table></body></html>""".format(msg["Subject"], os.path.dirname(os.path.dirname(os.path.dirname(os.path.join(config.cfg['Notification']['Email']['HostSite'],os.path.relpath(config.outputdir, config.cfg['WorkDir']))))), os.path.join(config.cfg['Notification']['Email']['HostSite'],os.path.relpath(config.outputdir, config.cfg['WorkDir'])),UTCDateTime().now().year)

    # Turn these into plain/html MIMEText objects
    # Add HTML/plain-text parts to MIMEMultipart message
    # The email client will try to render the last part first
    msg.attach(MIMEText(html, "html"))

    with smtplib.SMTP(config.cfg['Notification']['Email']['Smtp'].split(':')[0], config.cfg['Notification']['Email']['Smtp'].split(':')[1]) as server:
        server.ehlo()  # Can be omitted
        server.starttls(context=ssl.create_default_context())
        server.ehlo()  # Can be omitted
        server.login(config.cfg['Notification']['Email']['User'], config.cfg['Notification']['Email']['Pass'])
        for recipient in config.cfg['Notification']['Email']['Recipients'].replace(', ', ',').split(','):
            msg["To"] = recipient
            server.sendmail(config.cfg['Notification']['Email']['Sender'], recipient, msg.as_string())
        server.quit()


