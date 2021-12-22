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

# Example Python script to apply an action to all events in the MT directory (e.g. fix the event's id)

import shutil,sys
import glob, os
import web
import config,event
from obspy import read_events

year=sys.argv[1]

cfg={}

cfg['WorkDir']='/home/user/gisola/src/results/realtime'


config.cfg=config.read('./config.yaml')


for folder in sorted(os.listdir(os.path.join(cfg['WorkDir'],str(year))), reverse=True):
    # if not dir, ignore it
    if os.path.isfile(os.path.join(cfg['WorkDir'],str(year),folder)):
        continue

    try:
        config.outputdir=None
        dirs=[x for x in glob.glob(os.path.join(cfg['WorkDir'],str(year),folder, '*')) if os.path.isfile(os.path.join(x,'output','event.xml'))]
        rundir=max(dirs, key=os.path.getctime)

        config.evt=read_events(os.path.join(rundir,'output','event.xml'))[0]

        config.outputdir=os.path.join(rundir,'output')

        if os.path.isfile(os.path.join(rundir,'output','mt.revise.txt')):
            print('revise')

        print(rundir)

    except Exception as e:
        print(e)
        continue

