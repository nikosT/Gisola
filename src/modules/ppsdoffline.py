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

# Code snippset also needs ppsdoffline plugin (or a similar tool)

import sys, os, inspect, ast
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import config

def check(trace, dirname='modules/PPSD'):
    """
    Read PPSD results.
    trace: Obspy trace
    path: Path with .txt result files 

    Output:
    True/False
    """
    try:
        with open(os.path.join(dirname,trace.get_id()), 'r') as f:
            lines = f.readlines()

        return ast.literal_eval(lines[0].split()[0])
    except:
        config.logger.exception('PPSD OFFLINE check error occurred for ' +trace)
        return False

