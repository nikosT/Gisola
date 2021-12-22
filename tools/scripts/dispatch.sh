#!/usr/bin/bash
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

# example of how to dispatch a SC3ML generated by Gisola file into Gisola's SeisComP FDSNWS-event plugin
file="/home/user/event_sc.xml"
docker exec -it gisola-fdsnws bash -c "/home/sysop/seiscomp/bin/seiscomp exec scdispatch -i $file -O merge"

