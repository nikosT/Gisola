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

# Example script of how to autobuild Gisola's web suite

# cron-job example
# re-create once per year the new year list item
# @yearly /home/user/tools/years.sh > /home/user/tools/years_log

# create new year directory
workdir=/home/user/gisola/src/results/realtime

# remove last slash if exists
workdir=${workdir%*/}

# current UTC year
cur=`date -u +%Y`

echo "mkdir $workdir/$cur"
mkdir $workdir/$cur

# generate home again containing current year in list options
source /home/user/.bash_profile; cd /home/user/gisola/src

for f in $workdir/*; do
    [[ -e $f ]] || continue

    # folder must be year
    year=`basename $f`
    if ! [[ $year =~ ^[0-9]+$ ]] ; then
        continue
    fi

    # folder must NOT be current year
    if [[ "$year" == "$cur" ]] ; then
        continue
    fi

    echo $year
    python3 ./gisola.py --home $year

done

echo 'all'
python3 ./gisola.py --home all
