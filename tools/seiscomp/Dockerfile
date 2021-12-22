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

FROM ubuntu:20.04

RUN apt-get update -y && apt-get install -y sudo wget python3 libpython3.8 lsb-release

RUN groupadd sysop

RUN useradd -rm -d /home/sysop -s /bin/bash -g sysop -G sudo -u 1001 sysop -p sysop

RUN passwd -d sysop

USER sysop
WORKDIR /home/sysop

RUN wget 'https://www.seiscomp.de/downloader/seiscomp-4.6.0-ubuntu20.04-x86_64.tar.gz'

RUN wget 'https://www.seiscomp.de/downloader/seiscomp-maps.tar.gz'

RUN wget 'https://www.seiscomp.de/downloader/seiscomp-4.6.0-doc.tar.gz'

RUN tar -xf seiscomp-4.6.0-ubuntu20.04-x86_64.tar.gz && rm seiscomp-4.6.0-ubuntu20.04-x86_64.tar.gz

RUN tar -xf seiscomp-maps.tar.gz && rm seiscomp-maps.tar.gz

RUN tar -xf seiscomp-4.6.0-doc.tar.gz && rm seiscomp-4.6.0-doc.tar.gz

RUN sudo DEBIAN_FRONTEND="noninteractive" apt-get install tzdata -y

RUN sudo apt update -y

RUN sh -c /bin/echo -e "y" | sudo ./seiscomp/bin/seiscomp install-deps base

RUN sh -c /bin/echo -e "y" | sudo ./seiscomp/bin/seiscomp install-deps fdsnws

RUN sh -c /bin/echo -e "y" | sudo ./seiscomp/bin/seiscomp install-deps mariadb-server

RUN mkdir /home/sysop/.seiscomp
RUN sudo chown sysop:sysop -R /home/sysop/seiscomp /home/sysop/.seiscomp
RUN sudo chmod -R 777 /home/sysop/seiscomp /home/sysop/.seiscomp

RUN sed -i 's/localhost/gisola-db/g' /home/sysop/seiscomp/etc/defaults/scmaster.cfg
RUN echo 'serveDataSelect = false\nserveStation = false' >> /home/sysop/seiscomp/etc/defaults/fdsnws.cfg

COPY ./wait-for-it.sh /home/sysop
RUN sudo chmod -R 777 ./wait-for-it.sh

RUN sudo apt-get install python3-pyqt5 -y

RUN echo "\nexport SEISCOMP_ROOT=\"/home/sysop/seiscomp\"\nexport PATH=\"/home/sysop/seiscomp/bin:$PATH\"\nexport LD_LIBRARY_PATH=\"/home/sysop/seiscomp/lib:$LD_LIBRARY_PATH\"\nexport PYTHONPATH=\"/home/sysop/seiscomp/lib/python:$PYTHONPATH\"\nexport MANPATH=\"/home/sysop/seiscomp/share/man:$MANPATH\"\nexport LC_ALL=C\nsource \"/home/sysop/seiscomp/share/shell-completion/seiscomp.bash\"\n" >> /home/sysop/.bashrc

RUN sudo apt install -y openssh-server xauth

RUN sudo mkdir /var/run/sshd \
    && mkdir /home/sysop/.ssh \
    && chmod 700 /home/sysop/.ssh \
    && sudo ssh-keygen -t rsa -N '' -f /home/sysop/.ssh/id_rsa

RUN sudo sed -i "s/^.*PasswordAuthentication.*$/PasswordAuthentication yes/" /etc/ssh/sshd_config
RUN sudo sed -i 's/#PermitRootLogin prohibit-password/PermitRootLogin yes/' /etc/ssh/sshd_config
RUN sudo sed -i "s/^.*X11Forwarding.*$/X11Forwarding yes/" /etc/ssh/sshd_config
RUN sudo sed -i "s/^.*X11UseLocalhost.*$/X11UseLocalhost no/" /etc/ssh/sshd_config
RUN sudo grep "^X11UseLocalhost" /etc/ssh/sshd_config || echo "X11UseLocalhost no" >> /etc/ssh/sshd_config

RUN sudo apt install net-tools tightvncserver novnc -y

ENV DEBIAN_FRONTEND noninteractive

RUN sudo DEBIAN_FRONTEND=noninteractive apt-get install -y keyboard-configuration
RUN sudo DEBIAN_FRONTEND=noninteractive apt-get install -y lightdm

RUN sudo apt install -y xfce4 xfce4-goodies

RUN export USER=sysop

RUN mkdir $HOME/.vnc && echo "example\nexample\n" | vncpasswd > $HOME/.vnc/passwd && chmod 400 $HOME/.vnc/passwd

RUN echo "#!/bin/bash\nxrdb $HOME/.Xresources\nstartxfce4 &\nxfce4-terminal" > $HOME/.vnc/xstartup && chmod +x $HOME/.vnc/xstartup

ENV USER sysop

RUN echo 'alias scolv-gisola="scolv -d mysql://sysop:sysop@gisola-db/seiscomp"' >> ~/.bashrc

RUN wget http://old.kali.org/kali/pool/main/x/xfwm4/xfwm4_4.14.2-1_amd64.deb && sudo dpkg -i *.deb && rm *.deb

# Clean up (baseimage recommended)
RUN sudo apt-get clean && sudo rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
