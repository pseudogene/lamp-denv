#
# Copyright 2014-2016, Michaël Bekaert <michael.bekaert@stir.ac.uk>
#
# This file is part of lamp-denv.
#
# lamp-denv is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# lamp-denv is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with lamp-denvs. If not, see <http://www.gnu.org/licenses/>.
#
FROM ubuntu:15.10
MAINTAINER Michael Bekaert <michael.bekaert@stir.ac.uk>

LABEL description="LAMP-DENV Docker" version="1.0" Vendor="Institute of Aquaculture, University of Stirling"

USER root

RUN apt-get update
RUN mkdir /dengue
COPY denv_install.sh /denv_install.sh
COPY class_sequences.pl /dengue/class_sequences.pl
COPY collect_genomevirus.pl /dengue/collect_genomevirus.pl
COPY denv.R /dengue/denv.R
COPY map_lamp.pl /dengue/map_lamp.pl
RUN /bin/bash /denv_install.sh
RUN rm -f /denv_install.sh

WORKDIR /dengue
