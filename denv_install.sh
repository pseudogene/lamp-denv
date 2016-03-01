#!/bin/bash
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
echo "deb http://cran.rstudio.com/bin/linux/ubuntu wily/" >> /etc/apt/sources.list
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
gpg -a --export E084DAB9 | apt-key add -
apt-get update
DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential gfortran libblas-dev liblapack-dev make git bioperl primer3 wget unzip dos2unix libgd-dev libgd-perl r-base tk tk-dev

perl -MCPAN -e 'my $c = "CPAN::HandleConfig"; $c->load(doit => 1, autoconfig => 1); $c->edit(prerequisites_policy => "follow"); $c->edit(build_requires_install_policy => "yes"); $c->commit'
perl -MCPAN -e 'install Bio::Graphics'

wget http://heanet.dl.sourceforge.net/project/primer3/primer3/2.3.7/primer3-2.3.7.tar.gz
tar xfz primer3-2.3.7.tar.gz
mkdir -p /usr/local/share/primer3/
cp -r primer3-2.3.7/src/primer3_config /usr/local/share/primer3/primer3_config
rm -rf primer3-2.3.7.tar.gz primer3-2.3.7

wget http://bioinfo.unl.edu/downloads/GramAlign3_00.zip
unzip GramAlign3_00.zip 
cd GramAlign3_00/src
make clean
sed -i 's/CFLAGS = -O3 -funroll/CFLAGS = -O2 -funroll/' Makefile
make
mv GramAlign /usr/local/bin/
cd ../..
rm -rf  GramAlign3_00  GramAlign3_00.zip  __MACOSX

Rscript -e "install.packages('adegenet', repos='http://cran.rstudio.com', dependencies = TRUE);"
Rscript -e "install.packages('ape', repos='http://cran.rstudio.com', dependencies = TRUE);"

git clone https://github.com/pseudogene/lava-dna.git
cd lava-dna
perl Makefile.PL
make
make install
dos2unix *.pl
cd ..

chmod 755 /dengue/*.pl
chmod 644 /dengue/*.R
