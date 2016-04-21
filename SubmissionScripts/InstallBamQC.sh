#!/bin/bash
#
#PBS -l walltime=12:00:00
#

module add languages/java-jdk-1.8.0-66
cd ~/src
curl -O -L http://mirror.catn.com/pub/apache//ant/binaries/apache-ant-1.9.7-bin.zip
unzip apache-ant-1.9.7-bin.zip
ln -s ~/src/apache-ant-1.9.7/bin/ant ~/bin/ant
curl -L -O https://github.com/s-andrews/BamQC/archive/master.zip
unzip master.zip
cd BamQC-master/
ant
chmod 755 bin/bamqc
ln -s ~/src/BamQC-master/bin/bamqc ~/bin/bamqc

