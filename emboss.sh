#!/bin/bash
basedir=`pwd`
emboss_config_options='--prefix='$basedir
mkdir -p embossdownload
cd embossdownload
wget 'ftp://emboss.open-bio.org/pub/EMBOSS/emboss-latest.tar.gz'
echo 'decompressing files'
sleep 2
gunzip emboss-latest.tar.gz
tar xf emboss-latest.tar
emboss_dir=`ls -dt EMBOSS-*[^z]|head -1`
cd $emboss_dir
echo 'About to configure'
sleep 2
./configure $emboss_config_options
echo 'About to make'
sleep 2
make
echo 'About to install'
sleep 2
make install
cd $basedir




