#!/bin/bash
#VAPPER setup script
#-------------------
echo "Temporarily setting up the PATH variable"
export PATH=$PATH:$PWD/bin
echo "setting up virtual environment"
virtualenv VAPENV
echo "Activating virtual environment"
source ./VAPENV/bin/activate

