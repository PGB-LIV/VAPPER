#!/bin/bash
#VAPPER setup script
#-------------------
#Test command is there
iscmd() {
    command -v >&- "$@"
}

echo "Temporarily setting up the PATH variable"
export PATH=$PATH:$PWD/bin
echo "Make binaries executable"
chmod 755 bin/*

#Checking for dependent applications.
iscmd "virtualenv" || { 
            echo $"virtualenv is not found - installing it"
			sudo apt-get install python-virtualenv
			}
iscmd "transeq" || { 
            echo $"transeq is not found - installing it"
			sudo apt-get install emboss
			}
iscmd "bowtie2" || { 
            echo $"bowtie2 is not found - installing it"
			sudo apt-get install bowtie2
			}
echo "setting up virtual environment"
virtualenv VAPENV
echo "Activating virtual environment"
source ./VAPENV/bin/activate
echo "Installing python packages where required"
VAPENV/bin/pip install -r requirements.txt
