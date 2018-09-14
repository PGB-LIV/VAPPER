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


iscmd "transeq" || { 
            echo "transeq is not found - installing emboss"
			echo "this may take a while"
			sleep 2
			source emboss.sh
			}
echo "setting up virtual environment"
virtualenv VAPENV
echo "Activating virtual environment"
source ./VAPENV/bin/activate
echo "Installing python packages where required"
VAPENV/bin/pip install -r requirements.txt
echo install complete
