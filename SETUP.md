 
To ensure your system has all the required dependencies and before running 
the VAPPER code for the first time, 
please type: 
source install.sh

install.sh will:
Temporarily add a path to your system PATH variable 
Check for the installation of transeq, if not available it will download, configure, and compile the EMBOSS
suite of packages. 
It wll set up a python virtual environment VAPENV in which any required python packages may be installed.  

This only needs to be done once per installation - after this to set $PATH and the virtual environment again type: source setup.sh


We have provided a small test contig file (to keep file size manageable) 
To test the T.congolense contig pathway - type:
python Vap.py tc_test -con test_data/Tc_contigs.fa

The results will be found in results/tc_test.


