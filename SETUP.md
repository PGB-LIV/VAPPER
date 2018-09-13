 
To ensure your system has all the required dependencies and before running the VAPPER code for the first time, please type: source install.sh
You will need sudo priveleges
install.sh will:
Temporarily add a path to your system PATH variable 
Check for the installation of three programs and attempt to install them if absent. 
These are: virtualenv, transeq, bowtie2
Set up a python virtual environment VAPENV in which any required python packages may be installed.  

This only needs to be done once per installation - after this to set $PATH and the virtual environment again type: source setup.sh


We have provided some small test contigs (to keep file size manageable) 
To test the T.congolense contig pathway - type
python Vap.py tc_test -con test_data/Tc_contigs.fa
To test the T.vivax contig pathway
python Vap.py tv_test -s T.vivax -con test_data/Tv493_contigs.fa

The results will be found in results/tc_test and results/tv_test respectively
 


