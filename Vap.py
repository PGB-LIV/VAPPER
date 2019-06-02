"""
 * Copyright 2018 University of Liverpool
 * Author John Heap, Computational Biology Facility, UoL
 * Based on original scripts of Sara Silva Silva Pereira, Institute of Infection and Global Health, UoL
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 """
import os
import sys
#import pandas as pd
#import numpy as np
#import seaborn as sns
#import matplotlib.pyplot as plt
#from matplotlib.mlab import PCA
import Tryp_G
import Tryp_T
import Tryp_V
import Tryp_Multi
import Tryp_V_T


import argparse
#Entry .sort out the arguments
errorString = "**** ERROR ********\n"
pdf = False

# this dictionary is filled out by initialiseParameters() and passed to individual processes and methods
vapperDict = {
    'name':"",
    'species':"",
    'strain':"",
    'directory':"",
    'pdf': "",
    'pathway':"",
    'kmers': "",
    'inslen':"",
    'covcut':"",
    'forward': "",
    'reverse': "",
    'contigs':"",
    'html_file': "",
    'html_resource': "",
    'refFastq': ""
}



"""
Function: initialiseParamters(cargs):
Flushes out the global vapperDict with the provided command line parameters
"""
def initialiseParameters(cargs):
    global vapperDict
    vapperDict['name']=cargs['name']      #Prefix for results directory and files therein. Only positional argument
    vapperDict['species'] = 'T.congolense'  #default
    if cargs['s'] != 'T.congolense':
        vapperDict['species']='T.vivax'     #if not T.congolense then T.vivax

    vapperDict['strain'] = cargs['strain']  #defaults to Tc148
    if cargs['cdir']:
        vapperDict['directory']= cargs['cdir']
        vapperDict['contigs']="multiCon"    #just a token to ensure we go down the contigs pathway
    if cargs['dir']:                  #is there a directory? If so here we hold multipe samples
        vapperDict['directory']=cargs['dir']

    if cargs['p']:
        vapperDict['pdf'] = "PDF_Yes"
    else:
        vapperDict['pdf'] = "PDF_No"

    #set up pathway
    vapperDict['pathway']="Genomic" #default
    if cargs['t']:
        vapperDict['pathway'] = 'Transcriptomic'
        if cargs['ref']:
            vapperDict['refFastq'] = cargs['ref']
        # vapperDict['refFastq'] = "vivax_trans/Trinity.fasta"

    #assembly directives
    vapperDict['kmers']=str(cargs['k'])
    vapperDict['inslen']=str(cargs['i'])
    vapperDict['covcut']=str(cargs['cov'])

    #forward, reverse, or contig file
    if cargs['f']:
        vapperDict['forward']= cargs['f']
    if cargs['r']:
        vapperDict['reverse']= cargs['r']
    if cargs['con']:
        vapperDict['contigs']=cargs['con']

    #output directory and html file naem
    htmlfile = r"./results/" + cargs['name'] + "/" + cargs['name'] + ".html"
    htmldir = r"./results/" + cargs['name'] + "/"
    if not os.path.exists(htmldir):
        os.makedirs(htmldir)
    vapperDict['html_file']=htmlfile
    vapperDict['html_resource']=htmldir


"""
Entry point to Vap.py
"""
#set up command line argument parser
parser = argparse.ArgumentParser(description='Variant Antigen Profiler - the VAPPER.')
parser.add_argument('name', help = "Prefix for results directory and files therein")
parser.add_argument('-s', default="T.congolense", help= "Species: T.congolense (default) or T.vivax")
parser.add_argument('-con', help = "Contigs File (fasta)")
parser.add_argument('-t','-T', action = 'store_true', default = False, help = "Transcriptomic Pathway")
parser.add_argument('-p','-P', action = 'store_true', default = False, help = "Export PDFs to results directory, as well as .png")
parser.add_argument('-ref', default="/data/Reference/Trinity.fa",
                            help="Fastq reference file for T.vivax Transcriptomic pathway. (defaults to Trinity.fa)")
parser.add_argument('-strain',default = "Tc148", help = "strain for Transcriptomic pathway. (defaults to Tc148)")
parser.add_argument('-dir', help = "Directory that holds multiple paired NGS readfiles for analysis")
parser.add_argument('-cdir', help = "Directory that holds multiple pre-assembled contigs (fasta) files for analysis")
parser.add_argument('-f', help = "Forward NGS read file")
parser.add_argument('-r', help = "Reverse NGS Read File")
parser.add_argument('-k', type = int, default = 65, help = 'kmers (default = 65) as used in velvet')
parser.add_argument('-i',type = int, default = 400, help = 'Insert Length (default = 400) as used in velvet')
parser.add_argument('-cov', type = int, default = 5, help = 'Coverage cut off (default = 5) as used in velvet')
cargs = vars(parser.parse_args())  # cargs = list of commnand line arguments
initialiseParameters(cargs)         #flushes out vapperDict with command line arguments
# print (vapperDict)

#
if vapperDict['directory'] != '':
    # MULTIPLE SAMPLES
    print("Multiple samples indicated in directory: %s" % vapperDict['directory'])
    if vapperDict['species'] != "T.congolense":
        print("Assuming species = T.Vivax")
        if vapperDict['contigs'] != "":
            print("Looking for Contig files in %s" % vapperDict['directory'])
            Tryp_Multi.multi_V_Contigs(vapperDict)
        else:
            if vapperDict['pathway'] == 'Transcriptomic':
                print("Looking for paired NFS readfiles in %s for Transcriptomic Analysis" % vapperDict['directory'])
                Tryp_Multi.multi_V_T(vapperDict)
            else:
                print("Looking for paired NFS readfiles in %s" % vapperDict['directory'])
                Tryp_Multi.multi_V_Assembly(vapperDict)
    else:
        print("Assuming species = T.congolense")
        if vapperDict['pathway'] == 'Transcriptomic':
            print('Transcriptomic Pathway: searching for paired transcripts in %s' %vapperDict['directory'])
            Tryp_Multi.multi_T_process(vapperDict)
        else:
            print('T.congolense Genomic Pathway:')
            if vapperDict['contigs'] != "":
                print("Looking for Contig files in %s" % vapperDict['directory'])
                Tryp_Multi.multi_G_contigs(vapperDict)
            else:
                print("Looking for paired NFS readfiles in %s" % vapperDict['directory'])
                Tryp_Multi.multi_G_Assembly(vapperDict)
else:
    #single sample version of the above
    if vapperDict['species'] != "T.congolense":
        print("Assuming species = T.vivax")
        if vapperDict['pathway'] == 'Transcriptomic':
            print('T.vivax Transcriptomic Pathway:')
            print("Using paired NGS readfiles %s and %s" % (vapperDict['forward'], vapperDict['reverse']))
            Tryp_V_T.transcriptomicProcess(vapperDict)
            sys.exit()

        if vapperDict['contigs']!="":
            print("Looking for Contig file: %s" % vapperDict['contigs'])
            Tryp_V.vivax_contigs(vapperDict)
        else:
            print("Assuming full T.Vivax assembly")
            if vapperDict['forward'] == "" or vapperDict['reverse'] == "":
                print(
                    errorString + "For full assembly we require both forward and reverse NGS readfiles to be sepecified\n" + errorString)
                parser.print_help()
                sys.exit()
            print("Assembling paired NGS readfiles %s and %s" %(vapperDict['forward'],vapperDict['reverse']))
            Tryp_V.vivax_assemble(vapperDict)
    else:
        print("Assuming species = T.congolense")
        if vapperDict['pathway'] =='Transcriptomic':
            print('T.congolense Transcriptomic Pathway:')
            print("Using paired NGS readfiles %s and %s" %(vapperDict['forward'],vapperDict['reverse']))
            Tryp_T.transcriptomicProcess(vapperDict)
        else:
            print('T.congolense Genomic Pathway:')
            if vapperDict['contigs'] != "":
                print("Looking for Contig file: %s" % vapperDict['contigs'])
                Tryp_G.contigs(vapperDict)
            else:
                print("Assuming full T.congolense assembly")
                #need forwards and reverse file names
                if vapperDict['forward'] == "" or vapperDict['reverse'] == "":
                    print(
                            errorString + "For full assembly we require both forward and reverse NGS readfiles to be sepecified\n" + errorString)
                    parser.print_help()
                    sys.exit()
                print("Assembling paired NGS readfiles %s and %s" % (vapperDict['forward'], vapperDict['reverse']))
                Tryp_G.assemble(vapperDict)

sys.exit()








