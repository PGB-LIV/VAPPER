#import subprocess
#import re
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
import argparse
#Entry .sort out the arguments
errorString = "**** ERROR ********\n"
pdfExport = False
parser = argparse.ArgumentParser(description='Variant Antigen Profiler - the VAP.')

parser.add_argument('name', help = "Prefix for results directory/files")
parser.add_argument('-s', default= "T.congolense", help = "Species: T.congolense (default) or T.vivax")
parser.add_argument('-con', help = "Contigs File")
parser.add_argument('-t','-T', action = 'store_true', default = False, help = "Transciptomic Pathway")
parser.add_argument('-p','-P', action = 'store_true', default = False, help = "Export PDFs to results directory")
parser.add_argument('-strain',default = "Tc148", help = "strain required for Transcriptomic pathway")
parser.add_argument('-f', help = "Forward NGS read file")
parser.add_argument('-r', help = "Reverse NGS Read File")
parser.add_argument('-k', type = int, default = 65, help = 'kmers (default = 65)')
parser.add_argument('-i',type = int, default = 400, help = 'Insert Length (default = 400)' )
parser.add_argument('-cov', type = int, default = 5, help = 'Coverage cut off default = 5')

cargs = vars(parser.parse_args())        #cargs = list of commnand line arguments
print (cargs)

#we need to do some sanity checking.
htmlfile = r"./results/"+cargs['name']+"/"+cargs['name']+".html"
htmldir = r"./results/"+cargs['name']+"/"
if not os.path.exists(htmldir):
    os.makedirs(htmldir)



#then wrangle the command line arguments to be the same arrangement as the Galaxy server ones.
if cargs['p']:
    cargs['p']= "PDF_Yes"
else:
    cargs['p']="PDF_No"

if cargs['s'] != "T.congolense":
    print("Assuming Vivax")
    if cargs['con']!=None:
        print("Assuming Contig file available")
        arguments = [cargs['name'],(cargs['con']), htmlfile, htmldir,cargs['p']]
        argdict = {'name': 0, 'pdfexport':4 , 'contigs': 1, 'html_file': 2, 'html_resource': 3}
        Tryp_V.vivax_contigs(arguments, argdict)
    else:
        print("Assuming full T.Vivax assembly")
        if cargs['f'] == None or cargs['r'] == None:
            print(
                errorString + "For full assembly we require both forward and reverse NGS readfiles to be sepecified\n" + errorString)
            parser.print_help()
            sys.exit()

        argdict = {'name': 0, 'pdfexport': 8, 'kmers': 3, 'inslen': 4, 'covcut': 5, 'forward': 1, 'reverse': 2,
                   'html_file': 6, 'html_resource':7}
        arguments = [cargs['name'],cargs['f'],cargs['r'],str(cargs['k']),str(cargs['i']),str(cargs['cov']), htmlfile, htmldir,cargs['p']]
        Tryp_V.vivax_assemble(arguments, argdict)
else:
    print("Assuming T.congolense")
    if cargs['t']:
        print('T.congolense Transcriptomic Pathway')
        argdict = {'name': 0, 'pdfexport': 1, 'strain': 2, 'forward': 3, 'reverse': 4, 'html_file': 5,
                   'html_resource': 6}
        arguments = [cargs['name'],cargs['p'],cargs['strain'],cargs['f'],cargs['r'],htmlfile,htmldir]
        Tryp_T.transcriptomicProcess(arguments, argdict)
    else:
        print('T.congolense Genomic Pathway')
        if cargs['con'] != None:
            print("Assuming Contig file available")
            arguments = [cargs['name'], cargs['con'], htmlfile, htmldir, cargs['p']]
            argdict = {'name': 0, 'pdfexport': 4, 'contigs': 1, 'html_file': 2, 'html_resource': 3}
            Tryp_G.contigs(arguments, argdict)
        else:
            print("Assuming full T.congolense assembly")
            #need forwards and reverse file names
            if cargs['f']==None or cargs['r']==None:
                print (errorString+ "For full assembly we require both forward and reverse NGS readfiles to be sepecified\n"+errorString)
                parser.print_help()
                sys.exit()

            argdict = {'name': 0, 'pdfexport': 8, 'kmers': 3, 'inslen': 4, 'covcut': 5, 'forward': 1, 'reverse': 2,
                       'html_file': 6, 'html_resource': 7}
            arguments = [cargs['name'], cargs['f'], cargs['r'], str(cargs['k']), str(cargs['i']), str(cargs['cov']),
                         htmlfile, htmldir, cargs['p']]


            print(arguments)
            Tryp_G.assemble(arguments, argdict)

sys.exit()







#parser.add_argument('heatmapFile')
#parser.add_argument('PCAFile')
#parser.add_argument('devheatmapFile')
#args = parser.parse_args()

#we have numerous parameters....
#hard code it for differnt types?


#For congolense - we have..
"""
arguments = sys.argv
htmldir = arguments[len(arguments)-1]   #last argument is always html_resource
if not os.path.exists(htmldir):
    os.mkdir(htmldir)

if arguments[1] == 'g_assemble':
    argdict = {'name':2, 'pdfexport':3, 'kmers':4,'inslen':5, 'covcut':6, 'forward':7, 'reverse':8, 'html_file':9, 'html_resource':10}
    Tryp_G.assemble(arguments,argdict)
if arguments[1] == 'g_contigs':
    argdict = {'name':2, 'pdfexport':3, 'contigs':4, 'html_file':5, 'html_resource':6}
    Tryp_G.contigs(arguments,argdict)
if arguments[1] == 'transcipt':
    argdict = {'name':2, 'pdfexport': 3, 'strain': 4, 'forward': 5, 'reverse': 6, 'html_file': 7, 'html_resource': 8}
    Tryp_T.transcriptomicProcess(arguments,argdict)
if arguments[1] == 'v_assemble':
    argdict = {'name':2, 'pdfexport':3, 'kmers':4,'inslen':5, 'covcut':6, 'forward':7, 'reverse':8, 'html_file':9, 'html_resource':10}
    Tryp_V.vivax_assemble(arguments,argdict)
if arguments[1] == 'v_contigs':
    argdict = {'name':2, 'pdfexport':3, 'contigs':4, 'html_file':5, 'html_resource':6}
    Tryp_V.vivax_contigs(arguments,argdict)


sys.exit()

"""



