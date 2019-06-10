"""
 * Copyright 2019 University of Liverpool
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


import subprocess
import pandas as pd
import re
import os
import sys
import shutil
# import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np




# copies the user provided Fasta file to data/reference/file/file.fasta
def uploadUserReferenceFastq(refFastq):
    refBase = os.path.basename(refFastq)
    ref = os.path.splitext(refBase)[0]   # 'mydata/test.fasta' -> 'test'
    dir_path = os.path.dirname(os.path.realpath(__file__))  # directory of this file
    refDir = dir_path + "/data/Reference/" + ref            #propose putting file in '/data/reference/ref/
    if not os.path.isdir(refDir):  # if directory data/Reference/ref doesn't exist
        os.mkdir(refDir)
    refPath = refDir+"/"
    shutil.copy(refFastq, refPath + refBase)  #copy reference file into the directory
    argString = "bowtie2-build " + refPath + refBase+" "+refPath+ref
    print("Building the bowtie2 reference files.")
    subprocess.call(argString, shell=True)
    return

def transcriptMapping(inputname, refFastq, forwardFN, reverseFN):
    # where is our Reference data?
    refBase = os.path.basename(refFastq)
    ref = os.path.splitext(refBase)[0]
    dir_path = os.path.dirname(os.path.realpath(__file__))
    refDir = dir_path + "/data/Reference/" + ref + "/"
    refName = refDir + ref
    # now have reference file so we can proceed with the transcript mapping via bowtie2
    argString = "bowtie2 -x "+refName+" -1 "+forwardFN+" -2 "+reverseFN+" -S "+inputname+".sam"
    print(argString)
    subprocess.call(argString, shell=True)  #outputs a name.sam file
    return



def processSamFiles(inputname):
    cur_path = os.getcwd()
    samName = cur_path+"/"+inputname
    argString = "samtools view -bS "+inputname+".sam > "+samName+".bam"
    print(argString)
    subprocess.call(argString, shell=True)

    argString = "samtools sort "+samName+".bam -o "+samName+".sorted"
    print("argstring = "+argString)
    subprocess.call(argString, shell=True)

    argString = "samtools index "+samName+".sorted "+samName+".sorted.bai"
    print("argstring = " + argString)
    subprocess.call(argString, shell=True)
    return  #we have saved out the relevent name.bam, name.sorted and name.sorted.bai files

# we will not have the .gtf file so call cufflinks without -G option
def transcriptAbundance(inputname):
    argString = "cufflinks -o "+inputname+".cuff -u -p 8 "+inputname+".sorted"
    subprocess.call(argString, shell = True)
    os.remove(inputname+".sorted")  #remove name.sorted
    os.remove(inputname+".sorted.bai")
    os.remove(inputname+".bam")
    os.remove(inputname+".sam")
    return

def transcriptsForBlast(name, refFastq):
    # quick and dirty just to see.
    refBase = os.path.basename(refFastq)
    ref = os.path.splitext(refBase)[0]  # 'mydata/test.fasta' -> 'test'
    dir_path = os.path.dirname(os.path.realpath(__file__))  # directory of this file
    refPath = dir_path + "/data/Reference/" + ref + "/" + refBase   # eg refPath = data/Reference/Trinity/Trinity.fasta
    # used for dirty # refPath = 'Trinity.fasta' # dirty one
    track_df = pd.read_csv(name+'.cuff/genes.fpkm_tracking', sep='\t')
    names = track_df['locus']
    nlist = []
    for n in range(0,len(names)):
        i = names[n].find(':')
        nlist.append(names[n][:i])
    nameset = set(nlist)        #get unique.
    with open(refPath, 'r') as myRef:
        refData = myRef.read()
        refData= refData+'\n>'

    with open(name + '_for_blast.fa', 'w') as outfile:
        for trans_id in nameset:
            namepos = refData.find(trans_id)
            endpos = refData.find('>', namepos)
            outfile.write('>'+refData[namepos:endpos])
    pass

def blastContigs(test_name, database, old_name='old' ):
    if old_name == 'old':
        old_name = test_name
    db_path = database
    argString = "blastx -db "+db_path+" -query "+test_name+"_for_blast.fa -outfmt 10 -out "+test_name+"_blast.txt"
    print(argString)
    returncode = subprocess.call(argString, shell=True)
    if returncode != 0:
        return "Error in blastall"
    blast_df = pd.read_csv(""+test_name+"_blast.txt")
    blast_df.columns = ['qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue','bitscore']
    blastResult_df = blast_df[(blast_df['pident']>=70) & (blast_df['length'] > 100) & (blast_df['evalue'] <=0.001) ]
    blastResult_df = blastResult_df[['qaccver', 'saccver', 'pident', 'evalue', 'bitscore']]   #query accession.version, subject accession.version, Percentage of identical matches
    # need to allocate the transcripts (if allocated more than once to the phylotype with least error.
    transcripts = blastResult_df['qaccver']
    b_df = pd.DataFrame(columns=['qaccver', 'saccver', 'pident', 'evalue', 'bitscore'])
    transSet = set(transcripts)
    for t in transSet:
        temp_df = blastResult_df[(blastResult_df['qaccver'] == t)]
        # get one with smallest error value
        #print(t + ":")
        #print(temp_df)
        temp_df = temp_df.sort_values(by=['evalue'])
        b_df = b_df.append(temp_df.iloc[[0]])

    b_df.sort_values(by=['qaccver'])
    fdir = r"./results/" + old_name + "/"
    if not os.path.exists(fdir):
        os.makedirs(fdir)
    b_df.to_csv(fdir + test_name + '_transcript.csv')
    b_df.to_csv(test_name + '_transcript.csv')
    os.remove(test_name+"_for_blast.fa")  # remove name_for_blast_fa
    os.remove(test_name+"_blast.txt")  # remove name_for_blast_fa
    return b_df


def createMultiHTML(tdict,composite_df):
    labelList = composite_df.columns.tolist()
    htmlString = r"<html><title>T.vivax VAP (Transcriptomic Pathway(</title><body><div style='text-align:center'><h2><i>Trypanosoma vivax</i> Variant Antigen Profile</h2><h3>"
    htmlString += r"Sample name: "+tdict['name']
    htmlString += r"<br>Transcriptomic Analysis</h3></p>"
    htmlString += "<p style = 'margin-left:20%; margin-right:20%'>Legend: " \
                  "Variant Antigen Profile of a <i>Trypanosoma vivax</i> transcriptomes. " \
                  "Weighted Frequency reflects Phylotype abundance and is expressed as " \
                  "phylotype frequencies adjusted for the combined transcript abundance. " \
                  "Data was produced with VAPPER-Variant Antigen Profiler " \
                  "(Silva Pereira et al., 2019).</p> "
    htmlString += r"<style> table, th, tr, td {border: 1px solid black; border-collapse: collapse;}</style>"

    header = r"<table style='width:50%;margin-left:25%;text-align:center'><tr><th>Phylotype</th>"
    wLists = []

    for j in range(1,len(labelList)):
        wLists.append(composite_df[labelList[j]])
        header += r"<th>" + str(labelList[j]) + "</th>"

    htmlString += "</tr>\n" + header
    tabString = ""
    phyList = composite_df['Phylotype']



    for i in range(0, len(composite_df)):
        tabString += "<tr><td>" + str(phyList[i]) + "</td>"
        for j in range(0,len(labelList)-1):
            #print(j)
            f = format(wLists[j][i], '.4f')
            tabString += "<td>" + str(f) + "</td>"
        tabString += "</tr>\n"

    htmlString += tabString + "</table><br><br><br><br><br>"
    htmlString += r"<h3>Weighted Relative Frequencies of Detected Phylotypes.</h3>"
    imgString = r"<img src = '"+ tdict['name']+"_phylotypes.png' alt='Bar chart of phylotype variation' style='max-width:100%'><br><br>"
    htmlString += imgString

    with open(tdict['html_file'], "w") as htmlfile:
        htmlfile.write(htmlString)


def createHTML(tdict,sum_df):
    #assumes imgs are heatmap.png, dheatmap.png, vapPCA.png and already in htmlresource
    htmlString = r"<html><title>T.vivax VAP (Transcriptomic Pathway(</title><body><div style='text-align:center'><h2><i>Trypanosoma vivax</i> Variant Antigen Profile</h2><h3>"
    htmlString += r"Sample name: "+tdict['name']
    htmlString += r"<br>Transcriptomic Analysis</h3></p>"
    htmlString += "<p style = 'margin-left:20%; margin-right:20%'>Legend: " \
                  "Variant Antigen Profile of a <i>Trypanosoma vivax</i> transcriptome. " \
                  "Weighted Frequency reflects Phylotype abundance and is expressed as " \
                  "phylotype frequencies adjusted for the combined transcript abundance. " \
                  "Data was produced with VAPPER-Variant Antigen Profiler " \
                  "(Silva Pereira et al., 2019).</p> "
    htmlString += r"<style> table, th, tr, td {border: 1px solid black; border-collapse: collapse;}</style>"

    htmlString += r"<table style='width:50%;table-layout: auto; margin-left:25%;text-align:center'><tr><th>Phylotype</th>" \
                  r"<th>Combined FPKM</th><th>Weighted Frequency</th></tr>"
    tabString = ""
    # flush out table with correct values
    phySeries = sum_df['Phylotype']
    # sacSeries = sum_df['saccver']
    fSeries = sum_df['FPKM']
    total = fSeries.sum()
    # print("Total="+str(total))
    for i in range(0, len(sum_df)):
        # print(phySeries[i])
        f = format(fSeries[i], '.2f')
        w = format(fSeries[i]/total, '.2f')

        #w = format(weightList[i], '.4f')

        tabString += "<tr><td>" + str(phySeries[i]) + "</td><td>" + str(f) + "</td><td>"+str(w)+"</tr>"
    htmlString += tabString + "</table><br><br><br><br><br>"
    htmlString += r"<h3>Weighted Relative Frequencies of Detected Phylotypes.</h3>"
    imgString = r"<img src = '"+ tdict['name']+"_phylotypes.png' alt='Bar chart of phylotype variation' style='max-width:100%'><br><br>"
    htmlString += imgString

    with open(tdict['html_file'], "w") as htmlfile:
        htmlfile.write(htmlString)



def getPhyloNumber(sac):
    i = sac.find('_')
    return int(sac[1:i])

def combineFPMK(tdict):
    fpkm_df = pd.read_csv(tdict['name']+'.cuff/genes.fpkm_tracking', sep='\t')

    #fpkm_df = pd.read_csv('genes.fpkm_tracking',sep='\t')
    #print(fpkm_df.head())
    fpkm_df['locus'] = fpkm_df['locus'].apply(lambda names: names[:names.find(':')])
    #print(fpkm_df.head())
    reducedBlast_df = pd.read_csv(tdict['name']+'_transcript.csv')
    # reducedBlast_df = pd.read_csv('TrinityVT_transcript.csv')
    saccverSet = set(reducedBlast_df['saccver'])
    saccverList = list(saccverSet)
    saccverList.sort()
    # print(saccverList[:5])
    new_df = pd.DataFrame(columns=['qaccver','saccver','FPKM'])
    for sv in saccverList:
        #print(sv)
        temp_df = reducedBlast_df[reducedBlast_df['saccver'] == sv]
        qList = list(temp_df['qaccver'])
        for q in qList:
            f_df = fpkm_df[(fpkm_df['locus'] == q)]
            if len(f_df) > 1:
                print('WARNING MULTIPLE FPKM')
            new_fpkm=list(f_df['FPKM'])
            f = (new_fpkm[0])
            # print(f)
            new_df = new_df.append({'qaccver': q, 'saccver': sv, 'FPKM': f}, ignore_index=True)
    FPKMsum_df = new_df.groupby('saccver')['FPKM'].sum().reset_index()

    FPKMsum_df['Phylotype'] = FPKMsum_df.apply(lambda row: getPhyloNumber(row['saccver']), axis=1)
    FPKMsum_df = FPKMsum_df.sort_values(by=['Phylotype'])
    FPKMsum_df = FPKMsum_df.reset_index(drop=True)

    # print(FPKMsum_df)
    FPKMsum_df.to_csv('FPKM_sum.csv')
    FPKMsum2_df = FPKMsum_df.groupby('Phylotype')['FPKM'].sum().reset_index()
    FPKMsum2_df = FPKMsum2_df.sort_values(by=['Phylotype'])

    # print(FPKMsum2_df)
    FPKMsum2_df.to_csv('FPKM_sum2.csv') # in case more than one entry for a particular phylotype

    os.remove(tdict['name']+'_transcript.csv')

    return FPKMsum_df, FPKMsum2_df



def normalisef(f,max):
    return f/max

def getComposite_sum2(nameList,sum2_dfs):
    # lets get a composite sum2_df from all of the sum2_dfs
    phyList = []

    for i in range(0, len(sum2_dfs)):
        total = sum2_dfs[i]['FPKM'].sum()
        sum2_dfs[i]['w'] = sum2_dfs[i].apply(lambda row: normalisef(row['FPKM'], total), axis=1)
        pSeries = sum2_dfs[i]['Phylotype']
        for p in pSeries:
            phyList.append(p)  # get all the phylotypes in this one
    phyList = list(set(phyList))
    phyList.sort()
    composite_sum2_df = pd.DataFrame(phyList, columns=['Phylotype'])
    for i in range(0, len(sum2_dfs)):
        wList = []
        pindf = list(sum2_dfs[i]['Phylotype'])
        # print(pindf)
        for p in phyList:
            if p in pindf:
                df = sum2_dfs[i]
                w = df.loc[df['Phylotype'] == p, 'w'].iloc[0]
            else:
                w = 0
            wList.append(w)
        composite_sum2_df[nameList[i]] = wList
    #print(composite_sum2_df)
    #composite_sum2_df.to_csv('composite.csv')
    return composite_sum2_df


def doMultiBarChart(tdict, composite_df):       #array of multiple sum2_dfs
    labelList = composite_df.columns.tolist()
    sampnum = len(labelList)-1
# need to arrange bars
# number of phylotype = len(composite_df)
#number of bars = (len(labelist)-1) +1 for space
# ytick needs to ne

    cmap = plt.cm.get_cmap('tab10')
    palette = [cmap(i) for i in range(cmap.N)]
    title = "Legend: Variant Antigen Profile of a $\itTrypanosoma$ $\itvivax$ transcriptomes. " \
            "Phylotype abundance is expressed as phylotype frequencies adjusted " \
            "for combined transcript abundance. " \
            "Data was produced with VAPPER-Variant Antigen Profiler (Silva Pereira et al., 2019)."
    width = 0.6
    ind = np.arange(width*sampnum/2, len(composite_df)*width*(sampnum+1), width*(sampnum+1))
    print(ind)
    ysize = len(composite_df)*0.4

    fig, ax = plt.subplots(figsize=(10,ysize))


    for s in range(1, len(labelList)):
        ax.barh(ind, composite_df[labelList[s]], width, color=palette[s], label=labelList[s])
        ind = ind + width

    ax.set(yticks=np.arange(width*(sampnum+2)/2, len(composite_df)*width*(sampnum+1), width*(sampnum+1)), yticklabels=composite_df['Phylotype']) # , ylim=[(len(labelList)-1) * width - 1, len(composite_df)])
    ax.legend()


    ax.set_ylabel('Phylotype')
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel('Weighted Phylotype Frequency')

    # plt.text(-0.3, -0.15, title, va="top", wrap="True")
    #plt.tight_layout()

    plt.subplots_adjust(bottom=0.1, top=0.92, left=0.15, right=0.9)
    ax.set_title(title, x=0, wrap='True',ha='left',)

    plt.savefig(tdict['html_resource'] + tdict['name']+"_phylotypes.png")
    if tdict['pdf'] == 'PDF_Yes':
        plt.savefig(tdict['html_resource'] + tdict['name']+"phylotypes.pdf")
    plt.show()
    pass



def doBarChart(tdict, sum2_df):
    cmap = plt.cm.get_cmap('tab20')
    palette = [cmap(i) for i in range(cmap.N)]
    title = "Legend: Variant Antigen Profile of a $\itTrypanosoma$ $\itvivax$ transcriptome. " \
            "Phylotype abundance is expressed as phylotype frequencies adjusted " \
            "for combined transcript abundance. " \
            "Data was produced with VAPPER-Variant Antigen Profiler (Silva Pereira et al., 2019)."
       # get a list of phylotype, create equivalent of saccver, get a list of
    maxFPKM = sum2_df['FPKM'].max()
    total = sum2_df['FPKM'].sum()

    sum2_df['Normalised'] = sum2_df.apply(lambda row: normalisef(row['FPKM'], maxFPKM),axis=1)
    sum2_df['Weighted'] = sum2_df.apply(lambda row: normalisef(row['FPKM'], total),axis=1)
    pList = sum2_df['Phylotype']
    phList = []
    for p in pList:
        phList.append(str(p))

    fList = sum2_df['Weighted']
    ysize = len(phList)*0.3
    fig, ax = plt.subplots(figsize=(10,ysize))

    ax.barh(phList, fList, color=palette)
    ax.set_ylabel('Phylotype')
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel('Weighted Phylotype Frequency')

    # plt.text(-0.3, -0.15, title, va="top", wrap="True")
    #plt.tight_layout()
    plt.subplots_adjust(bottom=0.1, top=0.9, left=0.15, right=0.9)
    ax.set_title(title, x=0, wrap='True',ha='left',)

    plt.savefig(tdict['html_resource'] + tdict['name']+"_phylotypes.png")
    if tdict['pdf'] == 'PDF_Yes':
        plt.savefig(tdict['html_resource'] + tdict['name']+"phylotypes.pdf")
    # plt.show()
    pass


def transcriptomicProcess(tdict):
    uploadUserReferenceFastq(tdict['refFastq'])
    transcriptMapping(tdict['name'], tdict['refFastq'], tdict['forward'], tdict['reverse'])        #uses bowtie
    processSamFiles(tdict['name'])                              #uses samtools
    transcriptAbundance(tdict['name'])                          #uses cufflinks -> ?.cuff/*.*
    transcriptsForBlast(tdict['name'], tdict['refFastq'])       #creates name+4blast.fa
    blastContigs(tdict['name'],'data/vivax/Database/Phylotype_typeseqs.fas')
    sum_df, sum2_df = combineFPMK(tdict)
    shutil.rmtree(tdict['name'] + '.cuff')
    doBarChart(tdict, sum2_df)
    createHTML(tdict, sum_df)


if __name__ == "__main__":


    tdict = {}
    tdict['refFastq'] = "vivax_trans/Trinity.fasta"
    tdict['forward'] = 'vivax_trans/Tv1392_BF1.fq'
    tdict['reverse'] = 'vivax_trans/Tv1392_BF2.fq'
    tdict['name'] = 'TrinityVT'
    tdict['strain'] = 'TrinityVTstrain'
    tdict['vivax_trans_database'] = 'data/vivax/Database/Phylotype_typeseqs.fas'
    htmlfile = r"./results/" + tdict['name'] + "/" + tdict['name'] + ".html"
    htmldir = r"./results/" + tdict['name'] + "/"
    if not os.path.exists(htmldir):
        os.makedirs(htmldir)
    tdict['html_file'] = htmlfile
    tdict['html_resource'] = htmldir
    tdict['pdf']='PDF_No'

    sum_df, sum2_df = combineFPMK(tdict)
    doBarChart(tdict, sum2_df)
    createHTML(tdict, sum_df)
    exit()

    nameList = ['multi_One', 'multi_Two', 'multi_Three', 'multi_Four']

    sum2_df0 = pd.read_csv('FPKM_sum2.csv')
    sum2_df1 = pd.read_csv('FPKM_1_sum2.csv')
    sum2_df2 = pd.read_csv('FPKM_2_sum2.csv')
    sum2_df3 = pd.read_csv('FPKM_3_sum2.csv')

    sum2_dfs = [sum2_df0, sum2_df1, sum2_df2, sum2_df3]
    composite_df = getComposite_sum2(nameList,sum2_dfs)
    doMultiBarChart(tdict, composite_df)
    createMultiHTML(tdict,composite_df)
    exit()
    #sum_df, sum2_df = combineFPMK(tdict)
    #doBarChart(tdict, sum2_df)
    #createHTML(tdict, sum_df)

    # blastContigs(tdict['name'], tdict['vivax_trans_database'])
    exit()
    #
    # transcriptomicProcess(tdict)

    #transcriptsForBlast(tdict['name'], tdict['refFastq'])
    #transcriptAbundance(tdict['name'])                          #uses cufflinks -> ?.cuff/*.*
    # blastContigs(tdict['name'],tdict['vivax_trans_database'])



