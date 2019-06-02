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


import subprocess
import pandas as pd
import re
import os
import sys
import shutil
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

pList = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15']

def transcriptMapping(inputname, strain, forwardFN,reverseFN):
    #where is our Reference data -
    dir_path = os.path.dirname(os.path.realpath(__file__))
    refName = dir_path+"/data/Reference/Tc148" #default
    if strain == "Tc148":
        refName = dir_path+"/data/Reference/Tc148"
    if strain == "IL3000":
        refName = dir_path+"/data/Reference/IL3000"
    #now have reference file so we can proceed with the transcript mapping via bowtie2
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

def transcriptAbundance(inputname, strain):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    refName = dir_path + "/data/Reference/ORFAnnotation.gtf"  # default
    if strain == "Tc148":
        refName = dir_path + "/data/Reference/ORFAnnotation.gtf"    #still default
    if strain == "IL3000":
        refName = dir_path + "/data/Reference/IL3000.gtf"
    argString = "cufflinks -q -G "+refName+" -o "+inputname+".cuff -u -p 8 "+inputname+".sorted"
    subprocess.call(argString, shell = True)
    os.remove(inputname+".sorted")  #remove name.sorted
    os.remove(inputname+".sorted.bai")
    os.remove(inputname+".bam")
    return

def convertToFasta(inputName, strain):  #equivalent to Sara's awk script
    dir_path = os.path.dirname(os.path.realpath(__file__))
    refName = dir_path + "/data/Reference/ORFAnnotation.gtf"  # default
    if strain == "Tc148":
        refName = dir_path + "/data/Reference/148_prot.fasta"
    if strain == "IL3000":
        refName = dir_path + "/data/Reference/IL3000_prot.fasta"

    cuff_df = pd.read_csv(inputName+".cuff/genes.fpkm_tracking", sep='\t')
    cuff_df = cuff_df[(cuff_df['FPKM'] > 0)]
    cuff_df.to_csv("cuffTest.csv")
    gene_id_List = cuff_df['gene_id'].tolist()

    number = 0
    all = 0
    with open(inputName+"_6frame.fas", 'w') as outfile:
        ref = open(refName,'r')
        #ref = open(r"Reference/IL3000_prot.fasta",'r')
        n = 0
        line = ref.readline()
        while line:
            if line[0] == '>':
                all = all+1
                ln = line[1:]   #remove >
                ln = ln.rstrip()     #remove /n /r etc
                #print (ln)
                if ln in gene_id_List:
                    number = number+1
                    outfile.write(line)
                    line = ref.readline()
                    if line:
                        while line[0] != '>':
                            outfile.write(line)
                            line=ref.readline()
                            if not line:
                                break;
                else:
                    line = ref.readline()
            else:
                line =ref.readline()
    ref.close()
    return cuff_df

def HMMerMotifSearch(name, strain, cuff_df):
    motifs = ['1', '2a', '2b', '3', '4a', '4b', '4c', '5', '6', '7', '8a', '8b', '9a', '9b',
              '9c', '10a', '10b', '11a', '11b', '12', '13a', '13b', '13c', '13d', '14', '15a', '15b', '15c']
    dir_path = os.path.dirname(os.path.realpath(__file__))
    phylopath = dir_path + "/data/Motifs/Phylotype"
    lineCounts = []
    compoundList = []
    for m in motifs:
        argString = "hmmsearch "+phylopath + m + ".hmm " + name + "_6frame.fas > Phy" + m + ".out"
        print(argString)
        subprocess.call(argString, shell=True)
        hmmResult = open("Phy" + m + ".out", 'r')
        regex = r"Tc148[0-9]{1,8}"
        if strain == "Tc148":
            regex = r"Tc148[0-9]{1,8}"
        if strain == "IL3000":
            regex = r"TcIL3000_[0-9]{1,4}_[0-9]{1,5}"
        n = 0
        outList = []
        for line in hmmResult:
            ms = re.search(regex, line)
            if ms:
                outList.append(""+ms.group())
                n += 1
            if re.search(r"inclusion", line):
                print("inclusion threshold reached")
                break
        compoundList.append(outList)
        lineCounts.append(n)
        hmmResult.close()
        os.remove("Phy" + m + ".out")

    concatGroups = [1, 2, 1, 3, 1, 1, 1, 2, 3, 2, 2, 1, 4, 1, 3]
    countList = []
    weightList = []
    totalCount = 0
    totalWeigth = 0
    for c in concatGroups:
        a = []
        weight = []
        for n in range(0, c):
            a = a + compoundList.pop(0)
        t = set(a)
        countList.append(len(t))
        wa = 0
        for w in t:
            wt = cuff_df.loc[cuff_df['gene_id'] == w, 'FPKM'].iloc[0]
            #print(w)
            #print(wt)
            wa = wa+wt
        weightList.append(wa)
        totalWeigth+=wa
        totalCount += len(t)
    countList.append(totalCount)
    weightList.append(totalWeigth)
    os.remove(name + "_6frame.fas")
    shutil.rmtree(name+".cuff")
    return countList,weightList


def relativeFrequencyTableNoSave(countList):
    relFreqList = []
    c = float(countList[15])
    for i in range(0, 15):
        relFreqList.append(countList[i] / c)
    return relFreqList

def relativeFrequencyTable(countList, name, htmlresource):
    relFreqList = []
    c = float(countList[15])
    for i in range(0, 15):
        relFreqList.append(countList[i] / c)

    data = {'Phylotype': pList, 'Relative Frequency': relFreqList}
    relFreq_df = pd.DataFrame(data)
    j_fname = htmlresource+ "/" + name + "_t_relative_frequency.csv"
    relFreq_df.to_csv(j_fname)
    return relFreqList


def weightedFrequencyTableNoSave(countList):
    relFreqList = []
    c = float(countList[15])
    for i in range(0, 15):
        relFreqList.append(countList[i] / c)
    return relFreqList

def weightedFrequencyTable(countList, name, htmlresource):
    relFreqList = []
    c = float(countList[15])
    for i in range(0, 15):
        relFreqList.append(countList[i] / c)

    data = {'Phylotype': pList, 'Weighted Frequency': relFreqList}
    relFreq_df = pd.DataFrame(data)
    j_fname = htmlresource+ "/" + name + "_t_weighted_frequency.csv"
    relFreq_df.to_csv(j_fname)
    return relFreqList



def createStackedBar(name,freqList,strain,pdf,html_resource):

    VAP_148 = [0.072, 0.032, 0.032, 0.004, 0.007,
               0.005, 0.202, 0.004, 0.006, 0.014,
               0.130, 0.133, 0.054, 0.039, 0.265]

    VAP_IL3000 = [0.073, 0.040, 0.049, 0.018, 0.060,
                  0.055, 0.054, 0.025, 0.012, 0.060,
                  0.142, 0.100, 0.061, 0.078, 0.172]
    cmap = plt.cm.get_cmap('tab20')
    palette = [cmap(i) for i in range(cmap.N)]

    if strain == "Tc148":
        VAPtable = VAP_148
        VAPname='Tc148\nGenome VAP'
    if strain == "IL3000":
        VAPtable = VAP_IL3000
        VAPname= 'IL3000\nGenome VAP'
    width = 0.35  # the width of the bars: can also be len(x) sequence
    plots = []
    fpos = 0
    vpos = 0
    for p in range(0, 15):
        tp = plt.bar(0, freqList[p], width, color= palette[p], bottom = fpos)
        fpos +=freqList[p]

        tp = plt.bar(1, VAPtable[p], width, color= palette[p], bottom = vpos)
        vpos +=VAPtable[p]

        plots.append(tp)
    plt.xticks([0,1],[name,VAPname])
    plt.legend(plots[::-1],['p15','p14','p13','p12','p11','p10','p9','p8','p7','p6','p5','p4','p3','p2','p1'])
    title = "Figure Legend: The transcriptomic Variant Antigen Profile of $\itTrypanosoma$ $\itcongolense$ estimated as phylotype " \
            "proportion adjusted for transcript abundance and the reference genomic Variant Antigen Profile. " \
            "\nData was produced with VAPPER-Variant Antigen Profiler (Silva Pereira et al., 2019)."
    plt.text(-0.3, -0.15, title, va="top", wrap="True")
    plt.tight_layout(pad=1.5)
    plt.subplots_adjust(bottom = 0.3,top=0.99,left=0.125,right=0.9,hspace=0.2,wspace=0.2)

    plt.savefig(html_resource + "/stackedbar.png")
    if pdf == 'PDF_Yes':
        plt.savefig(html_resource + "/stackedbar.pdf")



def createHTML(name,htmlfn,htmlresource,freqList,weightList):
    #assumes imgs are heatmap.png, dheatmap.png, vapPCA.png and already in htmlresource
    htmlString = r"<html><title>T.congolense VAP</title><body><div style='text-align:center'><h2><i>Trypanosoma congolense</i> Variant Antigen Profile</h2><h3>"
    htmlString += name
    htmlString += r"<br>Transcriptomic Analysis</h3></p>"
    htmlString += "<p style = 'margin-left:20%; margin-right:20%'>Table Legend: Variant Antigen Profiles of a transcriptome of <i>Trypanosoma congolense</i> estimated as phylotype proportion. " \
                  "Weighted frequency refers to the phylotype proportion based transcript abundance. " \
                  "Data was produced with VAPPER-Variant Antigen Profiler (Silva Pereira et al., 2019).</p> "
    htmlString += r"<style> table, th, tr, td {border: 1px solid black; border-collapse: collapse;}</style>"

    htmlString += r"<table style='width:50%;margin-left:25%;text-align:center'><tr><th>Phylotype</th><th>Relative Frequency</th><th>Weighted Frequency</th></tr>"
    tabString = ""
    # flush out table with correct values
    for i in range(0, 15):
        f = format(freqList[i], '.4f')
        w = format(weightList[i], '.4f')
        tabString += "<tr><td>phy" + str(i + 1) + "</td><td>"  + f + "</td><td>" + w + "</td></tr>"
    htmlString += tabString + "</table><br><br><br><br><br>"
    htmlString += r"<p> <h3>Stacked Bar chart of Phylotype Frequency</h3> The 'weighted' relative frequency of each phylotype alongside the VAP of selected strain.</p>"
    imgString = r"<img src = 'stackedbar.png' alt='Stacked bar chart of phylotype variation' style='max-width:100%'><br><br>"
    htmlString += imgString

    with open(htmlfn, "w") as htmlfile:
        htmlfile.write(htmlString)


def transcriptomicProcess(dict):
    transcriptMapping(dict['name'], dict['strain'], dict['forward'], dict['reverse'])        #uses bowtie
    processSamFiles(dict['name'])                              #uses samtools
    transcriptAbundance(dict['name'],dict['strain'])                          #uses cufflinks -> ?.cuff/*.*
    cuff_df = convertToFasta(dict['name'],dict['strain'])
    countList, weightList = HMMerMotifSearch(dict['name'],dict['strain'], cuff_df)
    relFreqList = relativeFrequencyTable(countList,dict['name'],dict['html_resource'])
    relWeightList = weightedFrequencyTable(weightList,dict['name'],dict['html_resource'])
    createStackedBar(dict['name'],relWeightList, dict['strain'],dict['pdf'],dict['html_resource'])
    print("Placing results in " + dict['html_resource'])
    createHTML(dict['name'],dict['html_file'],dict['html_resource'], relFreqList, relWeightList)

if __name__ == "__main__":
    sys.exit()

