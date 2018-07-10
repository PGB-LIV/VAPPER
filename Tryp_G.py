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
import re
import os
import sys
import shutil
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.mlab import PCA
import seaborn as sns

# some globals for convenience
pList = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15']

#align and assemble NSG paired reads fastq1name and fastq2name
def assembleWithVelvet(name, kmers, inslen, covcut, fastq1name,fastq2name):
    argString = "velveth " + name + "_k"+ kmers+" "+ kmers + " -shortPaired -fastq " + fastq1name+" "+fastq2name
    print(argString)
    returncode = subprocess.call(argString, shell=True)
    if returncode != 0:
        return "Error in velveth"
    argString = "velvetg " + name + "_k"+kmers+" -exp_cov auto -ins_length "+inslen+" -cov_cutoff "+covcut+" -clean yes -ins_length_sd 50 -min_pair_count 20"
    print(argString)
    returncode = subprocess.call(argString, shell = True)
    if returncode != 0:
        return "Error in velvetg"
    shutil.copyfile(name + "_k"+kmers+"//contigs.fa",name + ".fa")  #copy contigs file from _k65 to working directory
    shutil.rmtree(name+"_k"+kmers)   #remove temporary directory eg Test1_1_k65
    return

#translate the contig file using EMBOSS.transeq (translate in all six forward and reverse frames)
def contigTranslation(name):
    argString = "transeq " + name + ".fa " + name + "_6frame.fas -frame=6 "
    print(argString)
    returncode = subprocess.call(argString, shell=True)
    #os.remove(name+".fa")   #remove unwanted fasta file
    return returncode

#use hmmsearch to search each of the profiles against name_6frams.fas (fasta sequence)
#then concatenate the counts for each phyloypte
def HMMerMotifSearch(name):
    motifs = ['1', '2a', '2b', '3', '4a', '4b', '4c', '5', '6', '7', '8a', '8b', '9a', '9b',
              '9c', '10a', '10b', '11a', '11b', '12', '13a', '13b', '13c', '13d', '14', '15a', '15b', '15c']
    lineCounts = []
    compoundList = []
    dir_path = os.path.dirname(os.path.realpath(__file__))
    phylopath = dir_path + "/data/Motifs/Phylotype"
    for m in motifs:
        argString = "hmmsearch " + phylopath + m + ".hmm " + name + "_6frame.fas > Phy" + m + ".out"
        print(argString)
        subprocess.call(argString, shell=True)

        hmmResult = open("Phy" + m + ".out", 'r')
        n = 0
        outList = []
        for l in range(0,14):
            hmmResult.readline()        #miss out the first 14 lines. data we want starts on line 15

        for line in hmmResult:
            if re.search(r"inclusion", line):   #inclusion reached
                break
            if len(line) <= 1:  #end of data
                break
            nod = line[60:-1]
            outList.append("" + nod + "\n")
            n += 1

        compoundList.append(outList)
        lineCounts.append(n)
        hmmResult.close()
        os.remove("Phy" + m + ".out")

    #concatenate
    concatGroups = [1, 2, 1, 3, 1, 1, 1, 2, 3, 2, 2, 1, 4, 1, 3]
    countList = []
    totalCount = 0
    for c in concatGroups:
        a = []
        for n in range(0, c):
            a = a + compoundList.pop(0)
        t = set(a)
        countList.append(len(t))
        totalCount += len(t)
    countList.append(totalCount)
    os.remove(name + "_6frame.fas")
    return countList

#calculates the relative frequency for multiple samples (hence no save0
def relativeFrequencyTableNoSave(countList):
    relFreqList = []
    c = float(countList[15])
    if c == 0:
        return [0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0]
    for i in range(0, 15):
        relFreqList.append(countList[i] / c)
    return relFreqList  # 0-14 = p1-p15 counts [15] = total counts

#calculates the relative frequency for single sample and svaes it as name_relative_frequency.csv
def relativeFrequencyTable(countList, name, htmlresource):
    relFreqList = []
    c = float(countList[15])
    if c == 0:
        return [0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0]
    for i in range(0, 15):
        relFreqList.append(countList[i] / c)
    data = {'Phylotype': pList, 'Relative Frequency': relFreqList}
    relFreq_df = pd.DataFrame(data)
    j_fname = htmlresource+"/" + name + "_relative_frequency.csv"
    relFreq_df.to_csv(j_fname)
    return relFreqList  # 0-14 = p1-p15 counts [15] = total counts

#calculates the deviaton from the mean for multiple samples (hence no save)
def getDeviationFromMeanNoSave(frequencyList):
    devList = []
    dir_path = os.path.dirname(os.path.realpath(__file__))
    j_fname = dir_path + "/data/congodata.csv"
    #j_fname = r"data/congodata.csv"
    congo_df = pd.read_csv(j_fname)  # we get the means from congo_df
    for p in range(0, 15):
        m = congo_df[pList[p]].mean()
        dev = -(m - frequencyList[p])
        devList.append(dev)
    return devList

#calculates the deviaton from the mean for individual samples and saves as name_deviation_from_mean.csv
def getDeviationFromMean(frequencyList, name, htmlresource):
    devList = []
    dir_path = os.path.dirname(os.path.realpath(__file__))
    j_fname = dir_path + "/data/congodata.csv"
    congo_df = pd.read_csv(j_fname)  # we get the means from congo_df
    for p in range(0, 15):
        m = congo_df[pList[p]].mean()
        dev = -(m - frequencyList[p])
        devList.append(dev)
    data = {'Phylotype': pList, 'Deviation from Mean': devList}
    dev_df = pd.DataFrame(data)
    j_fname = htmlresource+"/" + name + "_deviation_from_mean.csv"
    dev_df.to_csv(j_fname)
    return devList

#create the relative frequency heatmap
def relativeFrequencyHeatMap(name, freqList, pdf, htmlresource):
    localFreqList = freqList[:]
    localFreqList.insert(0, name)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    j_fname = dir_path+"/data/congodata.csv"
    congo_df = pd.read_csv(j_fname)
    congo_df.drop('Colour', axis=1, inplace=True)
    congo_df.loc[congo_df.index.max() + 1] = localFreqList
    congo_df.set_index('Strain', inplace=True)

    cg = sns.clustermap(congo_df, method='ward', cmap = "RdBu_r", col_cluster=False, yticklabels = congo_df.index.values)
    plt.setp(cg.ax_heatmap.yaxis.get_ticklabels(), rotation=0, fontsize=8)  # get y labels printed horizontally
    ax=cg.ax_heatmap
    title = "Variant Antigen Profiles of $\itTrypanosoma$ $\itcongolense$ estimated as the phylotype proportion across the\nsample cohort. "
    title += "Dendrogram reflects the relationships amongst the VSG repertoires of each strain. "
    title += "Strains\nwere isolated from multiple African countries as described in Silva Pereira et al. (2018)."
    title += "\nData was produced with the 'Variant Antigen Profiler' (Silva Pereira and Jackson, 2018)."
    ax.text(-0.15,-0.05, title,va = "top",wrap = "True", transform = ax.transAxes )

    plt.savefig(htmlresource+"/heatmap.png",bbox_inches='tight')
    if pdf == 'PDF_Yes':
        plt.savefig(htmlresource+"/heatmap.pdf", bbox_inches='tight')
    return
#create the deviation from the mean frequency heatmap
def deviationFromMeanHeatMap(name,devList, pdf, htmlresource):
    localDevList = devList[:]
    localDevList.insert(0, name)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    j_fname = dir_path+ "/data/congodata_deviationfromthemean.csv"
    congo_df = pd.read_csv(j_fname)
    congo_df.drop('Colour', axis=1, inplace=True)
    congo_df.loc[congo_df.index.max() + 1] = localDevList
    congo_df.set_index('Strain', inplace=True)
    cg = sns.clustermap(congo_df, method='ward',cmap = "RdBu_r", col_cluster=False, yticklabels = congo_df.index.values)
    plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize=8)  # get y labels printed horizontally
    ax = cg.ax_heatmap
    title = "Variant Antigen Profiles of $\itTrypanosoma$ $\itcongolense$ expressed as the deviation from the mean phylotypes "
    title +="\nproportions of the sample cohort. Dendrogram reflects the relationships amongst the VSG repertoires of "
    title +="each \nstrain. Strains were isolated from multiple African countries as described in Silva Pereira et al. (2018)."
    title +="\nData was produced with the 'Variant Antigen Profiler' (Silva Pereira and Jackson, 2018)."
    ax.text(-0.2, -0.05, title, va="top", transform=ax.transAxes, wrap="True")
    plt.savefig(htmlresource+"/dheatmap.png",bbox_inches='tight')
    if pdf == 'PDF_Yes':
        plt.savefig(htmlresource+"/dheatmap.pdf", bbox_inches='tight')
    return

#PCA analysis for sample freqList
def prepAndPlotPCA(name,freqList,pdf,htmlresource):
    localFreqList = freqList[:]
    localFreqList.insert(0, name)
    localFreqList.append(name)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    j_fname = dir_path + "/data/congodata.csv"
    congo_df = pd.read_csv(j_fname)
    congo_df.loc[congo_df.index.max() + 1] = localFreqList
    plotPCA(congo_df,pdf,htmlresource)

#do a qick PCA analysis. NB. PCA is now deprecated in 2.2 - may need updating soon
def plotPCA(congo_df,pdf,htmlresource):
    myColours = congo_df['Colour']
    myCountries = congo_df.drop_duplicates('Colour')['Colour'].tolist()
    congo_df.drop('Colour', axis=1, inplace=True)
    congo_df.set_index('Strain', inplace=True)
    dataArray = congo_df.as_matrix()
    pcaResult = PCA(dataArray)
    compoundList = []
    for i in myCountries:
        compoundList.append([])
    i = 0
    for item in pcaResult.Y:
        col = myCountries.index(myColours[i])
        compoundList[col].append(-item[0])
        compoundList[col].append(item[1])
        i += 1

    colormap = plt.cm.tab20  # nipy_spectral, Set1,Paired
    cols = [colormap(i) for i in np.linspace(0, 1, 20)]
    fig, ax = plt.subplots(figsize=(9, 6))
    i = 0
    for d in myCountries:
        a = compoundList[i]
        b = a[::2]
        c = a[1::2]
        ax.scatter(b, c, color=cols[i], label=myCountries[i])
        i = i + 1
    leg = ax.legend( bbox_to_anchor=(1.02,1.02), loc = "upper left")        #move legend out of plot
    title = "Principal Component Analysis of the Variant Antigen Profiles of $\itTrypanosoma$ $\itcongolense$. " \
            "The plot reflects the\nrelationships amongst the VSG repertoires of each strain. Strains are color-coded " \
            "by location of collection according\nto key. Strains were isolated from multiple African countries as described in Silva Pereira et al. (2018)."
    title +="\nData was produced with the 'Variant Antigen Profiler' (Silva Pereira and Jackson, 2018)."
    tx = ax.text(-0.1, -0.07, title, va="top", transform=ax.transAxes, wrap="True")
    fig.subplots_adjust(bottom = 0.3)
    fig.savefig(htmlresource+"/vapPCA.png", bbox_extra_artists=(leg,tx), bbox_inches='tight')
    if pdf == 'PDF_Yes':
        fig.savefig(htmlresource+"/vapPCA.pdf",bbox_extra_artists=(leg,tx), bbox_inches='tight')
    return

def createHTML(name,htmlfn,freqList,devList):
    #assumes imgs are heatmap.png, dheatmap.png, vapPCA.png and already in htmlresource
    #creates an html file in results/name/name.html - displays images and tables
    htmlString = r"<html><title>T.congolense VAP</title><body><div style='text-align:center'><h2><i>Trypanosoma congolense</i> Variant Antigen Profile</h2><h3>"
    htmlString += name
    htmlString += r"<br/>Genomic Analysis</h3>"
    htmlString += "<p style = 'margin-left:23%; margin-right:23%'>Table Legend: Variant Antigen Profiles of <i>Trypanosoma congolense</i> estimated as the phylotype proportion and as the deviation from the mean across the sample cohort.<br>" \
                  "Data was produced with the 'Variant Antigen Profiler' (Silva Pereira and Jackson, 2018).</p>"
    htmlString += r"<style> table, th, tr, td {border: 1px solid black; border-collapse: collapse;}</style>"

    htmlString += r"<table style='width:50%;margin-left:25%;text-align:center'><tr><th>Phylotype</th><th>Relative Frequency</th><th>Deviation from Mean</th></tr>"
    tabString = ""
    # flush out table with correct values
    for i in range(0, 15):
        f= format(freqList[i],'.4f')
        d= format(devList[i],'.4f')
        tabString += "<tr><td>phy" + str(i + 1) + "</td><td>"  + f + "</td><td>" + d + "</td></tr>"
    htmlString += tabString + "</table><br><br><br><br><br>"

    htmlString += r"<h3>The Variation Heat Map and Dendrogram</h3><p>The absolute phylotype variation in the sample compared to model dataset.</p>"
    imgString = r"<img src = 'heatmap.png' alt='Variation Heatmap' style='max-width:100%'><br><br>"
    htmlString += imgString

    htmlString += r"<br><br><br><br><h3>The Deviation Heat Map and Dendrogram</h3><p>The phylotype variation expressed as the deviation from your sample mean compared to the model dataset</p>"
    imgString = r"<img src = 'dheatmap.png' alt='Deviation Heatmap' style='max-width:100%'><br><br>"
    htmlString += imgString

    htmlString += r"<br><br><br><br><h3>The Variation PCA plot</h3><p>PCA analysis corresponding to absolute variation. Colour coded according to location</p>"
    imgString = r"<img src = 'vapPCA.png' alt='PCA Analysis' style='max-width:100%'><br><br>"
    htmlString += imgString + r"</div></body></html>"

    with open(htmlfn, "w") as htmlfile:
        htmlfile.write(htmlString)


def assemble(dict):
    #assembly from NGS paired read files
    assembleWithVelvet(dict['name'],dict['kmers'], dict['inslen'],dict['covcut'], dict['forward'],dict['reverse'])
    contigTranslation(dict['name'])
    myCountList = HMMerMotifSearch(dict['name'])
    myFreqList = relativeFrequencyTable(myCountList, dict['name'],dict['html_resource'])  # saves out inputname_relative_frequncy.csv
    myDevList = getDeviationFromMean(myFreqList, dict['name'], dict['html_resource'])  # saves out inputname_deviation_from_mean.csv
    relativeFrequencyHeatMap(dict['name'], myFreqList,dict['pdf'], dict['html_resource'])
    deviationFromMeanHeatMap(dict['name'], myDevList,dict['pdf'], dict['html_resource'])
    prepAndPlotPCA(dict['name'], myFreqList,dict['pdf'], dict['html_resource'])
    createHTML(dict['name'], dict['html_file'], myFreqList, myDevList)  # assumes imgs are heatmap.png, dheatmap.png, vapPCA.png and already in htmlresource

def contigs(dict):
    #uses pre-assembled contigs
    shutil.copyfile(dict['contigs'],dict['name']+".fa")
    contigTranslation(dict['name'])
    myCountList = HMMerMotifSearch(dict['name'])
    myFreqList = relativeFrequencyTable(myCountList, dict['name'],dict['html_resource'])
    myDevList = getDeviationFromMean(myFreqList, dict['name'],dict['html_resource'])  # saves out inputname_deviation_from_mean.csv
    relativeFrequencyHeatMap(dict['name'], myFreqList, dict['pdf'], dict['html_resource'])
    deviationFromMeanHeatMap(dict['name'], myDevList, dict['pdf'], dict['html_resource'])
    prepAndPlotPCA(dict['name'], myFreqList,dict['pdf'],dict['html_resource'])
    createHTML(dict['name'],dict['html_file'], myFreqList,myDevList)  # assumes imgs are heatmap.png, dheatmap.png, vapPCA.png and already in htmlresource


if __name__ == "__main__":
    print("Error: Must be called from Vap.py")
    sys.exit()
