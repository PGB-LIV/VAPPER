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

quietString = ""	#" >>"+os.path.dirname(os.path.realpath(__file__))+"/log/Vap_log.txt 2>&1"

def assembleWithVelvet(name, kmers, inslen, covcut, fastq1name,fastq2name):
    #argString = "velveth " + name + "_k65 65 -shortPaired -fastq " + name + "_R1.fastq " + name + "_R2.fastq"
    argString = "velveth " + name + "_k"+ kmers+" "+ kmers + " -shortPaired -fastq " + fastq1name+" "+fastq2name+quietString
    print(argString)
    returncode = subprocess.call(argString, shell=True)
    if returncode != 0:
        return "Error in velveth"
    argString = "velvetg " + name + "_k"+kmers+" -exp_cov auto -ins_length "+inslen+" -cov_cutoff "+covcut+" -clean yes -ins_length_sd 50 -min_pair_count 20"+quietString
    #argString = "velvetg " + name + "_k65 -exp_cov auto -ins_length 400 -cov_cutoff 5 -clean yes -ins_length_sd 50 -min_pair_count 20"+quietString
    print(argString)
    returncode = subprocess.call(argString, shell = True)
    if returncode != 0:
        return "Error in velvetg"
    shutil.copyfile(name + "_k"+kmers+"//contigs.fa",name + ".fa")  # my $namechange = "mv ".$input."_k65/contigs.fa ".$input.".fa";
    return "ok"

def contigTranslation(name):
    argString = "transeq " + name + ".fa " + name + "_6frame.fas -frame=6 " #+quietString
    print(argString)
    returncode = subprocess.call(argString, shell=True)
    #subprocess.call('ls -l *.fa', shell = True)
    #sys.exit(1)
    #if returncode != 0:
    #    return "Error in Transeq"
    #return 'ok'


def HMMerMotifSearch(name):
    motifs = ['1', '2a', '2b', '3', '4a', '4b', '4c', '5', '6', '7', '8a', '8b', '9a', '9b',
              '9c', '10a', '10b', '11a', '11b', '12', '13a', '13b', '13c', '13d', '14', '15a', '15b', '15c']
    lineCounts = []
    compoundList = []
    dir_path = os.path.dirname(os.path.realpath(__file__))
    phylopath = dir_path + "/data/Motifs/Phylotype"
    for m in motifs:
        argString = "hmmsearch " + phylopath + m + ".hmm " + name + "_6frame.fas > Phy" + m + ".out"  # +quietString
        # argString = "hmmsearch "+phylopath + m + ".hmm " + dir_path+"/data/Test_6frame.fas > Phy" + m + ".out"
        #print(argString)
        subprocess.call(argString, shell=True)

        hmmResult = open("Phy" + m + ".out", 'r')
        #tempout = open(dir_path + "/data/" + "Phy" + m + ".txt", 'w')
        #regex = r"NODE_[0-9]{1,7}_length_[0-9]{1,7}_cov_[0-9]{1,10}.[0-9]{1,7}_[0-9]{1,2}"
        n = 0
        outList = []
        for l in range(0,14):
            hmmResult.readline()        #hacky? miss out the first 14 lines. data we want starts on line 15


        for line in hmmResult:
            if re.search(r"inclusion", line):
                #print("inclusion threshold reached")
                break
            if len(line) <= 1:
                #print("end of data")
                break
            nod = line[60:-1]
            #print(m)
            #tempout.write(m.group() + "\n")
            outList.append("" + nod + "\n")
            n += 1
        compoundList.append(outList)
        lineCounts.append(n)
        hmmResult.close()
        os.remove("Phy" + m + ".out")



    print(lineCounts)
    motifGroups = [['1'], ['2a', '2b'], ['3'], ['4a', '4b', '4c'], ['5'], ['6'], ['7'], ['8a', '8b'], ['9a', '9b',
                                                                                                       '9c'],
                   ['10a', '10b'], ['11a', '11b'], ['12'], ['13a', '13b', '13c', '13d'], ['14'], ['15a', '15b', '15c']]
    concatGroups = [1, 2, 1, 3, 1, 1, 1, 2, 3, 2, 2, 1, 4, 1, 3]
    countList = []
    countIndex = 0
    totalCount = 0

    for c in concatGroups:
        a = []
        for n in range(0, c):
            a = a + compoundList.pop(0)
        t = set(a)
        countList.append(len(t))
        totalCount += len(t)
    countList.append(totalCount)
    #print(countList)
    #print("--------")
    os.remove(name + "_6frame.fas")
    return countList

"""
def HMMerMotifSearch(name):
    motifs = ['1', '2a', '2b', '3', '4a', '4b', '4c', '5', '6', '7', '8a', '8b', '9a', '9b',
              '9c', '10a', '10b', '11a', '11b', '12', '13a', '13b', '13c', '13d', '14', '15a', '15b', '15c']
    lineCounts = []
    compoundList = []
    dir_path = os.path.dirname(os.path.realpath(__file__))
    phylopath = dir_path+"/data/Motifs/Phylotype"
    for m in motifs:
        argString = "hmmsearch "+phylopath + m + ".hmm " + name + "_6frame.fas > Phy" + m + ".out"  #+quietString
        #argString = "hmmsearch "+phylopath + m + ".hmm " + dir_path+"/data/Test_6frame.fas > Phy" + m + ".out"
        print(argString)
        subprocess.call(argString, shell=True)

        hmmResult = open("Phy" + m + ".out", 'r')
        tempout = open(dir_path+"/data/"+"Phy" + m + ".txt", 'w')
        regex = r"NODE_[0-9]{1,7}_length_[0-9]{1,7}_cov_[0-9]{1,10}.[0-9]{1,7}_[0-9]{1,2}"
        n = 0
        outList = []
        for line in hmmResult:
            m = re.search(regex, line)
            if m:
                tempout.write(m.group() + "\n")
                outList.append(""+m.group()+"\n")
                n += 1
            if re.search(r"inclusion", line):
                print("inclusion threshold reached")
                break
        compoundList.append(outList)
        lineCounts.append(n)
        hmmResult.close()
        #tempout.close()
    print(lineCounts)
    motifGroups = [['1'], ['2a', '2b'], ['3'], ['4a', '4b', '4c'], ['5'], ['6'], ['7'], ['8a', '8b'], ['9a', '9b',
                                                                                                   '9c'],
               ['10a', '10b'], ['11a', '11b'], ['12'], ['13a', '13b', '13c', '13d'], ['14'], ['15a', '15b', '15c']]
    concatGroups = [1, 2, 1, 3, 1, 1, 1, 2, 3, 2, 2, 1, 4, 1, 3]
    countList = []
    countIndex = 0
    totalCount = 0

    for c in concatGroups:
        a = []
        for n in range(0, c):
            a = a + compoundList.pop(0)
        t = set(a)
        countList.append(len(t))
        totalCount += len(t)
    countList.append(totalCount)
    print(countList)
    print("--------")
    return countList
"""



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




def getDeviationFromMean(frequencyList, name, htmlresource):
    devList = []
    dir_path = os.path.dirname(os.path.realpath(__file__))
    j_fname = dir_path + "/data/congodata.csv"
    #j_fname = r"data/congodata.csv"
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


def relativeFrequencyHeatMap(name, freqList, pdf, htmlresource):
    localFreqList = freqList[:]
    localFreqList.insert(0, name)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    j_fname = dir_path+"/data/congodata.csv"
    #print(dir_path)
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

    #title = "Variant Antigen Profiles of Trypanosoma congolense estimated as the phylotype proportion across the sample cohort. Dendrogram reflects the relationships amongst the VSG repertoires of each strain. Strains were isolated from multiple African countries as described in Silva Pereira et al. (2018). Data was produced with the 'Variant Antigen Profiler' (Silva Pereira and Jackson, 2018)."
    #ax.set_title(title, ha = "center", va = "bottom",wrap = "True")
    #title = "Where is this!"
    ax.text(-0.15,-0.05, title,va = "top",wrap = "True", transform = ax.transAxes )




    # cg.dendrogram_col.linkage  # linkage matrix for columns
    # cg.dendrogram_row.linkage  # linkage matrix for rows
    #plt.savefig(r"results/" + name + "_heatmap.png")
    plt.savefig(htmlresource+"/heatmap.png",bbox_inches='tight')
    if pdf == 'PDF_Yes':
        plt.savefig(htmlresource+"/heatmap.pdf", bbox_inches='tight')
        #shutil.copyfile("heatmap.pdf",heatmapfn)  #
    #plt.show()

def deviationFromMeanHeatMap(name,devList, pdf, htmlresource):
    localDevList = devList[:]
    localDevList.insert(0, name)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    j_fname = dir_path+ "/data/congodata_deviationfromthemean.csv"
    #j_fname = r"data/congodata_deviationfromthemean.csv"
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
    #ax.set_title(title,ha = "center", va = "bottom",wrap = "True")
    ax.text(-0.2, -0.05, title, va="top", transform=ax.transAxes, wrap="True")
    plt.savefig(htmlresource+"/dheatmap.png",bbox_inches='tight')
    if pdf == 'PDF_Yes':
        plt.savefig(htmlresource+"/dheatmap.pdf", bbox_inches='tight')
        #shutil.copyfile("dheatmap.pdf",dhmapfn)
    #plt.show()


def plotPCA(name, freqList, pdf, htmlresource):
    localFreqList = freqList[:]
    localFreqList.insert(0, name)
    localFreqList.append(name)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    j_fname = dir_path + "/data/congodata.csv"
    #j_fname = r"data/congodata.csv"
    congo_df = pd.read_csv(j_fname)
    congo_df.loc[congo_df.index.max() + 1] = localFreqList
    #    print(congo_df.tail(2))
    myColours = congo_df['Colour']
    myCountries = congo_df.drop_duplicates('Colour')['Colour'].tolist()
    #    print(myCountries)
    congo_df.drop('Colour', axis=1, inplace=True)
    congo_df.set_index('Strain', inplace=True)
    dataArray = congo_df.as_matrix()
    pcaResult = PCA(dataArray)
    # pcaResult.center(0)
    # can't seem to find a simple way of prooducing a decent legend.
    # going to seperate items in to different countires.
    compoundList = []
    for i in myCountries:
        compoundList.append([])

    i = 0
    for item in pcaResult.Y:
        col = myCountries.index(myColours[i])
        compoundList[col].append(-item[0])
        compoundList[col].append(item[1])
        i = i + 1
    cols = ['r', 'g', 'b', 'c', 'm', 'y', 'grey', 'k']

    fig, ax = plt.subplots(figsize=(9, 6))
    #plt.figure(num=1,figsize=(12, 6))
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
    #plt.title(title, ha = "center", va = "bottom",wrap = "True")
    tx = ax.text(-0.1, -0.07, title, va="top", transform=ax.transAxes, wrap="True")
    #fig.add_axes([0,0.05,1.05,1.05])
    #fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.subplots_adjust(bottom = 0.3)

    fig.savefig(htmlresource+"/vapPCA.png", bbox_extra_artists=(leg,tx), bbox_inches='tight')
    #fig.savefig(htmlresource+"/vapPCA.png", bbox_extra_artists=(leg,))
    if pdf == 'PDF_Yes':
        fig.savefig(htmlresource+"/vapPCA.pdf",bbox_extra_artists=(leg,tx), bbox_inches='tight')
        #shutil.copyfile("vapPCA.pdf",PCAfn)  # my $namechange = "mv ".$input."_k65/contigs.fa ".$input.".fa";
    #plt.show()

def createHTML(name,htmlfn,freqList,devList):
    #assumes imgs are heatmap.png, dheatmap.png, vapPCA.png and already in htmlresource
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
        #tabString += "<tr><td>phy" + str(i + 1) + "</td><td>"  + str(freqList[i]) + "</td><td>" + str(devList[i]) + "</td></tr>"
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


def assemble(args,dict):
    #argdict = {'name': 2, 'pdfexport': 3, 'kmers': 4, 'inslen': 5, 'covcut': 6, 'forward': 7, 'reverse': 8, 'html_file': 9,'html_resource': 10}
    assembleWithVelvet(args[dict['name']],args[dict['kmers']], args[dict['inslen']],args[dict['covcut']], args[dict['forward']],args[dict['reverse']])
    contigTranslation(args[dict['name']])
    myCountList = HMMerMotifSearch(args[dict['name']])
    myFreqList = relativeFrequencyTable(myCountList, args[dict['name']],args[dict['html_resource']])  # saves out inputname_relative_frequncy.csv
    # myFreqList = [0.111670020120724, 0.103621730382294, 0.0784708249496982, 0.0110663983903421,
    #              0.0543259557344064, 0.0563380281690141, 0.0734406438631791, 0.0160965794768612,
    #              0.0110663983903421, 0.028169014084507, 0.126760563380282, 0.0583501006036217, 0.062374245472837,
    #              0.0372233400402414, 0.17102615694165]


    myDevList = getDeviationFromMean(myFreqList, args[dict['name']], args[dict['html_resource']])  # saves out inputname_deviation_from_mean.csv
    relativeFrequencyHeatMap(args[dict['name']], myFreqList,args[dict['pdfexport']], args[dict['html_resource']])
    deviationFromMeanHeatMap(args[dict['name']], myDevList,args[dict['pdfexport']], args[dict['html_resource']])
    plotPCA(args[dict['name']], myFreqList,args[dict['pdfexport']], args[dict['html_resource']])
    createHTML(args[dict['name']], args[dict['html_file']], myFreqList, myDevList)  # assumes imgs are heatmap.png, dheatmap.png, vapPCA.png and already in htmlresource

def contigs(args,dict):
    #argdict = {'name': 2, 'pdfexport': 3, 'contigs': 4, 'html_file': 5, 'html_resource': 6}

    shutil.copyfile(args[dict['contigs']], args[dict['name']]+".fa")

    

    contigTranslation(args[dict['name']])
    myCountList = HMMerMotifSearch(args[dict['name']])
    myFreqList = relativeFrequencyTable(myCountList, args[dict['name']],
                                        args[dict['html_resource']])  # saves out inputname_relative_frequncy.csv
    # myFreqList = [0.111670020120724, 0.103621730382294, 0.0784708249496982, 0.0110663983903421,
    #              0.0543259557344064, 0.0563380281690141, 0.0734406438631791, 0.0160965794768612,
    #              0.0110663983903421, 0.028169014084507, 0.126760563380282, 0.0583501006036217, 0.062374245472837,
    #              0.0372233400402414, 0.17102615694165]


    myDevList = getDeviationFromMean(myFreqList, args[dict['name']],
                                     args[dict['html_resource']])  # saves out inputname_deviation_from_mean.csv
    relativeFrequencyHeatMap(args[dict['name']], myFreqList, args[dict['pdfexport']], args[dict['html_resource']])
    deviationFromMeanHeatMap(args[dict['name']], myDevList, args[dict['pdfexport']], args[dict['html_resource']])
    plotPCA(args[dict['name']], myFreqList, args[dict['pdfexport']], args[dict['html_resource']])
    createHTML(args[dict['name']], args[dict['html_file']], myFreqList,
               myDevList)  # assumes imgs are heatmap.png, dheatmap.png, vapPCA.png and already in htmlresource


def genomicProcess(inputname, exportpdf, forwardFN, reverseFN, htmlfile, htmlresource):
    assembleWithVelvet(inputname,forwardFN,reverseFN)
    contigTranslation(inputname)
    myCountList = HMMerMotifSearch(inputname)
    myFreqList = relativeFrequencyTable(myCountList, inputname, htmlresource)  # saves out inputname_relative_frequncy.csv
    #myFreqList = [0.111670020120724, 0.103621730382294, 0.0784708249496982, 0.0110663983903421,
    #              0.0543259557344064, 0.0563380281690141, 0.0734406438631791, 0.0160965794768612,
    #              0.0110663983903421, 0.028169014084507, 0.126760563380282, 0.0583501006036217, 0.062374245472837,
    #              0.0372233400402414, 0.17102615694165]


    myDevList = getDeviationFromMean(myFreqList, inputname,htmlresource)  # saves out inputname_deviation_from_mean.csv

    relativeFrequencyHeatMap(inputname, myFreqList, exportpdf, htmlresource)
    deviationFromMeanHeatMap(inputname, myDevList, exportpdf, htmlresource)
    plotPCA(inputname, myFreqList, exportpdf, htmlresource)
    createHTML(inputname, htmlfile, myFreqList,myDevList)  # assumes imgs are heatmap.png, dheatmap.png, vapPCA.png and already in htmlresource
    return



if __name__ == "__main__":
    #contigTranslation('Tcongo')
    #contigTranslation('Test')
    #newHMMerMotifSearch('Test')
    #HMMerMotifSearch('Tcongo')
    #sys.exit()


    myFreqList = [0.111670020120724, 0.103621730382294, 0.0784708249496982, 0.0110663983903421,
                  0.0543259557344064, 0.0563380281690141, 0.0734406438631791, 0.0160965794768612,
                  0.0110663983903421, 0.028169014084507, 0.126760563380282, 0.0583501006036217, 0.062374245472837,
                  0.0372233400402414, 0.17102615694165]
    myDevList = [0.000790026,0.0073109,-0.001151769,-0.004502933,-0.013687421,-0.016159773,0.021689891,
                 0.007863809,-0.003133585,-0.001111709,-0.01313879,0.0036997,-0.00935284,0.005640693,0.015243802]

    relativeFrequencyHeatMap('test', myFreqList, "PDF_Yes","results")
    deviationFromMeanHeatMap('test', myDevList, "PDF_Yes","results")
    plotPCA('test',myFreqList,"PDF_Yes","results")

    createHTML('test',"results/test.html", myFreqList, myDevList)
    #contigTranslation("Test")
    #myCountList = HMMerMotifSearch("Test")


    sys.exit()
