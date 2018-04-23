"""
From Sara's VAP_TP.sh
readname=$1 #strain or input file name
res=$2 #forward transcript reads file file (fastq)
re=$3 #reverse transcript reads file file (fastq)

bowtie2 -x Reference/IL3000 -1 $2 -2 $3  -S $1.sam 2> log.txt

samtools view -bS $1.sam > $1.bam
samtools sort $1.bam $1.sorted
samtools index $1.sorted.bam $1.sorted.bai


cufflinks -G Reference/IL3000.gtf -o $1.cuff -u -p 8 $1.sorted.bam
"""

import subprocess
import pandas as pd
import re
import os
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

pList = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10', 'P11', 'P12', 'P13', 'P14', 'P15']
quietString = "" #"">> Vap_log.txt 2>&1"
def transcriptMapping(inputname, strain, forwardFN,reverseFN):
    #where is our Reference data -
    dir_path = os.path.dirname(os.path.realpath(__file__))
    refName = dir_path+"/data/Reference/Tc148" #default
    if strain == "Tc148":
        refName = dir_path+"/data/Reference/Tc148"
    if strain == "IL3000":
        refName = dir_path+"/data/Reference/IL3000"
    #argString = "bowtie2 -x Refe4rence/IL3000 -1 data/"+forwardFN+" -2 data/"+reverseFN+" -S "+inputname+".sam"    #>log.txt
    #argString = "bowtie2 -x Reference/Tc148 -1 data/"+forwardFN+" -2 data/"+reverseFN+" -S "+inputname+".sam"    #>log.txt
    argString = "bowtie2 -x "+refName+" -1 "+forwardFN+" -2 "+reverseFN+" -S "+inputname+".sam"+quietString    #>log.txt
    #print(argString)
    returncode = subprocess.call(argString, shell=True)

def processSamFiles(inputname):
    #debug use a mapping sam file we have already found
    #dir_path = os.path.dirname(os.path.realpath(__file__))
    #bugName = dir_path+"/data/T_Test" #defasult

    cur_path = os.getcwd()
    samName = cur_path+"/"+inputname

    #argString = "samtools view -bS "+bugName+" > "+inputname+".bam"
    argString = "samtools view -bS "+inputname+".sam > "+samName+".bam"+quietString
    #print(argString)
    returncode = subprocess.call(argString, shell=True)


    #argString = "samtools sort "+bugName+" -o "+inputname+".sorted"
    argString = "samtools sort "+samName+".bam -o "+samName+".sorted"+quietString
    #print("argstring = "+argString)
    returncode = subprocess.call(argString, shell=True)

    #argString = "samtools index "+bugName+".sorted "+inputname+".sorted.bai"
    argString = "samtools index "+samName+".sorted "+samName+".sorted.bai"+quietString
    #print("argstring = " + argString)
    returncode = subprocess.call(argString, shell=True)




def transcriptAbundance(inputname, strain):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    refName = dir_path + "/data/Reference/ORFAnnotation.gtf"  # defasult
    if strain == "Tc148":
        refName = dir_path + "/data/Reference/ORFAnnotation.gtf"
    if strain == "IL3000":
        refName = dir_path + "/data/Reference/IL3000.gtf"
    #argString = "cufflinks -G Reference/IL3000.gtf -o "+inputname+".cuff -u -p 8 "+inputname+".sorted"
    #argString = "cufflinks -G Reference/ORFAnnotation.gtf -o "+inputname+".cuff -u -p 8 "+inputname+".sorted"
    argString = "cufflinks -q -G "+refName+" -o "+inputname+".cuff -u -p 8 "+inputname+".sorted"+quietString
    returncode = subprocess.call(argString, shell = True)


def convertToFasta(inputName, strain):  #equivalent to Sara's awk scripte
    dir_path = os.path.dirname(os.path.realpath(__file__))
    refName = dir_path + "/data/Reference/ORFAnnotation.gtf"  # default
    if strain == "Tc148":
        refName = dir_path + "/data/Reference/148_prot.fasta"
    if strain == "IL3000":
        refName = dir_path + "data/Reference/IL3000_prot.fasta"

    cuff_df = pd.read_csv(inputName+".cuff/genes.fpkm_tracking", sep='\t')
    cuff_df = cuff_df[(cuff_df['FPKM'] > 0)]
    cuff_df.to_csv("cuffTest.csv")
    gene_id_List = cuff_df['gene_id'].tolist()

    #print(gene_id_List)
    #print ("Found from 8880="+str(found))

    # need to load in IL3000_prot.fasta
    # for each line with >TcIL3000_1_1940
    # search within cuff_df[gene_id] for match
    # add it to the outfile. (need to save it as used by hmmer later
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
                else:
                    line = ref.readline()
            else:
                line =ref.readline()
    ref.close()
    print(str(len(gene_id_List))+":"+str(number)+" from "+str(all))
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
            m = re.search(regex, line)
            if m:
                outList.append(""+m.group())
                n += 1
            if re.search(r"inclusion", line):
                print("inclusion threshold reached")
                break
        compoundList.append(outList)
        lineCounts.append(n)
        hmmResult.close()
        os.remove("Phy" + m + ".out")

    #print(lineCounts)

    #print(cuff_df)
    concatGroups = [1, 2, 1, 3, 1, 1, 1, 2, 3, 2, 2, 1, 4, 1, 3]
    countList = []
    weightList = []
    countIndex = 0
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
    #print(countList)
    #print("--------")
    #print(weightList)
    #print("--------")
    os.remove(name + "_6frame.fas")
    return countList,weightList

def relativeFrequencyTable(countList, name, htmlresource):
    relFreqList = []
    c = float(countList[15])
    for i in range(0, 15):
        relFreqList.append(countList[i] / c)

    data = {'Phylotype': pList, 'Relative Frequency': relFreqList}
    relFreq_df = pd.DataFrame(data)
    j_fname = htmlresource+ "/" + name + "_t_relative_frequency.csv"
    relFreq_df.to_csv(j_fname)
    return relFreqList  # 0-14 = p1-p15 counts [15] = total counts


def weightedFrequencyTable(countList, name, htmlresource):
    relFreqList = []
    c = float(countList[15])
    for i in range(0, 15):
        relFreqList.append(countList[i] / c)

    data = {'Phylotype': pList, 'Weighted Frequency': relFreqList}
    relFreq_df = pd.DataFrame(data)
    j_fname = htmlresource+ "/" + name + "_t_weighted_frequency.csv"
    relFreq_df.to_csv(j_fname)
    return relFreqList  # 0-14 = p1-p15 counts [15] = total counts



def createStackedBar(name,freqList,strain,pdf,html_resource):
    palette = ["#0000ff", "#6495ed", "#00ffff", "#caff70",
               "#228b22", "#528b8b", "#00ff00", "#a52a2a",
               "#ff0000", "#ffff00", "#ffa500", "#ff1493",
               "#9400d3", "#bebebe", "#000000", "#ff00ff"]

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
            "\nData was produced with the 'Variant Antigen Profiler' (Silva Pereira and Jackson, 2018)."
    #plt.title(title, wrap="True")
    #plt.text(-0.2, -0.05, title, va="top", transform=ax.transAxes, wrap="True")
    plt.text(-0.3, -0.15, title, va="top", wrap="True")
    plt.tight_layout(pad=1.5)
    plt.subplots_adjust(bottom = 0.3,top=0.99,left=0.125,right=0.9,hspace=0.2,wspace=0.2)

    plt.savefig(html_resource + "/stackedbar.png")
    if pdf == 'PDF_Yes':
        plt.savefig(html_resource + "/stackedbar.pdf")
    #plt.show()


def createHTML(name,htmlfn,htmlresource,freqList,weightList):
    #assumes imgs are heatmap.png, dheatmap.png, vapPCA.png and already in htmlresource
    htmlString = r"<html><title>T.congolense VAP</title><body><div style='text-align:center'><h2><i>Trypanosoma congolense</i> Variant Antigen Profile</h2><h3>"
    htmlString += name
    htmlString += r"<br>Transcriptomic Analysis</h3></p>"
    htmlString += "<p style = 'margin-left:20%; margin-right:20%'>Table Legend: Variant Antigen Profiles of a transcriptome of <i>Trypanosoma congolense</i> estimated as phylotype proportion. " \
                  "Weighted frequency refers to the phylotype proportion based transcript abundance. " \
                  "Data was produced with the 'Variant Antigen Profiler' (Silva Pereira and Jackson, 2018).</p> "
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

#    htmlString += r"<p><h3>The Deviation Heat Map and Dendogram</h3>The phylotype variation expressed as the deviation from your sample mean compared to the model dataset</p>"
#    imgString = r"<img src = 'dheatmap.png' alt='Deviation Heatmap' style='max-width:100%'><br><br>"
#    htmlString += imgString

#    htmlString += r"<p><h3>The Variation PCA plot</h3>PCA analysis corresponding to absolute variation. Colour coded according to location</p>"
#    imgString = r"<img src = 'vapPCA.png' alt='PCA Analysis' style='max-width:100%'><br><br>"
#    htmlString += imgString + r"</div></body></html>"

    with open(htmlfn, "w") as htmlfile:
        htmlfile.write(htmlString)

#argdict = {'name':2, 'pdfexport': 3, 'strain': 4, 'forward': 5, 'reverse': 6, 'html_file': 7, 'html_resource': 8}
def transcriptomicProcess(args,dict):
    transcriptMapping(args[dict['name']], args[dict['strain']], args[dict['forward']], args[dict['reverse']])        #uses bowtie
    processSamFiles(args[dict['name']])                              #uses samtools
    transcriptAbundance(args[dict['name']],args[dict['strain']])                          #uses cufflinks -> ?.cuff/*.*
    cuff_df = convertToFasta(args[dict['name']],args[dict['strain']])
    countList, weightList = HMMerMotifSearch(args[dict['name']],args[dict['strain']], cuff_df)
    relFreqList = relativeFrequencyTable(countList,args[dict['name']],args[dict['html_resource']])
    relWeightList = weightedFrequencyTable(weightList,args[dict['name']],args[dict['html_resource']])
    createStackedBar(args[dict['name']],relWeightList, args[dict['strain']],args[dict['pdfexport']],args[dict['html_resource']])
    createHTML(args[dict['name']],args[dict['html_file']],args[dict['html_resource']], relFreqList, relWeightList)

if __name__ == "__main__":
    #print("Commencing Transcript Mapping")
    #transcriptMapping("T_Test", "Transcripts.1","Transcripts.2")
    #print("Processimg Sam Files")
    #processSamFiles("T_Test")
    #print("Assessing Transcript Abundance")
    #transcriptAbundance("T_Test")
    #print ("Converting to Fasta Subset")
    #cuff_df = convertToFasta("T_Test")
    #print("Commencing HMMer search")
    #countList, weightList = HMMerMotifSearch("T_Test",cuff_df)
    #relativeFrequencyTable(countList,'T_Test')
    #weightedFrequencyTable(weightList,'T_Test')
    relFreqList = [0.111842105,0.059210526,0.026315789,0.013157895,
                   0.006578947,0.013157895,0.032894737,0.019736842,
                   0.039473684,0.046052632,0.217105263,0.065789474,
                   0.151315789,0.059210526,0.138157895]

    relWeightList = [0.07532571,0.05900545,0.009601452,0.042357532,0.01236219,0.001675663,0.04109726,
                     0.097464248,0.057491666,0.05826875,0.279457473,0.070004772,0.065329007,0.085361298,0.045197529]

    createStackedBar('T_Test',relWeightList, 'Tc148','PDF_Yes','results')
    createHTML("t_test","results/t_test.html","results",relFreqList,relWeightList)
