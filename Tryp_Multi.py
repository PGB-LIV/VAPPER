"""
 * Copyright 2018 University of Liverpool
 * Author: John Heap, Computational Biology Facility, UoL
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
import pandas as pd
import Tryp_G
import Tryp_T
import Tryp_V
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import glob

readFileList = []

def findpairedReadFiles(path):
    pairList = []
    fileList = glob.glob(path+"/*.*") # returns list
    print(fileList)
    #look in each filename to see if a 1 or 2 is near the end.
    for name in fileList:
        p1 = name.rfind('1')
        if p1 != -1:
            pairedname = name[:p1] + '2' + name[p1+1:]
            #print(pairedname)
            #check if in list
            if pairedname in fileList:
                pairList.append(name)
                pairList.append(pairedname)
    n = int(len(pairList) / 2)
    print("Assuming these are paired files:")
    for j in range(n):
        print(pairList[j * 2] + "," + pairList[j * 2 + 1])
    return pairList

def findContigFiles(path):
    contigList = []
    fileList = glob.glob(path+'/*.fa')
    for name in fileList:
        isFa = name.rfind('.fa')
        if isFa != -1:
            contigList.append(name)
    print(contigList)
    return contigList

def saveListsToCSV(name,freqLists,devLists,htmlresource):  #save freqLists and devLists to CSVs
    df_f = pd.DataFrame({'Phylotype': Tryp_G.pList})
    df_d = pd.DataFrame({'Phylotype': Tryp_G.pList})
    for j in range(len(freqLists)):
        df_f[freqLists[j][0]] = pd.Series(freqLists[j][1:16], index=df_f.index[0:15])
        df_d[devLists[j][0]] = pd.Series(devLists[j][1:16], index=df_f.index[0:15])

    j_fname = htmlresource + "/" + name + "_relative_frequency.csv"
    df_f.to_csv(j_fname, index=False)
    j_fname = htmlresource + "/" + name + "_deviation_from_mean.csv"
    df_d.to_csv(j_fname, index=False)

def saveTransListsToCSV(name,relFreqLists,relWeightLists,htmlresource):
    df_f = pd.DataFrame({'Phylotype': Tryp_G.pList})
    df_w = pd.DataFrame({'Phylotype': Tryp_G.pList})
    for j in range(len(relFreqLists)):
        df_f[relFreqLists[j][0]] = pd.Series(relFreqLists[j][1:16], index=df_f.index[0:15])
        df_w[relWeightLists[j][0]] = pd.Series(relWeightLists[j][1:16], index=df_f.index[0:15])
    j_fname = htmlresource + "/" + name + "_t_relative_frequency.csv"
    df_f.to_csv(j_fname, index=False)

    j_fname = htmlresource + "/" + name + "_t_weighted_frequecy.csv"
    df_w.to_csv(j_fname, index=False)



def createMultiHTML(name,htmlfn,htmlresource,freqLists,devLists):
    #assumes imgs are heatmap.png, dheatmap.png, vapPCA.png and already in htmlresource
    htmlString = r"<html><title>T.congolense VAP</title><body><div style='text-align:center'><h2><i>Trypanosoma congolense</i> Variant Antigen Profile</h2><h3>"
    htmlString += name
    htmlString += r"<br/>Genomic Analysis</h3>"
    htmlString += "<p style = 'margin-left:23%; margin-right:23%'>Table Legend: Variant Antigen Profiles of <i>Trypanosoma " \
                  "congolense</i> estimated as the phylotype proportion and as the deviation from the mean across the sample cohort.<br>" \
                  "Data was produced with the 'Variant Antigen Profiler' (Silva Pereira and Jackson, 2018).</p>"
    htmlString += r"<style> table, th, tr, td {border: 1px solid black; border-collapse: collapse;}</style>"

    header = r"<table style='width:50%;margin-left:25%;text-align:center'><tr><th>Phylotype</th>"
    for j in range(len(freqLists)):
        header += r"<th>" + str(freqLists[j][0]) + "</th>"
    htmlString += "</tr>\n" + header
    tabString = ""
    for i in range(15):
        tabString += "<tr><td>phy" + str(i + 1) + "</td>"
        for j in range(len(freqLists)):
            f = format(freqLists[j][i + 1], '.4f')
            tabString += "<td>" + f + "</td>"
        tabString += "</tr>\n"

    htmlString += tabString + "</table><br><br>\n"

    htmlString += "Deviation from the mean of the 15 phylotypes within the variant repertoire.\n"
    htmlString += r"<style> table, th, tr, td {border: 1px solid black; border-collapse: collapse;}</style>"
    header = r"<table style='width:50%;margin-left:25%;text-align:center'><tr><th>Phylotype</th>"
    for j in range(len(devLists)):
        header += r"<th>" + str(devLists[j][0]) + "</th>"
    htmlString += "</tr>\n" + header
    tabString = ""
    for i in range(15):
        tabString += "<tr><td>phy" + str(i + 1) + "</td>"
        for j in range(len(devLists)):
            f = format(devLists[j][i + 1], '.4f')
            tabString += "<td>" + f + "</td>"
        tabString += "</tr>\n"

    htmlString += tabString + "</table><br><br><br><br><br>\n"

    htmlString += r"<h3>The Variation Heat Map and Dendrogram</h3><p>The absolute phylotype variation in the sample compared to model dataset.</p>"
    imgString = r"<img src = 'heatmap.png' alt='Variation Heatmap' style='max-width:100%'><br><br>"
    htmlString += imgString

    htmlString += r"<br><br><br><br><h3>The Deviation Heat Map and Dendrogram</h3><p>The phylotype variation expressed as the deviation from your sample mean compared to the model dataset</p>"
    imgString = r"<img src = 'dheatmap.png' alt='Deviation Heatmap' style='max-width:100%'><br><br>"
    htmlString += imgString

    htmlString += r"<br><br><br><br><h3>The Variation PCA plot</h3><p>PCA analysis corresponding to absolute variation. Colour coded according to location</p>"
    imgString = r"<img src = 'vapPCA.png' alt='PCA Analysis' style='max-width:100%'><br><br>"
    htmlString += imgString + r"</div></body></html>"

    with open(htmlresource + '/' + htmlfn, "w") as htmlfile:
        htmlfile.write(htmlString)


def createHTML_T(name,htmlfn,htmlresource,weightLists):
    #assumes imgs are heatmap.png, dheatmap.png, vapPCA.png and already in htmlresource
    htmlString = r"<html><title>T.congolense VAP</title><body><div style='text-align:center'><h2><i>Trypanosoma congolense</i> Variant Antigen Profile</h2><h3>"
    htmlString += name
    htmlString += r"<br>Transcriptomic Analysis</h3></p>"
    htmlString += "<p style = 'margin-left:20%; margin-right:20%'>Table Legend: Variant Antigen Profiles of a transcriptome of <i>Trypanosoma congolense</i> estimated as phylotype proportion. " \
                  "Weighted frequency refers to the phylotype proportion based transcript abundance. " \
                  "Data was produced with the 'Variant Antigen Profiler' (Silva Pereira and Jackson, 2018).</p> "
    htmlString += r"<style> table, th, tr, td {border: 1px solid black; border-collapse: collapse;}</style>"


    header = r"<table style='width:50%;margin-left:25%;text-align:center'><tr><th>Phylotype</th>"
    for j in range (len(weightLists)):
        header += r"<th>" + str(weightLists[j][0]) + "</th>"
    htmlString += "</tr>\n" + header
    tabString = ""
    for i in range(15):
        tabString += "<tr><td>phy" + str(i + 1) + "</td>"
        for j in range(len(weightLists)):
            f = format(weightLists[j][i + 1], '.4f')
            tabString += "<td>" + f + "</td>"
        tabString += "</tr>\n"

    htmlString += tabString + "</table><br><br><br><br><br>"
    htmlString += r"<p> <h3>Stacked Bar chart of Phylotype Frequency</h3> The 'weighted' relative frequency of each phylotype alongside the VAP of selected strain.</p>"
    imgString = r"<img src = 'stackedbar.png' alt='Stacked bar chart of phylotype variation' style='max-width:100%'><br><br>"
    htmlString += imgString

    with open(htmlresource+'/'+htmlfn, "w") as htmlfile:
        htmlfile.write(htmlString)

def createMultiStackedBar(name,relfreqLists,strain,pdf,html_resource):
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
    fig, ax = plt.subplots()
    if strain == "Tc148":
        VAPtable = VAP_148
        VAPname='Tc148\nGenome VAP'
    if strain == "IL3000":
        VAPtable = VAP_IL3000
        VAPname= 'IL3000\nGenome VAP'
    n = len(relfreqLists)
    width = 2*(1.0/((n+1)*2))

    #width = 0.35  # the width of the bars: can also be len(x) sequence
    plots = []
    print(len(relfreqLists))
    vpos = 0
    for j in range(len(relfreqLists)):
        fpos = 0
        for p in range(0, 15):
            tp = plt.bar(j, relfreqLists[j][p+1], width, color= palette[p], bottom = fpos)
            fpos +=relfreqLists[j][p+1]
            #plots.append(tp)

    for p in range(0, 15):
        tp = plt.bar(len(relfreqLists), VAPtable[p], width, color=palette[p], bottom=vpos)
        vpos += VAPtable[p]
        plots.append(tp)
    tickList = []
    for j in range(len(relfreqLists)):
        tickList.append(relfreqLists[j][0])
    tickList.append(VAPname)
    plt.xticks(range(len(relfreqLists)+1),tickList)
    #plt.xticks([0,len(relfreqLists)+1],[name,VAPname])
    plt.legend(plots[::-1],['p15','p14','p13','p12','p11','p10','p9','p8','p7','p6','p5','p4','p3','p2','p1'], loc = (0.975,0.0125))
    title = "Figure Legend: The transcriptomic Variant Antigen Profile of $\itTrypanosoma$ $\itcongolense$ estimated as phylotype " \
            "proportion adjusted for transcript abundance and the reference genomic Variant Antigen Profile. " \
            "\nData was produced with the 'Variant Antigen Profiler' (Silva Pereira and Jackson, 2018)."
    #plt.title(title, wrap="True")
    #plt.text(-0.2, -0.05, title, va="top", transform=ax.transAxes, wrap="True")
    plt.text(-0.3, -0.15, title, va="top", wrap="True")
    plt.tight_layout(pad=4.5)
    plt.subplots_adjust(bottom = 0.3,top=0.99,left=0.125,right=0.9,hspace=0.2,wspace=0.2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.savefig(html_resource + "/stackedbar.png")
    if pdf == 'PDF_Yes':
        plt.savefig(html_resource + "/stackedbar.pdf")
    plt.show()





def dotheplots(tmpname,freqLists,devLists,pdf,htmlresource):

    dir_path = os.path.dirname(os.path.realpath(__file__))
    j_fname = dir_path + "/data/congodata.csv"
    congo_df = pd.read_csv(j_fname)

    j_fname = dir_path + "/data/congodata_deviationfromthemean.csv"
    congo_dev_df = pd.read_csv(j_fname)

    for j in range(len(freqLists)):
        congo_df.loc[congo_df.index.max() + 1] = freqLists[j]   #[0:16]
        congo_dev_df.loc[congo_df.index.max() + 1] = devLists[j]    #[0:16]

    Tryp_G.plotPCA(congo_df, pdf, htmlresource)

    ysize = len(congo_df) * 20 / 97.0  # make vertical size equivlanet 20' is ok for 97.
    congo_dev_df.set_index('Strain', inplace=True)

    #congo_df.drop('Colour', axis=1, inplace=True)       #congo_df now ready for inclusion of freqlist
    congo_dev_df.drop('Colour', axis=1, inplace=True)


    cg = sns.clustermap(congo_df, method='ward', cmap="RdBu_r", col_cluster=False, yticklabels=congo_df.index.values, figsize = (10,ysize))
    plt.setp(cg.ax_heatmap.yaxis.get_ticklabels(), rotation=0, fontsize=8)  # get y labels printed horizontally
    ax = cg.ax_heatmap
    title = "Variant Antigen Profiles of $\itTrypanosoma$ $\itcongolense$ estimated as the phylotype proportion across the\nsample cohort. "
    title += "Dendrogram reflects the relationships amongst the VSG repertoires of each strain. "
    title += "Strains\nwere isolated from multiple African countries as described in Silva Pereira et al. (2018)."
    title += "\nData was produced with the 'Variant Antigen Profiler' (Silva Pereira and Jackson, 2018)."
    ax.text(-0.15, -0.05, title, va="top", wrap="True", transform=ax.transAxes)
    plt.savefig(htmlresource + "/heatmap.png", bbox_inches='tight')
    if pdf == 'PDF_Yes':
        plt.savefig('htmlresource' + "/heatmap.pdf", bbox_inches='tight')



    cg = sns.clustermap(congo_dev_df, method='ward', cmap="RdBu_r", col_cluster=False, yticklabels=congo_df.index.values, figsize = (10,ysize))
    plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0, fontsize=8)  # get y labels printed horizontally
    ax = cg.ax_heatmap
    title = "Variant Antigen Profiles of $\itTrypanosoma$ $\itcongolense$ expressed as the deviation from the mean phylotypes "
    title += "\nproportions of the sample cohort. Dendrogram reflects the relationships amongst the VSG repertoires of "
    title += "each \nstrain. Strains were isolated from multiple African countries as described in Silva Pereira et al. (2018)."
    title += "\nData was produced with the 'Variant Antigen Profiler' (Silva Pereira and Jackson, 2018)."

    ax.text(-0.2, -0.05, title, va="top", transform=ax.transAxes, wrap="True")
    plt.savefig(htmlresource+ "/dheatmap.png", bbox_inches='tight')
    if pdf == 'PDF_Yes':
        plt.savefig(htmlresource + "/dheatmap.pdf", bbox_inches='tight')


def multi_G_Assembly(dict):
    #set up congodata.csv ready fro inclusion of freqlist
    freqLists= []   #list of freqlists
    devLists = []   #list p deviation lists

    readFileList = findpairedReadFiles(dict['directory'])  # get readfile pairs from mypath
    n = int(len(readFileList)/2)
    for j in range(n):
        #forward = os.path.dirname(args[dict['directory']])+'/'+readFileList[j*2]
        forward = readFileList[j*2]
        reverse = readFileList[j*2+1]
        #reverse = os.path.dirname(args[dict['directory']])+'/'+readFileList[j*2+1]
        tmpname = forward.split('.')[0]
        Tryp_G.assembleWithVelvet(tmpname,dict['kmers'], dict['inslen'],dict['covcut'], forward,reverse)
        Tryp_G.contigTranslation(tmpname)
        myCountList = Tryp_G.HMMerMotifSearch(tmpname)
        if '/' in tmpname:
            tmpname = tmpname.split('/')[1]
        tmpname = tmpname.split('.')[0]
        freqLists.append(Tryp_G.relativeFrequencyTableNoSave(myCountList))
        devLists.append(Tryp_G.getDeviationFromMeanNoSave(freqLists[j]))
        freqLists[j].insert(0,tmpname)
        freqLists[j].append(tmpname)
        devLists[j].insert(0, tmpname)
        devLists[j].append(tmpname)
    saveListsToCSV(dict['name'],freqLists,devLists,dict['html_resource'])
    dotheplots(dict['name'],freqLists,devLists,dict['pdf'],dict['html_resource'])
    createMultiHTML(dict['name'],dict['html_file'],dict['html_resource'], freqLists,devLists)  # assumes imgs are heatmap.png, dheatmap.png, vapPCA.png and already in htmlresource
    print("Placing results in " + dict['html_resource'])



def multi_G_contigs(dict):
    # argdict = {'name': 2, 'pdfexport': 3, 'directory': 4, 'html_file': 5, 'html_resource': 6}
    freqLists = []
    devLists = []
    contigList = findContigFiles(dict['directory'])  # get readfile pairs from mypath
    for j in range(len(contigList)):
        contig = contigList[j]
        tmpname = contig.split('/')[1]
        tmpname = tmpname.split('.')[0]
        #print (tmpname)
        Tryp_G.contigTranslation(dict['directory']+"/"+tmpname)
        myCountList = Tryp_G.HMMerMotifSearch(dict['directory']+"/"+tmpname)
        freqLists.append(Tryp_G.relativeFrequencyTableNoSave(myCountList))
        devLists.append(Tryp_G.getDeviationFromMeanNoSave(freqLists[j]))
        freqLists[j].insert(0, tmpname)
        freqLists[j].append(tmpname)
        devLists[j].insert(0, tmpname)
        devLists[j].append(tmpname)
        #print(freqLists[j])
    saveListsToCSV(dict['name'],freqLists,devLists,dict['html_resource'])
    dotheplots(dict['name'],freqLists,devLists,dict['pdf'],dict['html_resource'])
    createMultiHTML(dict['name'], dict['name']+".html",dict['html_resource'], freqLists,devLists)  # assumes imgs are heatmap.png, dheatmap.png, vapPCA.png and already in htmlresource
    print("Placing results in " + dict['html_resource'])

def multi_T_process(dict):
    relFreqLists = []
    relWeightLists = []
    readFileList = findpairedReadFiles(dict['directory'])  # get readfile pairs from mypath
    n = int(len(readFileList) / 2)
    for j in range(n):
        # forward = os.path.dirname(args[dict['directory']])+'/'+readFileList[j*2]
        forward = readFileList[j * 2]
        reverse = readFileList[j * 2 + 1]

        tmpname = forward.split('/')
        tmpname = tmpname[-1]
        tmpname = tmpname.split('.')[0]
        print('tmpname = %s' %(tmpname))

        Tryp_T.transcriptMapping(tmpname, dict['strain'],forward,reverse)        #uses bowtie
        Tryp_T.processSamFiles(tmpname)                              #uses samtools
        Tryp_T.transcriptAbundance(tmpname,dict['strain'])                       #uses cufflinks -> ?.cuff/*.*
        cuff_df = Tryp_T.convertToFasta(tmpname,dict['strain'])
        countList, weightList = Tryp_T.HMMerMotifSearch(tmpname,dict['strain'], cuff_df)
        relFreqLists.append(Tryp_T.relativeFrequencyTableNoSave(countList))
        relWeightLists.append(Tryp_T.weightedFrequencyTableNoSave(weightList))
        relFreqLists[j].insert(0,tmpname)
        relWeightLists[j].insert(0,tmpname)

    saveTransListsToCSV(dict['name'],relFreqLists,relWeightLists,dict['html_resource'])
    createMultiStackedBar(dict['name'],relWeightLists, dict['strain'],dict['pdf'],dict['html_resource'])
    createHTML_T(dict['name'],dict['name']+".html", dict['html_resource'], relWeightLists)
    print("Placing results in " + dict['html_resource'])

def save_V_ListsToCSV(name, nameList, cogList, htmlresource):
    df_cog = cogList[0]  # get the first...
    for j in range(1,len(cogList)):
        nm = nameList[j]
        c = cogList[j][nm]
        df_cog[nm]=c
    df_cog = df_cog.drop(df_cog.columns[0], axis = 1)   #remove the extra index we got
    print(df_cog)
    j_fname = htmlresource + "/" + name + "_cogspresent.csv"
    df_cog.to_csv(j_fname, index = False)




def addtoCurrentDatabase(nameList,cogList):
    j_fname = "" + r"data/vivax/TvDatabase.csv"
    current_tv_df = pd.read_csv(j_fname)  # get current TVDatabase (we add to this as we go through the multisample process
    for j in range(len(cogList)):
        nm = nameList[j]
        cogAsList = cogList[j][nm].tolist()
        current_tv_df.loc[:, nm] = cogAsList  # add to current database
    #print(current_tv_df)
    #current_tv_df.to_csv("tv_Test.csv")
    return current_tv_df

def create_V_MultiClusterMap(tv_df,name,html_path,pdf=False):
    Tryp_V.createClusterMap(tv_df,name,html_path,pdf=False)
    return

def multi_V_Contigs(dict):
    cog_presenceList = []  # list of 'cogPresence.df's
    tmpnameList = []
    name = dict['name']
    htmlpath = dict['html_resource']
    dir_path = os.path.dirname(os.path.realpath(__file__))
    #j_fname = dir_path + r"/data/vivax/TvDatabase.csv"
    #current_tv_df = pd.read_csv(j_fname)  # get current TVDatabase (we add to this as we go through the multisample process

    vivaxPath = os.path.dirname(os.path.realpath(__file__)) + r"/data/vivax"
    contigList = findContigFiles(dict['directory'])  # get readfile pairs from mypath
    n = int(len(contigList))
    for j in range(n):
        # forward = os.path.dirname(args[dict['directory']])+'/'+readFileList[j*2]
        contig = contigList[j]
        #forward = readFileList[j * 2]
        #reverse = readFileList[j * 2 + 1]
        # reverse = os.path.dirname(args[dict['directory']])+'/'+readFileList[j*2+1]
        tmpname = contig.split('.')[0]

        print(tmpname)

        #Tryp_V.assembleWithVelvet(tmpname, args[dict['kmers']], args[dict['inslen']], args[dict['covcut']], forward,reverse)
        blastResult_df = Tryp_V.blastContigs(tmpname, vivaxPath + r"/Database/COGs.fas")
        orthPresence_df = Tryp_V.getCogsPresent(blastResult_df, tmpname, vivaxPath + r"/Database/COGlist.txt")
        binBlastResult_df = Tryp_V.blastContigs(tmpname, vivaxPath + r"/Database/Bin_2.fas")
        binPresence_df = Tryp_V.getCogsPresent(binBlastResult_df, tmpname, vivaxPath + r"/Database/binlist.txt")
        cogPresence_df = orthPresence_df.append(binPresence_df, ignore_index=True)
        tmpnameList.append(tmpname)
        cog_presenceList.append(cogPresence_df)  # add in cogPresence to list

    tv_df = addtoCurrentDatabase(tmpnameList, cog_presenceList)
    Tryp_V.create_V_MultiClusterMap(tv_df, tmpnameList, name, htmlpath, pdf=False)
    save_V_ListsToCSV(name, tmpnameList, cog_presenceList, htmlpath)
    Tryp_V.create_V_MultiHTML(name, name + ".html", htmlpath)
    print("Placing results in " + dict['html_resource'])
    return


def multi_V_Assembly(dict):
    cog_presenceList = []  # list of 'cogPresence.df's
    tmpnameList = []
    name = dict['name']
    htmlpath = dict['html_resource']
    dir_path = os.path.dirname(os.path.realpath(__file__))
    j_fname = dir_path + r"/data/vivax/TvDatabase.csv"
    current_tv_df = pd.read_csv(
        j_fname)  # get current TVDatabase (we add to this as we go through the multisample process

    vivaxPath = os.path.dirname(os.path.realpath(__file__)) + r"/data/vivax"
    readFileList = findpairedReadFiles(dict['directory'])  # get readfile pairs from mypath
    n = int(len(readFileList) / 2)
    for j in range(n):
        # forward = os.path.dirname(args[dict['directory']])+'/'+readFileList[j*2]
        forward = readFileList[j * 2]
        reverse = readFileList[j * 2 + 1]
        # reverse = os.path.dirname(args[dict['directory']])+'/'+readFileList[j*2+1]
        tmpname = forward.split('.')[0]

        print(tmpname)

        Tryp_V.assembleWithVelvet(tmpname, dict['kmers'], dict['inslen'], dict['covcut'], forward,reverse)
        blastResult_df = Tryp_V.blastContigs(tmpname, vivaxPath + r"/Database/COGs.fas")
        orthPresence_df = Tryp_V.getCogsPresent(blastResult_df, tmpname, vivaxPath + r"/Database/COGlist.txt")
        binBlastResult_df = Tryp_V.blastContigs(tmpname, vivaxPath + r"/Database/Bin_2.fas")
        binPresence_df = Tryp_V.getCogsPresent(binBlastResult_df, tmpname, vivaxPath + r"/Database/binlist.txt")
        cogPresence_df = orthPresence_df.append(binPresence_df, ignore_index=True)
        tmpnameList.append(tmpname)
        cog_presenceList.append(cogPresence_df)  # add in cogPresence to list

    tv_df = addtoCurrentDatabase(tmpnameList, cog_presenceList)
    Tryp_V.create_V_MultiClusterMap(tv_df, tmpnameList, name, htmlpath, pdf=False)
    save_V_ListsToCSV(name, tmpnameList, cog_presenceList, htmlpath)
    Tryp_V.create_V_MultiHTML(name, name + ".html", htmlpath)
    print("Placing results in " + dict['html_resource'])
    return


if __name__ == "__main__":
    print("ERROR: Tryp_Multi.py should only be called from within VAp.py")
    sys.exit()

