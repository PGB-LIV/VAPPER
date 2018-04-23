import matplotlib as mpl
mpl.use('Agg')

import subprocess
import shutil
import re
import pandas as pd
import os

import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

def assembleWithVelvet(name, kmers, inslen, covcut, fastq1name,fastq2name):
    #argString = "velveth " + name + "_k65 65 -shortPaired -fastq " + name + "_R1.fastq " + name + "_R2.fastq"
    argString = "velveth " + name + "_k"+ kmers+" "+ kmers + " -shortPaired -fastq " + fastq1name+" "+fastq2name
    print(argString)
    returncode = subprocess.call(argString, shell=True)
    if returncode != 0:
        return "Error in velveth"
    argString = "velvetg " + name + "_k"+kmers+" -exp_cov auto -ins_length "+inslen+" -clean yes -ins_length_sd 50 -min_pair_count 20"
    #argString = "velvetg " + name + "_k"+kmers+" -exp_cov auto -ins_length "+inslen+" -cov_cutoff "+covcut+" -clean yes -ins_length_sd 50 -min_pair_count 20"
    #argString = "velvetg " + name + "_k65 -exp_cov auto -ins_length 400 -cov_cutoff 5 -clean yes -ins_length_sd 50 -min_pair_count 20"+quietString
    print(argString)
    returncode = subprocess.call(argString, shell = True)
    if returncode != 0:
        return "Error in velvetg"
    shutil.copyfile(name + "_k"+kmers+"//contigs.fa",name + ".fa")  # my $namechange = "mv ".$input."_k65/contigs.fa ".$input.".fa";
    return "ok"


def blastContigs(test_name,database):
    print(test_name)
    print(database)
    #db_path = os.path.dirname(os.path.realpath(__file__))+database
    db_path = database
    argString = "blastn -db "+db_path+" -query "+test_name+".fa -outfmt 10 -out "+test_name+"_blast.txt"
    print (argString)
    returncode = subprocess.call(argString, shell = True)
    if returncode != 0:
        return "Error in blastall"
    blast_df = pd.read_csv(""+test_name+"_blast.txt")
    #print (blast_df)
    #if ($temp[2] >= 98 & & $temp[3] > 100 & & $temp[10] < 0.001){
    #'qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    blast_df.columns = ['qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue','bitscore']
    blastResult_df = blast_df[(blast_df['pident']>=98) & (blast_df['length'] > 100) & (blast_df['evalue']<0.001) ]
    blastResult_df = blastResult_df[['qaccver', 'saccver', 'pident']]   #query accession.version, subject accession.version, Percentage of identical matches

    return blastResult_df



def getCogsPresent(blastResult_df,strain,cogOrBinList):
    blastResult_df.sort_values('pident',axis = 0, ascending=False, inplace=True)
    nodeList = blastResult_df['qaccver'].tolist()
    cogList = blastResult_df['saccver'].tolist()
    cogSet = set(cogList)   #get unique values
    cogList = sorted(cogSet)    #sort them

    #print (cogList)
    #print (len(cogList))

    thereList = []
    dataList = []
    #dir_path = os.path.dirname(os.path.realpath(__file__))
    fname = cogOrBinList
    cnt = 0
    with open (fname) as f:
        for line in f:
            dataList.append(line.rstrip('\n\r '))
            if line.rstrip('\n\r ') in cogList:
                thereList.append('1')
                cnt = cnt+1
            else:
                thereList.append('0')

    #print (thereList)
    #print (cnt)
    data = {'Cog': dataList, strain: thereList}
    presence_df = pd.DataFrame(data)
    #print (presence_df)
    return presence_df

def addToCurrentData(cog_df, name):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    j_fname = dir_path + r"/data/vivax/TvDatabase.csv"
    tv_df = pd.read_csv(j_fname)

    cogList = cog_df[name].tolist()
    #cogList.insert(0,'Test')
    #print (len(tv_df))
    #print(len(cogList))

    #print(cogList)
    tv_df.loc[:,name]=cogList
    return tv_df




def createClusterMap(tv_df,name,html_path,pdfExport):
    #Retrieve Data
    dir_path = os.path.dirname(os.path.realpath(__file__))
    j_fname = dir_path+r"/data/vivax/geoTv.csv"
    geo_df = pd.read_csv(j_fname)
    geo_df.loc[len(geo_df)] =[name,name,'k']
    print(geo_df)

    myStrains = tv_df.columns.values.tolist()   #beware first entry is COG
    myStrains = myStrains[1:]
    print(myStrains)
    myPal = []
    for s in myStrains:
        col = geo_df[(geo_df['Strain'] == s)]['colour'].tolist()
        myPal.append(col[0])
    print(myPal)
    mycogmap = ['skyblue', 'orangered']  # blue absent,red present
    tv_df.set_index('COG', inplace=True)
    tv_df = tv_df[tv_df.columns].astype(float)

    cg = sns.clustermap(tv_df, method='ward', col_colors=myPal, cmap=mycogmap, yticklabels = 1500, row_cluster=False, linewidths = 0)
    #cg = sns.clustermap(tv_df, method='ward', row_cluster=False, linewidths = 0)
    ax = cg.ax_heatmap
    #xasix ticks and labels.
    ax.xaxis.tick_top()     #set ticks at top
    newlabs = []

    labs = ax.xaxis.get_ticklabels()
    for i in range(0, len(labs)):
        print(labs[i])
        # labs[i].set_text("       "+labs[i].get_text())  #make enough room so label sits atop of col_color bars
        newlabs.append("       " + labs[i].get_text())
    ax.xaxis.set_ticklabels(newlabs)

    #labs = ax.xaxis.get_ticklabels()
    #for i in range(0, len(labs)):
    #    print(labs[i])
    #    labs[i].set_text("       "+labs[i].get_text())  #make enough room so label sits atop of col_color bars
    #    print(labs[i])
    #ax.xaxis.set_ticklabels(labs)
    plt.setp(cg.ax_heatmap.xaxis.get_ticklabels(), rotation=90)  # get x labels printed vertically

    cg.cax.set_visible(False)
    ax = cg.ax_heatmap
    ax.set_yticklabels("")
    ax.set_ylabel("")
    ax = cg.ax_heatmap
    # ax.set_title("Variant antigen profiles of T. vivax genomes.\nDendrogram reflects the VSG repertoire relationships of each strain inferred by the presence and absence of non-universal T. vivax VSG orthologs.", va = "top", wrap = "True")
    b = len(tv_df)
    print(b)
    title = "Figure Legend: The Variant Antigen Profiles of $\itTrypanosoma$ $\itvivax$ " \
            "showing the \ncombination of present and absent diagnostic cluster of VSG orthologs " \
            "across the sample cohort. \nDendrogram reflects the relationships amongst the VSG" \
            " repertoires of each strain. " \
            "Strains were isolated \nfrom multiple African countries as shown in the key.\nData was produced with the " \
            "'Variant Antigen Profiler' (Silva Pereira and Jackson, 2018)."

    ax.text(-1.5, len(tv_df) + 8,
            title,
            ha="left", va="top", wrap="True")
    col = cg.ax_col_dendrogram.get_position()
    cg.ax_col_dendrogram.set_position([col.x0, col.y0*1.08, col.width, col.height*1.1])


    legend_elements = [Patch(facecolor='orangered', label='COG Present'),
                       Patch (facecolor='skyblue', label='COG Absent'),
                       Patch(facecolor='w', label=''),
                       Patch (facecolor='b', label='Nigeria'),
                       Patch(facecolor = 'g', label='Uganda'),
                       Patch (facecolor='c', label='Gambia'),
                       Patch (facecolor='r', label='Ivory Coast'),
                       Patch(facecolor='m', label='Brazil'),
                       Patch(facecolor='k', label=name)]
    #legend_test = [[Patch(facecolor='orangered'),Patch(facecolor='r')],["test","test2"]]
    ax.legend(handles = legend_elements, bbox_to_anchor=(-0.3,1.2),loc = 'upper left')







    #plt.setp(cg.ax_heatmap.yaxis.get_ticklabels(), rotation=0 )  # get y labels printed horizontally
    # cg.dendrogram_col.linkage  # linkage matrix for columns
    # cg.dendrogram_row.linkage  # linkage matrix for rows
    # plt.savefig(r"results/" + name + "_heatmap.png")
    #plt.savefig(htmlresource + "/heatmap.png")
    #if pdf == 'PDF_Yes':
    #    plt.savefig(htmlresource + "/heatmap.pdf")
        # shutil.copyfile("heatmap.pdf",heatmapfn)  #
    #plt.legend()
    fname = html_path+"/"+name+"_clustermap.png"
    cg.savefig(fname)
    if pdfExport == 'PDF_Yes':
        fname = html_path + "/" + name + "_clustermap.pdf"
        cg.savefig(fname)
    #plt.show()

def createHTML(name,htmlfn,htmlPath):
    #assumes imgs are heatmap.png, dheatmap.png, vapPCA.png and already in htmlresource
    htmlString = r"<html><title>T.vivax VAP</title><body><div style='text-align:center'><h2><i>Trypanosoma vivax</i> Variant Antigen Profile</h2><h3>"
    htmlString += name


    htmlString += r'<p> <h3>The Heat Map and Dendrogram</h3></p>'
    imgString = r"<img src = '"+name+"_clustermap.png' alt='Cog Clustermap' style='max-width:100%'><br><br>"
    htmlString += imgString
    print(htmlString)

    with open(htmlfn, "w") as htmlfile:
        htmlfile.write(htmlString)


def vivax_assemble(args,dict):
    #argdict = {'name': 2, 'pdfexport': 3, 'kmers': 4, 'inslen': 5, 'covcut': 6, 'forward': 7, 'reverse': 8, 'html_file': 9,'html_resource': 10}
    #assembleWithVelvet("V2_Test", '65', '400', '5', "data/TviBrRp.1.clean", "data/TviBrRp.2.clean")

    vivaxPath = os.path.dirname(os.path.realpath(__file__))+r"/data/vivax"
    assembleWithVelvet(args[dict['name']], args[dict['kmers']], args[dict['inslen']], args[dict['covcut']],
                       args[dict['forward']], args[dict['reverse']])
    blastResult_df = blastContigs(args[dict['name']], vivaxPath+r"/Database/COGs.fas")
    orthPresence_df = getCogsPresent(blastResult_df, args[dict['name']], vivaxPath+r"/Database/COGlist.txt")

    binBlastResult_df = blastContigs(args[dict['name']], vivaxPath+r"/Database/Bin_2.fas")
    binPresence_df = getCogsPresent(binBlastResult_df, args[dict['name']], vivaxPath+r"/Database/binlist.txt")
    cogPresence_df = orthPresence_df.append(binPresence_df, ignore_index=True)

    fname = args[dict['html_resource']] +args[dict['name']]+"_cogspresent.csv"
    cogPresence_df.to_csv(fname)
    current_df = addToCurrentData(cogPresence_df,args[dict['name']])  # load in Tvdatabase and add cogPresence column to it.
    createClusterMap(current_df, args[dict['name']],args[dict['html_resource']],args[dict['pdfexport']])
    createHTML(args[dict['name']],args[dict['html_file']],args[dict['html_resource']])

def test_cluster(args,dict):
    print ("name: %s",args[dict['name']])
    cogPresence_df = pd.read_csv("test_cogspresent.csv")
    print(cogPresence_df)
    current_df = addToCurrentData(cogPresence_df,args[dict['name']])  # load in Tvdatabase and add cogPresence column to it.
    createClusterMap(current_df, args[dict['name']], args[dict['html_resource']], args[dict['pdfexport']])

def vivax_contigs(args,dict):
# argdict = {'name': 2, 'pdfexport': 3, 'contigs': 4, 'html_file': 5, 'html_resource': 6}

    #test_cluster(args,dict)
    #return;

    #subprocess.call('echo $PATH',shell = True)
    #sys.exit(1)

    vivaxPath = os.path.dirname(os.path.realpath(__file__))+r"/data/vivax"
    shutil.copyfile(args[dict['contigs']], args[dict['name']]+".fa")
    blastResult_df = blastContigs(args[dict['name']], vivaxPath+r"/Database/COGs.fas")
    orthPresence_df = getCogsPresent(blastResult_df, args[dict['name']], vivaxPath+r"/Database/COGlist.txt")

    binBlastResult_df = blastContigs(args[dict['name']], vivaxPath+r"/Database/Bin_2.fas")
    binPresence_df = getCogsPresent(binBlastResult_df, args[dict['name']], vivaxPath+r"/Database/binlist.txt")
    cogPresence_df = orthPresence_df.append(binPresence_df, ignore_index=True)

    fname = args[dict['html_resource']]+r"/"+ args[dict['name']]+"_cogspresent.csv"
    cogPresence_df.to_csv(fname)
    current_df = addToCurrentData(cogPresence_df,args[dict['name']])  # load in Tvdatabase and add cogPresence column to it.
    createClusterMap(current_df, args[dict['name']], args[dict['html_resource']], args[dict['pdfexport']])
    createHTML(args[dict['name']],args[dict['html_file']],args[dict['html_resource']])

if __name__ == "__main__":
    #assembleWithVelvet("V2_Test",'65','400', '5',"data/TviBrRp.1.clean","data/TviBrRp.2.clean")
    #assembleWithVelvet("V2_Test",'65','400', '5',"data/Tv493.1","data/Tv493.2")
    #blastResult_df=blastContigs("V2_Test",r"/Database/COGs.fas")
    cogPresence_df = pd.read_csv("test_cogspresent.csv")
    print(cogPresence_df)
    current_df = addToCurrentData(cogPresence_df,"vTest")  # load in Tvdatabase and add cogPresence column to it.
    createClusterMap(current_df, "vTest", "sausages","no")
    createHTML("vTest","vTest.html","sausages")
    sys.exit()

    blastResult_df=blastContigs("Tv493",r"/Database/COGs.fas")
    orthPresence_df = getCogsPresent(blastResult_df,"Tv493",r"/Database/COGlist.txt")

    #binBlastResult_df=blastContigs("V2_Test",r"/Database/Bin_2.fas")

    binBlastResult_df=blastContigs("Tv493",r"/Database/Bin_2.fas")
    binPresence_df = getCogsPresent(binBlastResult_df,"Tv493",r"/Database/binlist.txt")
    cogPresence_df = orthPresence_df.append(binPresence_df, ignore_index=True)
    #now do the next bit?
    current_df = addToCurrentData(cogPresence_df)  # load in Tvdatabase and add cogPresence column to it.
    createClusterMap(current_df,'Test',dict['html_resource'] ,dict['pdfexport'])


    #print(cogPresence_df)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    fname = dir_path+r"/results/V2_TestPresence.csv"
    #fnameb = dir_path+r"/results/V2_Test_blastOrth.csv"
    #fnameb_bin = dir_path+r"/results/V2_Test_blastBin.csv"
    #binBlastResult_df.to_csv(fnameb_bin)
    #blastResult_df.to_csv(fnameb)

    #cogPresence_df.to_csv(fname)
    cogPresence_df = pd.read_csv(fname)

    current_df = addToCurrentData(cogPresence_df)       #load in Tvdatabase and add cogPresence column to it.




