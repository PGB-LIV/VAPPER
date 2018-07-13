# Trypanosoma-VAP
Variant Antigen Profiling for Trypanosoma congolense and Trypanosoma vivax

Introduction:
-------------
Trypanosomes are important human and veterinary parasites that cause potentially lethal blood 
infections and a chronic wasting disease (African trypanosomiasis). These organisms use antigenic 
variation to evade the host immune response, for which their genomes contain many hundreds of 
variable antigen genes. Making sense of this hypervariable repertoire is a major challenge in genome 
analysis and a bottleneck in vaccine development.

We have devised a method for quantitative analysis of antigenic diversity in systems data (genomes, 
transcriptomes, and proteomes) called Variant Antigen Profiling (VAP). VAP has great potential for 
understanding how antigenic diversity relates to clinical outcome, how antigen genes may be used 
as epidemiological markers of virulence, and in measuring gene expression during experimental 
infections.

We would like Variant Antigen Profiling to be the standard approach to study this phenomenon and 
we are making our methods available to the research community to maximize its impact.


Currently two species of trypanosoma are considered, namely T. congolense and T. vivax. The 
antigen variability of these species differ and are examined using different methods. 

T. congolense:
The Trypanosoma congolense variant antigen repertoire is divided into 15 clades or phylotypes. 
These phylotypes are present in any T. congolense isolate, but their relative abundance varies 
between strains. The purpose of the VAP is to accurately quantify antigen diversity in any T. 
congolense isolate by calculating the relative frequency of each phylotype. 

T.congolense can be examined at the Genomic or Transcriptomic level. (optional [-t] parameter to 
flag the transcriptomic pathway.) 

Genomic VAP 
At the Genomic level the tool takes raw paired NGS reads as input, assembles them de-novo, 
searches for evidence of each phylotype based on hidden Markov models (HMM), and then 
calculates their relative abundances. Should the genome be already assembled the contigs file can 
be indicated using the optional [-con] command line parameter.
The results and their respective visualisations are stored in a results directory and include a table 
with each phylotype and their relative frequencies as proportions of the full repertoire in the given 
genome; a heat map with dendrogram showing either absolute VAP variation or deviation from the 
mean, using our pilot dataset; and a Principal Component Analysis (PCA) plot showing variation 
distribution in the given sample compared to our pilot dataset containing isolates described by Silva Pereira et al. (2018) and Tihon et al. (2017).

Transcriptomic VAP
At the transcriptomic level the tool takes two NGS paired reads, maps the transcripts and estimates 
their abundance. It then searches for evidence of each phylotype based on hidden Markov models 
(HMM) and calculates the relative abundance per phylotype. The result is provided in both tabular 
form and a bar chart for comparison.

T.vivax:
The approach for T. vivax is quite different, it relies on the presence/absence of clusters of orthologs 
(COGs). It takes paired sequencing read files in fastq format and outputs a binary matrix of the 
presence/absence of each COG/gene for a given sample.
The results compare this matrix with a database of 27+ isolates; a heatmap and dendrogram are 
provided for comparison.


Instructions:
-------------

Vap.py 	- parses command line parameters and selects pathways accordingly 
	imports files 
		Tryp_G.py	- the T.congolense genomic pathway   	 
		Tryp_T.py	- the T.congolense transcriptomic pathway
		Tryp_V.py	- the T.vivax analysis via COG assessment
		Tryp_Multi.py	- manages multiple samples of the above three parthways
        
	Requires the Data directory for strain comparisons and geographical origins
		data/Motifs	- the hmm files for the 15 phylotypes searched for by hmmer
		data/Reference - fasta files for different strains used by T.congolense transcriptomic pathway
		data/vivax - Database and geo tags for T.vivax strains
		data/congodata.csv - the relative frequency of phylotypes appearing in T.congolense strains so far 
		data/congodata_deviationfromthemean.csv - as above but holding the deviation from the mean frequency 
		
	The python program Vap.py uses the following packages to analyze the isolates.
	Please ensure that these are installed and available to the python environment
	
	package				version used		website
	velvet				1.2.10				https://www.ebi.ac.uk/~zerbino/velvet/
	EMBOSS transeq		6.6.0.0         	http://emboss.open-bio.org/
	HMMER				3.1b2       		http://hmmer.org/
	bowtie2				2.2.6				http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
	samtools			1.6					http://www.htslib.org/
	cufflinks           2.2.1				http://cole-trapnell-lab.github.io/cufflinks/
	blast				2.7.1               https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
	
	As well as the usual python libraries Vap.py requires seaborn version 0.8.0 for clustermaps
	
	USAGE
	
	python Vap.py --help		lists the command line arguments accepted as below
	
	usage: Vap.py [-h] [-s S] [-con CON] [-t] [-p] [-strain STRAIN] [-dir DIR]
              [-cdir CDIR] [-f F] [-r R] [-k K] [-i I] [-cov COV]
              name

	Variant Antigen Profiler - the VAPPER.

	positional arguments:
	  name            Prefix for results directory and files therein

	optional arguments:
	  -h, --help      show this help message and exit
	  -s S            Species: T.congolense (default) or T.vivax
	  -con CON        Contigs File (fasta)
	  -t, -T          Transcriptomic Pathway
	  -p, -P          Export PDFs of images to results directory as well as .pngs 
	  -strain STRAIN  Strain for Transcriptomic pathway (defaults to Tc148)
	  -dir DIR        Directory that holds multiple paired NGS readfiles for analysis
	  -cdir CDIR      Directory that holds multiple pre-assembled contigs (fasta)files for analysis
	  -f F            Forward NGS read file
	  -r R            Reverse NGS Read File
	  -k K            kmers (default = 65) as used in velvet
	  -i I            Insert Length (default = 400) as used in velvet
	  -cov COV        Coverage cut off (default = 5) as used in velvet

	
	
	Example of use.
	
	T.congolense Genomic pathway:
	Single sample of T.congolense from paired NGS read files. 
	python Vap.py sgtest -f Test1.fastq -r Test2.fastq  
	Result images, csv files  and html file will be found in directory results/sgtest/
	
	Multiple sample of T.congolense from several sets of paired NGS read files place in directory /mydata/
	Each set of paired files should have the same name except for trailing 1 or 2 (eg Test1.fastq, Test2.fastq) 
	python Vap.py mgtest -dir mydata
	Result images, csv files  and html file will be found in directory results/mgtest/
	
	Single sample of T.congolense from a contigs file 
	python Vap.py sctest -con Test.fa 
	Result images, csv files  and html file will be found in directory results/sctest/
	
	Multiple sample of T.congolense from several contigs file (*.fa) placed in directory mycdata 
	python Vap.py mctest -cdir mycdata 
	Result images, csv files  and html file will be found in directory results/mctest/

	T.congolense Transcriptomic pathway
	
	Single sample of T.congolense, transcriptomic pathway from paired Transcript read files 
	python Vap.py sttest -t -f Transcripts.1 -r Transcripts.2 
	Result images, csv files  and html file will be found in directory results/sttest/
	
	Multiple sample of T.congolense from several sets of paired transcript read files place in directory /mytdata/
	Each set of paired files should have the same name except for trailing 1 or 2 (eg Transcripts.1, Transcripts.2) 
	python Vap.py mttest -t -dir mytdata
	Result images, csv files  and html file will be found in directory results/mttest/

	T.vivax: 
	Single sample of T.vivax from paired NGS read files. 
	python Vap.py svtest -s T.vivax -f Test1.fastq -r Test2.fastq  
	Result images, csv files  and html file will be found in directory results/svtest/
	
	Multiple sample of T.vivax from several sets of paired NGS read files place in directory /myvdata/
	Each set of paired files should have the same name except for trailing 1 or 2 (eg Test1.fastq, Test2.fastq) 
	python Vap.py mvtest -s T.vivax -dir myvdata
	Result images, csv files  and html file will be found in directory results/mvtest/
	
	Single sample of T.vivax from a contigs file 
	python Vap.py scvtest -s T.vivax -con Test.fa 
	Result images, csv files  and html file will be found in directory results/scvtest/
	
	Multiple sample of T.vivax from several contigs file (*.fa) placed in directory mycdata 
	python Vap.py mcvtest -s T.vivax -cdir mycdata 
	Result images, csv files  and html file will be found in directory results/mcvtest/
	
Examples:
-------------
	The directory "Example_data" contain examples of the outputs that should be expected. For  T. congolense, this includes 
	two PDF and PNG heatmaps/dendrograms; a PCA plot; and two CSV files containing the VAP of a test sample, expressed as the
	phylotype relative frequncy and variation (the deviation from the mean). For T. vivax, this includes a cluster map in the 
	form of heatmap/dendrogram, and a CSV file with a binary matrix representing a VAP of a test sample.
	
  
References:
-------------

Silva Pereira, S. et al. (2018) Variant antigen repertoires in Trypanosoma congolense populations and experimental infections can be profiled from deep sequence data with a set of universal protein motifs. 

Tihon, E. et al. (2017) Discovery and genomic analyses of hybridization between divergent lineages of Trypanosoma congolense, causative agent of Animal African Trypanosomiasis. Mol Ecol. 26(23):6524-6538. doi: 10.1111/mec.14271. Epub 2017 Aug 24.
