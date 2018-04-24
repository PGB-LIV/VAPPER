# Trypanosoma-VAP
Variant Antigen Profiling for Trypanosoma congolense and Trypanosoma vivax

Introduction:
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
genome; a heat map with dendogram showing either absolute VAP variation or deviation from the 
mean, using our pilot dataset; and a Principal Component Analysis (PCA) plot showing variation 
distribution in the given sample compared to our pilot dataset.

Transcriptomic VAP
At the transcriptomic level the tool takes two NGS paired reads, maps the transcripts and estimates 
their abundance. It then searches for evidence of each phylotype based on hidden Markov models 
(HMM) and calculates the relative abundance per phylotype. The result is provided in both tabular 
form and a bar chart for comparison.

T.vivax:
The approach for T. vivax is quite different, it relies on the presence/absence of clusters of orthologs 
(COGs). It takes paired sequencing read files in fastq format and outputs a binary matrix of the 
presence/absence of each COG/gene for a given sample.
The results compare this matrix with a database of 20 isolates and a heatmap and dendrogram are 
provided for comparison.


System requirements

The python script require the following tools. 
Velvet for assembly		https://www.ebi.ac.uk/~zerbino/velvet/
EMBOSS transeq			http://emboss.open-bio.org/
HMMER for HMM search    http://hmmer.org/
bowtie2           		http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
samtools          		http://www.htslib.org/
cufflinks           	http://cole-trapnell-lab.github.io/cufflinks/
Blast     				https://blast.ncbi.nlm.nih.gov/Blast.cgi

Usage:

usage: Vap.py name [-h] [-s S] [-con CON] [-t] [-p] [-strain STRAIN] [-f F] [-r R]
              [-k K] [-i I] [-cov COV]
              name

Variant Antigen Profiler - the VAP.

positional arguments:
  name            Prefix for results directory/files

optional arguments:
  -h, --help      show this help message and exit
  -s S            Species: T.congolense (default) or T.vivax
  -con CON        Contigs File
  -t, -T          Transciptomic Pathway
  -p, -P          Export PDFs to results directory
  -strain STRAIN  strain required for Transcriptomic pathway
  -f F            Forward NGS read file
  -r R            Reverse NGS Read File
  -k K            kmers (default = 65)
  -i I            Insert Length (default = 400)
  -cov COV        Coverage cut off default = 5

  
  
