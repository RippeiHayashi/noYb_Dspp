# noYb_Dspp
code for analysing High-throughput sequencing data for studying the piRNA pathway in a variety of Drosophila species.

### Drosophila_genome-assemblies_2022-Feb.txt
Figure 1 and Supplementary Table 1
356 genome assemblies from NCBI and 15 Nanopore genome assemblies from Miller 2018 (PMID: 30087105) that were used to search the yb homologs
<genome assembly name> <species name>

### annotating_Yb.sh
Figure 1 and Supplementary Table 1
The code was used to find the Yb homologs in the Drosophila genome assemblies.

### Repeatmasker-and-annotating_gypsy-env.sh
Supplementary Tables 2 and 3
The code was used to find env-carrying gypsy retrotransposon insertions in the Drosophila genome assemblies.

### Dmel.122TEs_Senti2015.fasta
122 transposon sequences from D. melanogaster used in this study.
It was originally annotated in Senti 2015 (PMID: 26302790)

### Dpse.non-redundant.top109.TEs.fasta
Top 109 TEs that are non-redundant and produce abundant piRNAs in D. pseudoobscura

### Deug.non-redundant.top94.TEs.fasta
Top 109 TEs that are non-redundant and produce abundant piRNAs in D. pseudoobscura

### Dmelanogaster.sh
The code was used to analyse the sequencing data for D. melanogaster.

### Dpseudoobscura.sh
The code was used to analyse the sequencing data for D. pseudoobscura.

### Deugracilis.sh
The code was used to analyse the sequencing data for D. eugracilis.

### Dbifasciata.sh
The code was used to analyse the sequencing data for D. bifasciata.

### Dmojavensis.sh
The code was used to analyse the sequencing data for D. mojavensis.

### mouse.sh
The code was used to analyse the sequencing data for mouse.

### in-trans-pingpong_simulation.sh
The code was used to perform the simulation for the in-trans pingpong analysis.

### plotting_TE-piRNAs.R
The code was used to plot 5' and 3' ends of TE mapping piRNA reads.
It was used for Figure 6 and Extended Data Figures 7 and 8.

### pingpong_linkage.R
The code was used to measure the linkage value of canonical 5' to 5' pingpong.
It was used for Figure 6.

### phasing_linkage.R
The code was used to measure the linkage value of 3' and 5' ends of TE antisense piRNAs during phasing.
It was used for Figure 6 and Extended Data Figure 8.

### jellyfish_linkage_Dspp_mouse.R
The code was used to measure the linkage value of in trans pingpong in the whole ovary small RNA libraries and the mouse pachytene stage small RNA library.

### jellyfish_linkage_IP.R
The code was used to measure the linkage value of in trans pingpong between the whole ovary small RNA library and Piwi and Aubergine IP small RNA libraries from D. eugracilis.

### tile_analysis_plots.R
The code was used to plot scatter plots of the piRNA coverage of 0.5kb genomic tiles.
It was used in Figure 5 and Extended Data Figure 6.

### TE-piRNAs_analysis_plots.R
The code was used to plot a variety of analyses on TE mapping piRNAs.
It was used in Figure 6 and Extended Data Figures 6 to 8.

### in-trans-pingpong_plots.R
The code was used to plot the results of the simulation of in trans pingpong in Extended Data Figure 9.
