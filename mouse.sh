### the code was used to analyse the small RNA sequencing data from mouse for the in-trans pingpong analysis.

### software used
### fastx_toolkit
http://hannonlab.cshl.edu/fastx_toolkit/
### bowtie-1.2.3-linux-x86_64
https://github.com/BenLangmead/bowtie/releases/tag/v1.2.3
### samtools/1.10
https://github.com/samtools/samtools/releases/tag/1.10
### bedtools/2.28.0
https://github.com/arq5x/bedtools2/releases/tag/v2.28.0
### kent-ucsc tools
https://github.com/ucscGenomeBrowser/kent
### fastq_pair
https://github.com/linsalrob/fastq-pair
### STAR-2.7.3a
https://github.com/alexdobin/STAR/releases/tag/2.7.3a
### R/4.0.0
https://cran.r-project.org/src/base/R-4/R-4.0.0.tar.gz
### jellyfish 2.3.0
https://github.com/gmarcais/Jellyfish/releases/tag/v2.3.0

### define the variables
references="<directory where the genome sequence file and the gtf file are kept>"
raw_data="<directory where the fastq files are kept>"
analysis="<directory where the output files are kept>"
lib_sRNA="<name of the small RNA seq library>"
mkdir -p ${analysis}/${lib_sRNA}

### processing and analysing the small RNA sequencing libraries --- START ---
mkdir -p ${references}/mm10_GRCm39

### step 1: build bowtie index for miscRNA --- START ---
### fetch annotated sequences of the mouse infra-structural RNAs from ensembl
wget --directory-prefix="${references}/mm10_GRCm39" http://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz
wget --directory-prefix="${references}/mm10_GRCm39" http://ftp.ensembl.org/pub/release-106/gtf/mus_musculus/Mus_musculus.GRCm39.106.chr.gtf.gz
wget --directory-prefix="${references}/mm10_GRCm39" http://ftp.ensembl.org/pub/release-106/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz

zcat ${references}/mm10_GRCm39/Mus_musculus.GRCm39.ncrna.fa.gz |\
fasta_formatter - -t | awk '{split($5,a,":"); if($5 ~ "rRNA" || $5 ~ "miRNA" || $5 ~ "tRNA" ||  $5 ~ "snoRNA" || $5 ~ "snRNA") print $1"@"a[2],$NF
}' > ${references}/mm10_GRCm39/mouse_rRNA-tRNA-snoRNA-snRNA-miRNA.GRCm39.RebBase.tab

### Fetch RepBase entries for mouse repeats including rRNAs and tRNAs
fasta_formatter -i ${references}/mm10_GRCm39/mouse_ancestral_RepBase_2022-June.fasta -t |\
awk '{if($0 ~ "tRNA") print $1"@tRNA",$NF; else if($0 ~ "rRNA") print $1"@rRNA",$NF}' >> ${references}/mm10_GRCm39/mouse_rRNA-tRNA-snoRNA-snRNA-miRNA.GRCm39.RebBase.tab

### pick up rRNA, snRNA, snoRNA from the gtf file
### gtf files use 1-based coordinates
zcat ${references}/mm10_GRCm39/Mus_musculus.GRCm39.106.chr.gtf.gz |\
awk '{if($0~"rRNA" || $0~"snRNA" || $0~"snoRNA" || $0~"tRNA") print}' |\
awk -v ORS="" 'BEGIN{FS=";"} {print $1; for(i=1;i<=NF;i++) {if($i ~ "gene_biotype") print $i}; print "\n"}' | tr -d '"' |\
awk '{print $1,$4-1,$5,$10"@"$12,"1",$7}' | sort | uniq | tr ' ' '\t' > ${references}/mm10_GRCm39/Mus_musculus.GRCm39.106.rRNA-snRNA-snoRNA-Mt-tRNA.bed

### extract miscellaneous infra-structural RNA sequences from the major chromosomes and the mitochondrial genome
zcat ${references}/mm10_GRCm39/Mus_musculus.GRCm39.dna.toplevel.fa.gz |\
fasta_formatter - -t | awk '{if(NR<23) print ">"$1"\n"$NF}' > ${references}/mm10_GRCm39/Mus_musculus.GRCm39.nosmallchr.fasta
bedtools getfasta -tab -name -fi ${references}/mm10_GRCm39/Mus_musculus.GRCm39.nosmallchr.fasta \
-bed ${references}/mm10_GRCm39/Mus_musculus.GRCm39.106.rRNA-snRNA-snoRNA-Mt-tRNA.bed >> ${references}/mm10_GRCm39/mouse_rRNA-tRNA-snoRNA-snRNA-miRNA.GRCm39.RebBase.tab

### remove duplicates
awk '{print $1,toupper($2)}' ${references}/mm10_GRCm39/mouse_rRNA-tRNA-snoRNA-snRNA-miRNA.GRCm39.RebBase.tab |\
awk '!seen[$2]++' | awk '{print ">"$1"\n"$2}' > ${references}/mm10_GRCm39/mouse_rRNA-tRNA-snoRNA-snRNA-miRNA.GRCm39.RebBase.nonredundant.fasta

### make a bowtie index
mkdir -p ${references}/mm10_GRCm39/indices/
bowtie-build ${references}/mm10_GRCm39/mouse_rRNA-tRNA-snoRNA-snRNA-miRNA.GRCm39.RebBase.nonredundant.fasta ${references}/mm10_GRCm39/indices/mouse_rRNA-tRNA-snoRNA-snRNA-miRNA.GRCm39.RebBase.nonredundant
### step 1: build bowtie index for miscRNA --- END ---


### step 2: build bowtie index for the mouse genome --- START ---
### only consider major chromosomes, not including the mitochondrial genome
cat ${references}/mm10_GRCm39/Mus_musculus.GRCm39.nosmallchr.fasta | head -n 42 > ${references}/mm10_GRCm39/Mus_musculus.GRCm39.nosmallchr.noMT.fasta
bowtie-build ${references}/mm10_GRCm39/Mus_musculus.GRCm39.nosmallchr.noMT.fasta ${references}/mm10_GRCm39/indices/Mus_musculus.GRCm39.nosmallchr.noMT

### make a chromosome size file for later uses
fasta_formatter -i ${references}/mm10_GRCm39/Mus_musculus.GRCm39.nosmallchr.noMT.fasta -t |\
awk '{print $1,length($2)}' | tr ' ' '\t' > ${references}/mm10_GRCm39/Mus_musculus.GRCm39.nosmallchr.sizes
### step 2: build bowtie index for the mouse genome --- END ---


### step 3: Run RepeatMasker 4.1.0 on individual chromosomes by specifying "mouse" --- START ---
mkdir -p ${references}/mm10_GRCm39/RepeatMasker_4.1.0
fasta_formatter -i ${references}/mm10_GRCm39/Mus_musculus.GRCm39.nosmallchr.noMT.fasta -t |\
awk '{print ">"$1"\n"$2 > "${references}/mm10_GRCm39/RepeatMasker_4.1.0/Mus_musculus.GRCm39.nosmallchr.Chr"$1".fasta"}'

fastafile="${references}/mm10_GRCm39/RepeatMasker_4.1.0/Mus_musculus.GRCm39.nosmallchr.${CHR}.fasta"
RepeatMasker -species "mouse" -parallel 10 -xsmall $fastafile

### make a bed file from the RepeatMasker output
for CHR in Chr1 Chr2 Chr3 Chr4 Chr5 Chr6 Chr7 Chr8 Chr9 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17 Chr18 Chr19 ChrX ChrY; do
awk '{if($11 != "Simple_repeat" && $11 != "Low_complexity" && $11 != "rRNA" && $11 != "tRNA") print
}' ${references}/mm10_GRCm39/RepeatMasker_4.1.0/Mus_musculus.GRCm39.nosmallchr.${CHR}.fasta.out |\
awk '{if($9=="+") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9; else if($9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"
}' | tr ' ' '\t' >> ${references}/mm10_GRCm39/RepeatMasker_4.1.0/Mus_musculus.GRCm39.nosmallchr.allTE.bed
done
### step 3: Run RepeatMasker 4.1.0 on individual chromosomes by specifying "mouse" --- END ---


### step 4: generate a bowtie index for the RebBase TE sequences
### collected 429 mouse autonomous transposable elements (no simple repeats, no non-autonomous repeats) as of 2022-June
fasta_formatter -i ${references}/mm10_GRCm39/mouse_RepBase_transposable-elements_autonomous_2022-June.fasta -t |\
awk '{print ">"$1"\n"toupper($NF)}' > ${references}/mm10_GRCm39/mouse-autonomous-RepBase_429TEs.concise-name.fasta
bowtie-build ${references}/mm10_GRCm39/mouse-autonomous-RepBase_429TEs.concise-name.fasta ${references}/mm10_GRCm39/indices/mouse-autonomous-RepBase_429TEs.concise-name


### step 5: process the pachytene cells small RNA library from SRA sequence files from Basik, G&D, 2015 (PMID: 26115953) --- START ---
Illumina Genome Analyzer II sequencing; GSM1303504: het_RNF17_pachytene_cells_smallRNA; Mus musculus; miRNA-Seq	ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/003/SRR1104823/SRR1104823.fastq.gz

### step 5-0: download the fastq file
wget --directory-prefix="${raw_data}" ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/003/SRR1104823/SRR1104823.fastq.gz

### step 5-1: clip adapters for small RNA seq
lib_sRNA="SRR1104823"
mkdir -p ${analysis}/${lib_sRNA}
ADAPTER="CTGTAGGCACCATCAATC"
fastq_to_fasta -Q33 -i <(zcat ${raw_data}/${lib_sRNA}.fastq.gz) |\
# -c: discard non-clipped sequences, -l 18: discard sequences shorter than 18nt
fastx_clipper -c -a ${ADAPTER} -l 18 -i - | fasta_formatter - -t |\
awk '{if(length($NF)>17 && length($NF)<41) READS[$NF]++} END {
for(var in READS) print ">"var"@"READS[var]"\n"var}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa


### step 5-2: remove miscellaneous RNAs
misc_index="${references}/mm10_GRCm39/indices/mouse_rRNA-tRNA-snoRNA-snRNA-miRNA.GRCm39.RebBase.nonredundant"
bowtie -f -v 1 -a -S ${misc_index} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa \
--un ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-mapped.bed

### step 5-3: map reads to GRCm39 genome using bowtie, allowing up to 1MM, unique mappers
index_genome="${references}/mm10_GRCm39/indices/Mus_musculus.GRCm39.nosmallchr.noMT"
bowtie -f -v 1 -m 1 -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - | sort -k1,1 -k2,2n > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.mm10_GRCm39.unique.bed

### step 5-4: bowtie mapping to TEs, allowing up to 3MM, best and strata
TE_index="${references}/mm10_GRCm39/indices/mouse-autonomous-RepBase_429TEs.concise-name"
bowtie -f -v 3 --all --best --strata -S ${TE_index} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -b -S - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_TE_3MM_mappers.bed

### step 5-5: generate bigwig files of genome unique mappers --- START ---
### generate bedgraph files
TOTAL=`(awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]} END {print count}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.mm10_GRCm39.unique.bed)`
### bedgraph output per 1Mio miRNA mappers (plus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.mm10_GRCm39.unique.bed |\
bedtools genomecov -strand + -split -bg -i - -g ${references}/mm10_GRCm39/Mus_musculus.GRCm39.nosmallchr.sizes |\
awk -v TOTAL=${TOTAL} '{print $1,$2,$3,$4/TOTAL*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.mm10_GRCm39.unique.plus.bg
### bedgraph output per 1Mio miRNA mappers (minus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.mm10_GRCm39.unique.bed |\
bedtools genomecov -strand - -split -bg -i - -g ${references}/mm10_GRCm39/Mus_musculus.GRCm39.nosmallchr.sizes |\
awk -v TOTAL=${TOTAL} '{print $1,$2,$3,-$4/TOTAL*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.mm10_GRCm39.unique.minus.bg

### generate bigWig files
for strand in plus minus; do
bedGraphToBigWig ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.mm10_GRCm39.unique.${strand}.bg \
${references}/mm10_GRCm39/Mus_musculus.GRCm39.nosmallchr.sizes ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.mm10_GRCm39.unique.${strand}.bw
done
### step 5-5: generate bigwig files of genome unique mappers --- END ---


### step 5-6: map reads to GRCm39 genome using bowtie, allowing up to 1MM, all mappers, --all --best --strata option
index_genome="${references}/mm10_GRCm39/indices/Mus_musculus.GRCm39.nosmallchr.noMT"
bowtie -f -v 1 --all --best --strata -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.mm10_GRCm39.all.bed


### step 6: in-trans ping-pong analysis of genome unique and all piRNA mappers (>22nt)
mkdir -p ${analysis}/${lib_sRNA}/jellyfish
for CLASS in unique all; do

### step 6-1: extract g1g9, g2g10_revComp and last9 sequences from piRNAs
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.mm10_GRCm39.${CLASS}.bed |\
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr($4,1,9)}'  > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_g1g9.fasta
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.mm10_GRCm39.${CLASS}.bed |\
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr($4,2,9)}' | fastx_reverse_complement - > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_g2g10_revComp.fasta
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.mm10_GRCm39.${CLASS}.bed |\
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr(a[1],length(a[1])-8,9)}' > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_last9.fasta

### step 6-2: count the occurrences of 9mers by jellyfish
for TYPE in g2g10_revComp g1g9 last9; do
### jellyfish
jellyfish count -m 9 -s 100M -t 10 -o ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}_m9 ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}.fasta
jellyfish dump ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}_m9 > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}_m9_dump.fa
### collect all counts
TOTAL=`(fasta_formatter -i ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}_m9_dump.fa -t | awk '{count+=$1} END {print count}')`
### get rid of simple repeats
fasta_formatter -i ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_${TYPE}_m9_dump.fa -t |\
awk -v TYPE=${TYPE} -v TOTAL=${TOTAL} '{if($2!~"AAAAAA" && $2!~"CCCCCC" && $2!~"GGGGGG" && $2!~"TTTTTT") print $2,$1/TOTAL*1000,TYPE}' >> ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_all_m9_dump.tab
done

### step 6-3: measure linkage using R
Rscript jellyfish_linkage_Dspp_mouse.R DIRECTORY="${analysis}" LIB=${lib_sRNA} CLASS=${CLASS}
done

### step 6-4: collecting the in-trans pingpong linkage values across libraries
for CLASS in unique all; do
for lib_sRNA in w1118_ovary_oxidised_sRNA Dpse_ovary_oxidised_sRNA Deug_ovary_oxidised_sRNA_rep1 Deug_ovary_oxidised_sRNA_rep2 SRR1746887 SRR1104823; do
awk -v LIB=${lib_sRNA} -v CLASS=${CLASS} '{if(NR>1) print $1,$2,LIB"@"CLASS
}' ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_all_m9.linkage.txt >> ${analysis}/in-trans-pingpong.stats
done
done
### in-trans-pingpong_plots.R was used to make a barchart in Figure 6 and Extended Data Figure 9.
