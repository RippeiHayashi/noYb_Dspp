### the code was used to analyse the sequencing data for Drosophila eugracilis.

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
### weblogo 3.7.8
https://github.com/WebLogo/weblogo

### define the variables
references="<directory where the genome sequence file and the gtf file are kept>"
raw_data="<directory where the fastq files are kept>"
analysis="<directory where the output files are kept>"
lib_sRNA="<name of the small RNA seq library>"
lib_pA="<name of the poly-A+ RNA seq library>"
lib_Q="<name of the Quant seq library>"
mkdir -p ${analysis}/${lib_sRNA}
mkdir -p ${analysis}/${lib_pA}
mkdir -p ${analysis}/${lib_Q}


### processing and analysing the small RNA sequencing libraries --- START ---

### below are the small RNA libraries processed in this part of the script:
Deug_ovary_oxidised_sRNA_rep1
Deug_ovary_oxidised_sRNA_rep2
Deug_embryo_oxidised_sRNA
Deug_Aub-IP_oxidised_sRNA
Deug_Piwi-IP_oxidised_sRNA
Deug_ovary_unoxidised_sRNA


### preparing for the analyses, do this before processing individual libraries --- START ---
### step 0-1: generating a bowtie index for miscellaneous RNAs
### step 0-1-1: download miscRNA sequences from the FlyBase
mkdir -p ${references}/dm6/indices
wget --directory-prefix="${references}/dm6" http://ftp.flybase.net/releases/FB2019_06/dmel_r6.31/fasta/dmel-all-miscRNA-r6.31.fasta.gz
wget --directory-prefix="${references}/dm6" http://ftp.flybase.net/releases/FB2019_06/dmel_r6.31/fasta/dmel-all-miRNA-r6.31.fasta.gz
wget --directory-prefix="${references}/dm6" http://ftp.flybase.net/releases/FB2019_06/dmel_r6.31/fasta/dmel-all-tRNA-r6.31.fasta.gz

### step 0-1-2: making fasta including miRNA, rRNA, snRNA, snoRNA, and tRNA sequences
zcat ${references}/dm6/dmel-all-miscRNA-r6.31.fasta.gz ${references}/dm6/dmel-all-miRNA-r6.31.fasta.gz ${references}/dm6/dmel-all-tRNA-r6.31.fasta.gz |\
fasta_formatter - -t | tr -d ";" |\
awk '{split($5,a,"="); if($2~"miRNA" || $2~"rRNA" || $2~"snRNA" || $2~"snoRNA" || $2~"tRNA") print ">"a[2]"\n"$NF
}' > ${references}/dm6/dmel-all-miRNA-rRNA-snRNA-snoRNA-tRNA-r6.31.fasta

### step 0-1-3: include spikes in the misc RNA and make a bowtie index
cat ${references}/dm6/dmel-all-miRNA-rRNA-snRNA-snoRNA-tRNA-r6.31.fasta ${references}/dm6/19mer.35mer.spikes.fa > ${references}/dm6/dmel-misc_spikes.fa
bowtie-build ${references}/dm6/dmel-misc_spikes.fa ${references}/dm6/indices/dmel-all-miRNA-rRNA-snRNA-snoRNA-tRNA-r6.31_spikes
### ${references}/dm6/19mer.35mer.spikes.fa looks as follows:
>19mer
CGTACGCGGGTTTAAACGA
>35mer
CTCATCTTGGTCGTACGCGGAATAGTTTAAACTGT


### step 0-2: generating a bowtie index for the genome
### step 0-2-1: download the genomic sequence of the GCF_018153835.1_ASM1815383v1 assembly of Drosophila eugracilis
mkdir -p ${references}/Deug/indices
wget --directory-prefix="${references}/Deug" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/153/835/GCF_018153835.1_ASM1815383v1/GCF_018153835.1_ASM1815383v1_genomic.fna.gz
gunzip ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna.gz

### step 0-2-2: make a bowtie index
bowtie-build ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna ${references}/Deug/indices/Deug_UStan

### step 0-2-3: generate the size file
fasta_formatter -i ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna -t |\
awk '{print $1,length($NF)}' | tr ' ' '\t' > ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.sizes


### step 0-3: generating a bed file for extended genomic regions of tRNAs
### step 0-3-1: download the the gtf file of the GCF_018153835.1_ASM1815383v1 genome assembly of Drosophila eugracilis
wget --directory-prefix="${references}/Deug" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/153/835/GCF_018153835.1_ASM1815383v1/GCF_018153835.1_ASM1815383v1_genomic.gtf.gz

### step 0-3-2: make a bed file of tRNA annotations extended by 100 bp both upstream and downstream
zcat ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.gtf.gz | grep 'gbkey "tRNA"' |\
awk '{print $1,$4-1-100,$5+100,$10,".",$7}' | awk '{if($2<0) $2=0; print}' |\
sort -k1,1 -k2,2n | uniq | tr -d '";"' | tr ' ' '\t' > ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.tRNA-extended.bed


### step 0-4: generating a bed tile for genome unique 0.5kb tiles from the GCF_018153835.1_ASM1815383v1 genome assembly of Drosophila eugracilis
### step 0-4-1: obtain every 25mer from the genome, only sense
### convert lowercase to uppercase first
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna |\
fasta_formatter - -t | awk '{for(i=0;i<=length($NF)-25;i++) print ">"$1"@"i"_+:"substr($NF,i+1,25)"\n"substr($NF,i+1,25)
}' > ${references}/Deug/Deug_UStan_25mers.plus.fasta

### step 0-4-2: obtain genome unique mappers
index_genome="${references}/Deug/indices/Deug_UStan"
bowtie -f -v 0 -m 1 -S ${index_genome} ${references}/Deug/Deug_UStan_25mers.plus.fasta |\
samtools view -@ 16 -bS - | bamToBed -i - > ${references}/Deug/Deug_UStan_25mers.plus.unique-mappers.bed

### step 0-4-3: count the unique mappers in 0.5kb bins
### take if the centre of 25mer sits between the coordinate of 1 to 500 as the first 0.5kb bin
### take as the second 0.5kb bin if it sits between 501 to 1000
awk '{split($4,a,"@"); split(a[2],b,"_"); TILE[$1"_"int((b[1]+12)/500)]++
} END {for(var in TILE) print var,TILE[var]}' ${references}/Deug/Deug_UStan_25mers.plus.unique-mappers.bed > ${references}/Deug/Deug_UStan_0.5kbtiles.counts

### step 0-4-4: make a bed file that contains all kb tiles that are gt 85% mappability
awk '{split($1,a,"_"); if($2>425) print a[1]"_"a[2],a[3]*500,(a[3]+1)*500,a[1]"_"a[2]":"a[3]*500"-"(a[3]+1)*500"_+ . +""\n"a[1]"_"a[2],a[3]*500,(a[3]+1)*500,a[1]"_"a[2]":"a[3]*500"-"(a[3]+1)*500"_- . -"
}' ${references}/Deug/Deug_UStan_0.5kbtiles.counts |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/Deug/Deug_UStan_85percent_0.5kbtiles.25mers.bed


### step 0-5: preparing for the 0.2kb genomic tile analyses of the all genome mappers --- START ---
### step 0-5-1: make a bed file of 0.2kb genomic tiles
awk '{for(i=0;i<=int($2/200)-1;i++) {print $1,i*200,(i+1)*200,$1":"i*200"-"(i+1)*200":- . -""\n"$1,i*200,(i+1)*200,$1":"i*200"-"(i+1)*200":+ . +"
} {print $1,int($2/200)*200,$2,$1":"int($2/200)*200"-"$2":- . -""\n"$1,int($2/200)*200,$2,$1":"int($2/200)*200"-"$2":+ . +"}
}' ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.sizes |\
tr ' ' '\t' > ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.200nt.window.bed

### step 0-5-2: add annotations of piRNA clusters, gypsy and TE insertions, and exons
mkdir -p ${analysis}/Deug
### step 0-5-2-1: piRNA clusters
bedtools intersect -wo -s -a ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.200nt.window.bed \
-b ${references}/Deug/Deug_UStan_piRNA-clusters.bed |\
awk '{print $4,$NF/200,$10}' >> ${analysis}/Deug/Deug_0.2kb_tiles.counts
### ${references}/Deug/Deug_UStan_piRNA-clusters.bed looks as follows:
NW_024573031.1  271670  282825  Cluster3031_uni .       +
NW_024572449.1  1       50026   Cluster2449_uni .       +
NW_024573703.1  1452    20698   Cluster3703_uni .       -
NW_024573378.1  18967   36000   Cluster3378_dual        .       -
NW_024573378.1  18967   36000   Cluster3378_dual        .       +

### step 0-5-2-2: Run RepeatMasker 4.1.0 by specifying "drosophila"
fastafile="${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna"
RepeatMasker -species "drosophila" -parallel 10 -xsmall ${fastafile}

### step 0-5-2-3: extract gypsy insertions --- START ---
cat ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna.out | grep -i "gypsy" |\
awk '{if($9=="+") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9; else if($9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"
}' | tr ' ' '\t' > ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna.out.gypsy.bed

### extract sense and antisense gypsy insertions, these bed files are used in the IGV browser tracks
cat ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna.out | grep -i "gypsy" |\
awk '{if($9=="+" ) print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9}' |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna.out.gypsy.plus.bed
cat ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna.out | grep -i "gypsy" |\
awk '{if($9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"}' |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna.out.gypsy.minus.bed

### annotate tiles that overlap with gypsy insertions at the sense orientation
bedtools intersect -f 0.5 -s -a ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.200nt.window.bed \
-b ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna.out.gypsy.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"gypsy_S"}' >> ${analysis}/Deug/Deug_0.2kb_tiles.counts

### annotate tiles that overlap with gypsy insertions at the antisense orientation
bedtools intersect -f 0.5 -S -a ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.200nt.window.bed \
-b ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna.out.gypsy.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"gypsy_AS"}' >> ${analysis}/Deug/Deug_0.2kb_tiles.counts
### step 0-5-2-3: extract gypsy insertions --- END ---

### step 0-5-2-4: extract TE insertions --- START ---
### remove "Simple_repeat", "Low_complexity", "rRNA" and "tRNA"
awk '{if($11 != "Simple_repeat" && $11 != "Low_complexity" && $11 != "rRNA" && $11 != "tRNA") print}' ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna.out |\
awk '{if($9=="+") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9; else if($9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"
}' | tr ' ' '\t' > ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna.out.allTE.bed

### annotate tiles that overlap with TE insertions at the sense orientation
bedtools intersect -f 0.5 -s -a ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.200nt.window.bed \
-b ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna.out.allTE.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"TE_S"}' >> ${analysis}/Deug/Deug_0.2kb_tiles.counts

### annotate tiles that overlap with TE insertions at the antisense orientation
bedtools intersect -f 0.5 -S -a ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.200nt.window.bed \
-b ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna.out.allTE.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"TE_AS"}' >> ${analysis}/Deug/Deug_0.2kb_tiles.counts
### step 0-5-2-4: extract TE insertions --- END ---

### step 0-5-2-5: add exon annotations
zcat ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.gtf.gz | tr -d '";' |\
awk 'BEGIN{OFS="\t"} {if(NR>4 && $3=="exon") print $1,$4,$5,$10":"$12,".",$7}' | sort -k1,1 -k2,2n > ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.exon.bed
bedtools intersect -wo -s -a ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.200nt.window.bed \
-b ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.exon.bed |\
awk '{exon[$4]+=$NF} END {for(var in exon) print var,exon[var]/200,"exon"}' >> ${analysis}/Deug/Deug_0.2kb_tiles.counts
### step 0-5: preparing for the 0.2kb genomic tile analyses of the all genome mappers --- END ---


### step 0-6: make a bowtie index for all Drosophila autonomous transposons downloaded from the RepBase in March 2022 (3053 entries)
mkdir -p ${references}/TEs/indices
fasta_formatter -i ${references}/TEs/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.fasta -t |\
awk '{print ">"$1"\n"$NF}' > ${references}/TEs/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.concise.fasta
bowtie-build ${references}/TEs/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.concise.fasta ${references}/TEs/indices/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.concise


### step 0-7: obtain TEs that are non-redundant and produce abundant piRNAs --- START ---
### We used the whole ovary oxidised small RNA library (Deug_ovary_oxidised_sRNA_rep1) to obtain 150 TEs that produce most abundant piRNAs.
lib_sRNA="Deug_ovary_oxidised_sRNA_rep1"
### step 0-7-1: take all reads that mapped to the D.eugracilis genome, and use them to map against transposon sequences
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); print ">"$4"\n"a[1]}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.fa

### step 0-7-2: bowtie mapping to TEs, --all --best --strata option, allowing up to three mismatches
TE_index="${references}/TEs/indices/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.concise"
bowtie -f -v 3 --all --best --strata -S ${TE_index} ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_Deug_UStan_mapped_TE_3MM_mappers.bed

### step 0-7-3: We manually removed 56 TEs that resemble other TEs within the top 150 TEs, which left 94 non-redundant TEs.

### step 0-7-4: check the redundancy between the remaining 94 TEs by counting the piRNA reads that mapped to multiple TE sequences.
mkdir -p ${analysis}/${lib_sRNA}/TE_mapping.94TEs

### step 0-7-4-1: print out all the piRNA reads (gt22nt) mapping to individual TEs
cat ${references}/Deug/Deug.non-redundant.top94.TEs | while read TE READS; do
awk -v TE=${TE} '{if($1==TE && $3-$2>22) print $4}' ${analysis}/${lib_sRNA}/${lib_sRNA}_Deug_UStan_mapped_TE_3MM_mappers.bed |\
sort | uniq > ${analysis}/${lib_sRNA}/TE_mapping.94TEs/${lib_sRNA}_${TE}.mappers.nuc
done

### step 0-7-4-2: collect all TEs mapping piRNA reads that mapped to multiple TEs
cat ${analysis}/${lib_sRNA}/TE_mapping.94TEs/*.mappers.nuc |\
awk '{READS[$1]++} END {for(var in READS) {if(READS[var]>1) print var}}' |\
sort > ${analysis}/${lib_sRNA}/TE_mapping.94TEs/${lib_sRNA}_top94.redundant-mappers.nuc

### step 0-7-4-3: per each TE, measure how many of piRNA mappers mapped to other TEs
cat ${references}/Deug/Deug.non-redundant.top94.TEs | while read TE READS; do
printf $TE"\n""total reads count gt22nt: "$READS"\n" >> ${analysis}/${lib_sRNA}/TE_mapping.94TEs/${lib_sRNA}_top94.stats
join -1 1 -2 1 ${analysis}/${lib_sRNA}/TE_mapping.94TEs/${lib_sRNA}_${TE}.mappers.nuc ${analysis}/${lib_sRNA}/TE_mapping.94TEs/${lib_sRNA}_top94.redundant-mappers.nuc |\
awk '{split($1,a,"@"); count+=a[2]} END {print "redundant reads count gt22nt: "count}' >> ${analysis}/${lib_sRNA}/TE_mapping.94TEs/${lib_sRNA}_top94.stats
done
### This analysis confirmed that the non-redundant 94 TEs have more than 10% of redundant piRNA mappers.
### step 0-7: obtain TEs that are non-redundant and produce abundant piRNAs --- END ---


### step 0-8: generate pseudocounts of 1 at every coordinate of non-redundant TEs
mkdir -p ${references}/TEs/TE_pseudo.piRNAs
### run it per each TE
cat ${references}/Deug/Deug.non-redundant.top94.TEs | while read TE READS; do
### generate pseudocounts of 1 at every coordinate
length=`(fasta_formatter -i ${references}/TEs/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.fasta -t | awk -v TE=${TE} '{if($1==TE) print length($NF)}')`
if [[ ! -f ${references}/TEs/TE_pseudo.piRNAs/${TE}_pseudo.piRNAs.txt ]]; then
### generate pseudocounts of 1 at every coordinate
for ((i=0; i<=${length}-1; i++)); do
j=$((${i}+1))
printf "${TE} ${i} ${j} AAAAAAAAAAAAAAAAAAAAAAA@1 . +""\n""${TE} ${i} ${j} AAAAAAAAAAAAAAAAAAAAAAA@1 . -""\n" |\
tr ' ' '\t' >> ${references}/TEs/TE_pseudo.piRNAs/${TE}_pseudo.piRNAs.txt
done
fi

### preparing for the analyses, do this before processing individual libraries --- END ---


### common analyses for all small RNA libraries --- START ---

### step 1: trim the adapter sequence and collapse reads by sequence
### we take R1 reads from the paired end sequencing library.
fastq_to_fasta -Q33 -i <(zcat ${raw_data}/${lib_sRNA}.fastq.gz) |\
# -c: discard non-clipped sequences, -l 18: discard sequences shorter than 18nt
# collapse reads of the same sequence to one entry while retaining the number of reads
fastx_clipper -c -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -l 18 -i - | fasta_formatter - -t |\
awk '{if(length($NF)>25 && length($NF)<49) READS[substr($NF,5,length($NF)-8)]++} END {
for(var in READS) print ">"var"@"READS[var]"\n"var}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa
### ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa looks as follows:
>TCCAGATGAACGGTTAAGTGTCCAAAAAG@12
TCCAGATGAACGGTTAAGTGTCCAAAAAG
>TACTTGAAAGAATCAGGGGCCAACCAT@5
TACTTGAAAGAATCAGGGGCCAACCAT
...


### step 2: map to the Drosophila melanogaster miscellaneous RNA (misc RNA) allowing up to 1MM
### run bowtie to map reads to the miscRNA and take unmapped reads
misc_index="${references}/dm6/indices/dmel-all-miRNA-rRNA-snRNA-snoRNA-tRNA-r6.31_spikes"
bowtie -f -v 1 -a -S ${misc_index} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa \
--un ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-mapped.bed


### step 3: obtain genome unique mappers
### run bowtie to map reads to the genome, allowing up to 1MM, unique mappers
index_genome="${references}/Deug/indices/Deug_UStan"
bowtie -f -v 1 -m 1 -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - | sort -k1,1 -k2,2n > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Deug_UStan.bed


### step 4: remove tRNA reads by intersection
### exclude reads that intersect with tRNA annotations from the genome unique mappers
bedtools intersect -v -s -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Deug_UStan.bed \
-b ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.tRNA-extended.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Deug_UStan.tRNA-excluded.bed


### step 5: generate the bigwig files for the genome unique mappers --- START ---
### step 5-1: make bedgraph files
### only count reads that are greater than 22nt
### for IP libraries and oxidised libraries --- START ---
### normalise read counts to one million genome unique mappers
TOTAL=`(awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]} END {print count
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Deug_UStan.tRNA-excluded.bed)`
### bedgraph output per 1Mio miRNA mappers (plus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Deug_UStan.tRNA-excluded.bed |\
bedtools genomecov -strand + -split -bg -i - -g ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.sizes |\
awk -v TOTAL=${TOTAL} '{print $1,$2,$3,$4/TOTAL*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_plus.Deug_UStan.tRNA-excluded.bg
### bedgraph output per 1Mio miRNA mappers (minus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Deug_UStan.tRNA-excluded.bed |\
bedtools genomecov -strand - -split -bg -i - -g ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.sizes |\
awk -v TOTAL=${TOTAL} '{print $1,$2,$3,-$4/TOTAL*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_minus.Deug_UStan.tRNA-excluded.bg
### for IP libraries and oxidised libraries --- END ---

### for non-oxidised libraries --- START ---
### normalise read counts to one million microRNA reads
MIRNA=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-mapped.bed | grep mir | tr '@' ' ' | awk '{count+=$5} END {print count}')`
### bedgraph output per 1Mio miRNA mappers (plus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Deug_UStan.tRNA-excluded.bed |\
bedtools genomecov -strand + -split -bg -i - -g ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.sizes |\
awk -v MIRNA=${MIRNA} '{print $1,$2,$3,$4/MIRNA*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_plus.Deug_UStan.tRNA-excluded.bg
### bedgraph output per 1Mio miRNA mappers (minus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Deug_UStan.tRNA-excluded.bed |\
bedtools genomecov -strand - -split -bg -i - -g ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.sizes |\
awk -v MIRNA=${MIRNA} '{print $1,$2,$3,-$4/MIRNA*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_minus.Deug_UStan.tRNA-excluded.bg
### for non-oxidised libraries --- END ---

### step 5-2: make bigwig files
for strand in plus minus; do
bedGraphToBigWig ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_${strand}.Deug_UStan.tRNA-excluded.bg \
${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.sizes ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_${strand}.Deug_UStan.tRNA-excluded.bw
done
### step 5: generate the bigwig files for the genome unique mappers --- END ---


### step 6: 0.5kb tile analysis for the genome unique mappers
### count piRNA reads (greater than 22nt) that uniquely mapped to the 0.5kb tiles
READS=`(awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]} END {print count
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Deug_UStan.tRNA-excluded.bed)`
bedtools intersect -wo -s -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Deug_UStan.tRNA-excluded.bed \
-b ${references}/Deug/Deug_UStan_85percent_0.5kbtiles.25mers.bed |\
awk -v READS=${READS} -v LIB=${lib_sRNA} '{split($4,a,"@"); if(length(a[1])>22) TILE[$10]+=a[2]} END {for(var in TILE) print var,TILE[var]/READS*1000000,LIB
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Deug_UStan.tRNA-excluded.gt22nt.0.5kbtiles.counts


### step 7: 0.2kb tile analysis for the genome all mappers --- START ---
### step 7-1: map reads to GCF_018153835.1_ASM1815383v1 genome assembly of Drosophila eugracilis, allowing up to 1MM, all best strata
index_genome="${references}/Deug/indices/Deug_UStan"
bowtie -f -v 1 --all --best --strata -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.bed

### step 7-2: add mapping counts to the fasta entries. Mapping counts will be used to distribute the read coverage across tiles.
awk '{READS[$4]++} END {for(var in READS) {split(var,a,"@"); print ">"var":"READS[var]"\n"a[1]}
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.mapping-counts.fa

### step 7-3: map again all best strata
index_genome="${references}/Deug/indices/Deug_UStan"
bowtie -f -v 1 --all --best --strata -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.mapping-counts.fa |\
samtools view -bS - | bamToBed -i - | sort -k1,1 -k2,2n > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.bed

### step 7-4: exclude reads that intersect with tRNA annotations
bedtools intersect -v -s -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.bed \
-b ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.tRNA-extended.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.bed

### step 7-5: collect read counts per 200nt windows from all-best-strata
TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if($3-$2>22) count+=b[1]} END {print count}')`
awk 'BEGIN{OFS="\t"} {if($3-$2 > 22) print}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.bed |\
bedtools intersect -s -wo -F 0.5 -a ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.200nt.window.bed -b - |\
awk -v TOTAL=${TOTAL} -v LIB=${lib_sRNA} '{split($10,a,"@"); split(a[2],b,":"); TILE[$4]+=b[1]/b[2]} END {for(var in TILE) print var,TILE[var]/TOTAL*1000000,LIB"_abs"
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.gt22nt.0.2kbtiles.counts
### step 7: 0.2kb tile analysis for the genome all mappers --- END ---


### step 8: analyse TE mapping piRNAs --- START ---
### step 8-1: take all reads that mapped to the D.eugracilis genome, and use them to map against transposon sequences
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); print ">"$4"\n"a[1]}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.fa

### step 8-2: bowtie mapping to TEs, --all --best --strata option, allowing up to three mismatches
TE_index="${references}/TEs/indices/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.concise"
bowtie -f -v 3 --all --best --strata -S ${TE_index} ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_Deug_UStan_mapped_TE_3MM_mappers.bed

### step 8-3: making TE mapping plots and measuring pingpong and phasing signatures --- START ---
mkdir -p ${analysis}/${lib_sRNA}/End.plots.TEs.linkage
### step 8-3-1: normalise the counts of 3' end 5' ends per 1 million total genome mappers (>22nt)
TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if(length(a[1])>22) count+=b[1]} END {print count}')`

### step 8-3-2: count 5end and 3end of plus and minus mapping reads (>22nt)
### run it per each TE
cat ${references}/Deug/Deug.non-redundant.top94.TEs | while read TE READS; do
### use pseudocounts created in step 0-8
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_Deug_UStan_mapped_TE_3MM_mappers.bed ${references}/TEs/TE_pseudo.piRNAs/${TE}_pseudo.piRNAs.txt |\
awk -v TE=${TE} -v TOTAL=${TOTAL} '{split($4,a,"@"); if($6=="+" && length(a[1])>22 && $1==TE) {PLUS5[$2]+=a[2]; PLUS3[$3-1]+=a[2]
} else if($6=="-" && length(a[1])>22 && $1==TE) {MINUS5[$3-1]+=a[2]; MINUS3[$2]+=a[2]
}} END {for(var in PLUS5) print var,(PLUS5[var]-1)/TOTAL*1000000,"plus_5end""\n"var,(PLUS3[var]-1)/TOTAL*1000000,"plus_3end""\n"var,(MINUS5[var]-1)/TOTAL*1000000,"minus_5end""\n"var,(MINUS3[var]-1)/TOTAL*1000000,"minus_3end"
}' > ${analysis}/${lib_sRNA}/End.plots.TEs.linkage/${lib_sRNA}_${TE}_3MM_mappers.counts
### use R scripts to make plots and measure linkages
DIRECTORY="${analysis}/${lib_sRNA}/End.plots.TEs.linkage"
Rscript plotting_TE-piRNAs.R DIRECTORY=${DIRECTORY} LIB=${lib_sRNA} TE=${TE} LENGTH=${length}
Rscript pingpong_linkage.R DIRECTORY=${DIRECTORY} LIB=${lib_sRNA} TE=${TE}
Rscript phasing_linkage.R DIRECTORY=${DIRECTORY} LIB=${lib} TE=${TE}
done
### step 8-3: making TE mapping plots and measuring pingpong and phasing signatures --- END ---

### step 8-4: examine the nucleotide composition around the 3' ends of TE mapping piRNAs using weblogo
mkdir -p ${analysis}/${lib_sRNA}/weblogo
TE_fasta="${references}/TEs/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.concise.fasta"
### run it per each TE
cat ${references}/Deug/Deug.non-redundant.top94.TEs | while read TE READS; do
LENGTH=`(cat ${TE_fasta} | fasta_formatter - -t | awk -v TE=${TE} '{if($1==TE) print length($NF)}')`
### extract sequences of the eleven nucleotides around the piRNA 3' ends
awk -v TE=${TE} '{if($1==TE && $3-$2>22) print}' ${analysis}/${lib_sRNA}/${lib_sRNA}_Deug_UStan_mapped_TE_3MM_mappers.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); for(i=1;i<=b[1];i++) {if($6=="+") print $1,$3-5,$3+6,$4,$5,$6; else if($6=="-" && $1==TE && $3-$2>22) print $1,$2-6,$2+5,$4,$5,$6}}' |\
awk -v LENGTH=${LENGTH} '{if($2>=0 && $3<=LENGTH) print}' | tr ' ' '\t' |\
bedtools getfasta -s -fi ${TE_fasta} -tab -bed - | awk '{print ">"$1"\n"toupper($2)}' | tr 'T' 'U' > ${analysis}/${lib_sRNA}/weblogo/${lib_sRNA}_${TE}_3end_11nt_window.fasta
### generate the weblogo
weblogo -U probability -A rna -f ${analysis}/${lib_sRNA}/weblogo/${lib_sRNA}_${TE}_3end_11nt_window.fasta -F pdf -n 50 -c classic -o ${analysis}/${lib_sRNA}/weblogo/${lib_sRNA}_${TE}_3end_11nt_window_gt22_logo_prob.pdf
done

### step 8-5: size profile and nucleotide composition of TE mapping reads
TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_Deug_UStan_mapped_TE_3MM_mappers.bed | awk '{split($4,a,"@"); split(a[2],b,":"); count+=b[1]} END {print count}')`
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_Deug_UStan_mapped_TE_3MM_mappers.bed |\
awk -v TOTAL=${TOTAL} -v LIB=${lib_sRNA} '{split($4,a,"@"); if($6=="-") SIZE[$3-$2" "substr(a[1],1,1)]+=a[2]/TOTAL*1000} END {for(var in SIZE) print var,SIZE[var],LIB
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_1st-nucleotide_length.TE_AS.table
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_Deug_UStan_mapped_TE_3MM_mappers.bed |\
awk -v TOTAL=${TOTAL} -v LIB=${lib_sRNA} '{split($4,a,"@"); if($6=="+") SIZE[$3-$2" "substr(a[1],1,1)]+=a[2]/TOTAL*1000} END {for(var in SIZE) print var,SIZE[var],LIB
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_1st-nucleotide_length.TE_S.table
### TE-piRNAs_analysis_plots.R was used to make a plot used in Extended Data Figure 8.
### step 8: analyse TE mapping piRNAs --- END ---


### step 9: measure in-trans ping-pong linkage --- START ---
mkdir -p ${analysis}/${lib_sRNA}/jellyfish

### step 9-1: extract g1g9, g2g10_revComp and last9 sequences from piRNA genome unique mappers
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr($4,1,9)
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Deug_UStan.tRNA-excluded.bed > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_unique_g1g9.fasta
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr($4,2,9)
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Deug_UStan.tRNA-excluded.bed | fastx_reverse_complement - > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_unique_g2g10_revComp.fasta
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr(a[1],length(a[1])-8,9)
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Deug_UStan.tRNA-excluded.bed > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_unique_last9.fasta

### step 9-2: extract g1g9, g2g10_revComp and last9 sequences from piRNA genome all mappers
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if($3-$2>22) for(i=1;i<=b[1];i++) print ">"$4":"i"\n"substr($4,1,9)}'  > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_all_g1g9.fasta
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if($3-$2>22) for(i=1;i<=b[1];i++) print ">"$4":"i"\n"substr($4,2,9)}' | fastx_reverse_complement - > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_all_g2g10_revComp.fasta
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if($3-$2>22) for(i=1;i<=b[1];i++) print ">"$4":"i"\n"substr(a[1],length(a[1])-8,9)}' > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_all_last9.fasta

### step 9-3: count the occurrences of 9mers by jellyfish
for CLASS in unique all; do
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
done

### step 9-4: measure linkage using R
for CLASS in unique all; do
Rscript jellyfish_linkage_Dspp_mouse.R DIRECTORY="${analysis}" LIB=${lib_sRNA} CLASS=${CLASS}
done
### step 9: measure in-trans ping-pong linkage --- END ---

### common analyses for all small RNA libraries --- END ---


### post-processing of the analysed data for small RNA libraries --- START ---

### step 10: 0.5kb tile analysis of the genome unique mappers
### step 10-1: add annotations of piRNA cluster tiles to the tiles
mkdir -p ${analysis}/Deug
bedtools intersect -wo -s -a ${references}/Deug/Deug_UStan_85percent_0.5kbtiles.25mers.bed \
-b ${references}/Deug/Deug_UStan_piRNA-clusters.bed |\
awk '{print $4,$NF/500,$10}' > ${analysis}/Deug/deug_0.5kb_tiles.counts

### step 10-2: add counts from small RNA libraries
for lib_sRNA in <oxidied whole ovary small RNA library> <oxidied embryonic small RNA library> <oxidied Aubergine IP library> <oxidied Piwi IP library>; do
#for lib_sRNA in Deug_ovary_oxidised_sRNA_rep1 Deug_embryo_oxidised_sRNA Deug_Aub-IP_oxidised_sRNA Deug_Piwi-IP_oxidised_sRNA; do
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Deug_UStan.tRNA-excluded.gt22nt.0.5kbtiles.counts >> ${analysis}/Deug/deug_0.5kb_tiles.counts
done
### step 10-3: tile_analysis_plots.R was used to make a scatter plot in Figure 5.


### step 11: 0.2kb tile analysis of the genome all mappers --- START ---
### step 11-1: add counts from small RNA libraries
for lib_sRNA in <oxidied whole ovary small RNA library> <oxidied embryonic small RNA library>; do
#for lib_sRNA in Deug_ovary_oxidised_sRNA_rep1 Deug_embryo_oxidised_sRNA; do
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.gt22nt.0.2kbtiles.counts >> ${analysis}/Deug/Deug_0.2kb_tiles.counts
done

### step 11-2: make a table using reshape2 in R and count piRNA reads from individual annotations
#R/4.0.0
library(reshape2)
table <- read.table("${analysis}/Deug/Deug_0.2kb_tiles.counts", header=F)
table_d=dcast(table,V1~V3,value.var="V2")
table_d[is.na(table_d)] <- 0
write.table(table_d, file="${analysis}/Deug/Deug_0.2kb_tiles.counts.table", quote=F, col.names=T, row.names = F)
colnames(table_d)

[1] "V1"                            "Cluster2449_uni"
[3] "Cluster3031_uni"               "Cluster3378_dual"
[5] "Cluster3703_uni"               "Deug_embryo_oxidised_sRNA_abs"
[7] "Deug_ovary_oxidised_sRNA_rep1_abs" "exon"
[9] "gypsy_AS"                      "gypsy_S"
[11] "TE_AS"                         "TE_S"
#R/4.0.0

### count total somatic reads
soma_total=`(awk '{if(NR>1 && $7>0) print}' ${analysis}/Deug/Deug_0.2kb_tiles.counts.table |\
awk '{if($6/$7 < 0.1) print}' | awk '{count+=$7} END {print count}')`

### count somatic uni-stranded cluster reads
uni_total=`(awk '{if(NR>1 && $7>0) print}' ${analysis}/Deug/Deug_0.2kb_tiles.counts.table |\
awk '{if($6/$7 < 0.1) print}' | awk '{if($2>0 || $3>0 || $5>0) count+=$7} END {print count}')`

### count exonic piRNA reads
### exclude tiles that overlap with uni clusters and TE annotations
exons=`(awk '{if(NR>1 && $7>0) print}' ${analysis}/Deug/Deug_0.2kb_tiles.counts.table |\
awk '{if($6/$7 < 0.1) print}' | awk '{if($2==0 && $3==0 && $5==0 && $11==0 && $12==0 && $8>0) count+=$7} END {print count}')`

### count somatic gypsy_AS reads outside the clusters
no_uni_gypsy_AS=`(awk '{if(NR>1 && $7>0) print}' ${analysis}/Deug/Deug_0.2kb_tiles.counts.table |\
awk '{if($6/$7 < 0.1) print}' | awk '{if($2==0 && $3==0 && $5==0 && $9>0 && $10==0) count+=$7} END {print count}')`

### count somatic gypsy_S reads outside the clusters
no_uni_gypsy_S=`(awk '{if(NR>1 && $7>0) print}' ${analysis}/Deug/Deug_0.2kb_tiles.counts.table |\
awk '{if($6/$7 < 0.1) print}' | awk '{if($2==0 && $3==0 && $5==0 && $10>0 && $9==0) count+=$7} END {print count}')`

### collect the stats
others=`echo "$soma_total - $uni_total - $exons - $no_uni_gypsy_AS - $no_uni_gypsy_S" | bc`
printf $soma_total" Deug soma_total""\n"$uni_total" Deug uni_total""\n"$exons" Deug exons""\n"$no_uni_gypsy_AS" Deug no_uni_gypsy_AS""\n"$no_uni_gypsy_S" Deug no_uni_gypsy_S""\n"$others" Deug others""\n" >> ${analysis}/0.2kb_tiles.counts.stats
### step 11: 0.2kb tile analysis of the genome all mappers --- END ---


### step 12: comparing the abundance of somatic and germline TE mapping piRNA reads (>22nt) --- START ---
for lib_sRNA in <oxidied whole ovary small RNA library> <oxidied embryonic small RNA library>; do
#for lib_sRNA in Deug_ovary_oxidised_sRNA_rep1 Deug_embryo_oxidised_sRNA; do
### normalise the counts to 1 million total genome mappers (>22nt)
TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if(length(a[1])>22) count+=b[1]} END {print count}')`

### step 12-1: count piRNA mappers (>22nt) both sense and antisense
cat ${references}/Deug/Deug.non-redundant.top94.TEs | while read TE READS; do
## some TE sequences have repeated sequences within and the same reads can map multiple times
awk '!seen[$1" "$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_Deug_UStan_mapped_TE_3MM_mappers.bed |\
awk -v TE=${TE} -v TOTAL=${TOTAL} -v LIB=${lib_sRNA} '{split($4,a,"@"); split(a[2],b,":"); if($3-$2>22 && $1==TE) count+=b[1]} END {print TE,count/TOTAL*1000000,LIB
}' >> ${analysis}/${lib_sRNA}/${lib_sRNA}_94TE_mappers.counts
done

### step 12-2: combine counts from the whole ovary and embryonic small RNA libraries
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_94TE_mappers.counts >> ${analysis}/Deug_94TE_mappers.counts
done

### step 12-3: TE-piRNAs_analysis_plots.R was used to make a scatter plot in Extended Data Figure 6.

### step 12-4: identify TEs that are enriched more than three times in the ovary sRNAseq library than in the embryonic sRNAseq library
#R/4.0.0
library(reshape2)
table <- read.table("${analysis}/Deug_94TE_mappers.counts", header=F)
table.d <- dcast(table,V1~V3,value.var="V2")
write.table(table.d, file="${analysis}/Deug_94TE_mappers.counts.table", quote=F, row.names=F)
#R/4.0.0

awk '{if(NR>1 && $3>$2*3) print $1}' ${analysis}/Deug_94TE_mappers.counts.table
Gypsy-1_DBi-I
Gypsy-28_DBp-I
Gypsy-35_DEl-I
Gypsy-6_DEu-I
### step 12: comparing the abundance of somatic and germline TE mapping piRNA reads (>22nt) --- END ---


### step 13: measuring the strand bias of the piRNA reads (>22nt) mapping to non-redundant 94 TEs
### exclude somatic TEs as well as Mariner DNA transposons. Mariner transposons have terminal inverted repeats (TIR). It is not possible to determine the strandness of piRNAs mapping to the TIRs.
### we used the oxidised whole ovary small RNA library, replicate 1
lib_sRNA="<oxidied whole ovary small RNA library>"
lib_sRNA="Deug_ovary_oxidised_sRNA_rep1"
awk '{if($1 != "Gypsy-1_DBi-I" && $1 != "Gypsy-28_DBp-I" && $1 != "Gypsy-35_DEl-I" && $1 != "Gypsy-6_DEu-I" && $1 != "Mariner-19_DK" && $1 != "Mariner-19_DAn" && $1 != "Mariner-2_Dro" && $1 != "Mariner-3_DAn" && $1 != "Mariner-5_DEu" && $1 != "Mariner-18_DAn") print
}' ${references}/Deug/Deug.non-redundant.top94.TEs | while read TE READS; do
### we consider 84 TEs in total for D. eugracilis
awk '!seen[$1" "$4" "$6]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_Deug_UStan_mapped_TE_3MM_mappers.bed |\
awk -v TE=${TE} -v LIB=${lib_sRNA} '{split($4,a,"@"); split(a[2],b,":"); {if($3-$2>22 && $1==TE && $6=="+") count_S+=b[1]; else if($3-$2>22 && $1==TE && $6=="-") count_AS+=b[1]
}} END {print TE,count_AS/(count_S+count_AS),LIB}' >> ${analysis}/TE_S_AS_bias.stats
done
## TE-piRNAs_analysis_plots.R was used to make a box plot in Figure 6.


### step 14: collecting the canonical pingpong linkage values of TE mapping piRNA reads (>22nt)
### we only considered TEs that produce both sense and antisense piRNAs by at least 10% of the total TE piRNAs
awk '{if($2>0.1 && $2<0.9) print}' ${analysis}/TE_S_AS_bias.stats | while read TE AS_ratio lib_sRNA; do
awk -v TE=${TE} -v LIB=${lib_sRNA} '{if(NR==11) print TE,$1,LIB}' ${analysis}/${lib_sRNA}/End.plots.TEs.linkage/${lib_sRNA}_${TE}_3MM_mappers.pingpong.plus_5end.minus_5end.table >> ${analysis}/pingpong.stats.selected
done
### TE-piRNAs_analysis_plots.R was used to make a box plot in Figure 6.


### step 15: collecting the phasing (minus_3end.minus_5end) linkage values of TE mapping piRNA reads (>22nt)
for lib_sRNA in Deug_ovary_oxidised_sRNA_rep1 Deug_Aub-IP_oxidised_sRNA Deug_Piwi-IP_oxidised_sRNA; do
awk '{if($2>0.1) print}' ${analysis}/TE_S_AS_bias.stats | grep Deug | while read TE AS_ratio lib; do
awk -v TE=${TE} -v LIB=${lib_sRNA} '{if(NR==12) print TE,$1,LIB}' ${analysis}/${lib_sRNA}/End.plots.TEs.linkage/${lib_sRNA}_${TE}_3MM_mappers.pingpong.minus_3end.minus_5end.table >> ${analysis}/phasing.stats.selected
done
### TE-piRNAs_analysis_plots.R was used to make a box plot in Figure 6.


### step 16: quantify miRNA reads, TE AS and S mapping piRNAs (gt22nt) from non-oxidised sRNAseq libraries
lib_sRNA="<unoxidied whole ovary small RNA library>"
lib_sRNA="Deug_ovary_unoxidised_sRNA"
miRNA=`(cat ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-mapped.bed | grep mir | awk '!seen[$4]++' |\
awk '{split($4,a,"@"); count+=a[2]} END {print count}')`

TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); {if($3-$2>22) count+=b[1]}} END {print count}')`

### quantify sense and antisense TE piRNAs
TE_S_all=`(bedtools intersect -s -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.bed \
-b ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna.out.allTE.bed | awk '!seen[$4]++' |\
awk '{split($4,a,"@"); split(a[2],b,":"); {if($3-$2>22) count+=b[1]}} END {print count}')`

TE_AS_all=`(bedtools intersect -S -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.bed \
-b ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna.out.allTE.bed | awk '!seen[$4]++' |\
awk '{split($4,a,"@"); split(a[2],b,":"); {if($3-$2>22) count+=b[1]}} END {print count}')`

TE_all=`(bedtools intersect -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Deug_UStan.all-best-strata.mapping-counts.tRNA-excluded.bed \
-b ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna.out.allTE.bed | awk '!seen[$4]++' |\
awk '{split($4,a,"@"); split(a[2],b,":"); {if($3-$2>22) count+=b[1]}} END {print count}')`

TE_S=`echo "$TE_all - $TE_AS_all" | bc`
TE_AS=`echo "$TE_all - $TE_S_all" | bc`
TE_S_AS=`echo "$TE_S_all + $TE_AS_all - $TE_all" | bc`
others=`echo "$TOTAL - $TE_all" | bc`

awk -v miRNA=${miRNA} -v TE_S=${TE_S} -v TE_AS=${TE_AS} -v TE_S_AS=${TE_S_AS} -v others=${others} -v LIB=${lib_sRNA} '{if(NR==1) print "TE_S",TE_S/miRNA,LIB"\n""TE_AS",TE_AS/miRNA,LIB"\n""TE_S_AS",TE_S_AS/miRNA,LIB"\n""others",others/miRNA,LIB
}' ${analysis}/TE_S_AS_bias.stats >> ${analysis}/miRNA_TE-mappers.stats
### TE-piRNAs_analysis_plots.R was used to make a box plot in Figure 6.


### step 17: assessing Piwi and Aub IP for the enrichment of somatic piRNAs
### combine the tile coverage data table into a single table using R
#R/4.0.0
library(reshape2)
table <- read.table("${analysis}/Deug/deug_0.5kb_tiles.counts", header=F)
table_d <- dcast(table,V1~V3,value.var="V2")
table_d[is.na(table_d)] <- 0
write.table(table_d, file="${analysis}/Deug/deug_0.5kb_tiles.counts.table", quote=F, row.names=F)
colnames(table_d)

[1] "V1"                        "Cluster2449_uni"
[3] "Cluster3031_uni"           "Cluster3378_dual"
[5] "Cluster3703_uni"           "Deug_Aub-IP_oxidised_sRNA"
[7] "Deug_embryo_oxidised_sRNA"             "Deug_ovary_oxidised_sRNA_rep1"
[9] "Deug_Piwi-IP_oxidised_sRNA"
#R/4.0.0

### we only included the tiles that are greater than 10 CPM in the whole ovary small RNA library
awk '{if(NR>1 && $8 > 10) print $1,$6/$8,"Aub-IP-per-ovary_total""\n"$1,$9/$8,"Piwi-IP-per-ovary_total";
if(NR>1 && $8 > 10 && $2+$3+$5>0) print $1,$6/$8,"Aub-IP-per-ovary_uni""\n"$1,$9/$8,"Piwi-IP-per-ovary_uni"
}' ${analysis}/Deug/deug_0.5kb_tiles.counts.table > ${analysis}/Deug/deug_0.5kb_tiles.counts.Aub-Piwi-IP.table
### TE-piRNAs_analysis_plots.R was used to make a box plot in Extended Data Figure 8.


### step 18: measure in-trans pingpong linkage between the whole ovary library and the Piwi and Aubergine IP libraries
for CLASS in unique all; do
for lib_sRNA in Deug_ovary_oxidised_sRNA_rep1 Deug_ovary_oxidised_sRNA_rep2 Deug_Aub-IP_oxidised_sRNA Deug_Piwi-IP_oxidised_sRNA; do
awk -v LIB=${lib_sRNA} '{print $0,LIB
}' ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_all_m9_dump.tab >> ${analysis}/Deug/Deug_${CLASS}_all_m9_dump.tab
done
done
### jellifish_linkage_IP.R was used to measure the linkage values


### step 19: collecting the in-trans pingpong linkage values across libraries
for CLASS in unique all; do
for lib_sRNA in w1118_ovary_oxidised_sRNA Dpse_ovary_oxidised_sRNA Deug_ovary_oxidised_sRNA_rep1 Deug_ovary_oxidised_sRNA_rep2 SRR1746887 SRR1104823; do
awk -v LIB=${lib_sRNA} -v CLASS=${CLASS} '{if(NR>1) print $1,$2,LIB"@"CLASS
}' ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_all_m9.linkage.txt >> ${analysis}/in-trans-pingpong.stats
done
done

for CLASS in unique all; do
for PIWI in Piwi-IP Aub-IP; do
awk -v PIWI=${PIWI} -v CLASS=${CLASS} '{print $1,$2,"Deug-ovary-vs"PIWI"@"CLASS
}' ${analysis}/Deug/Deug_${CLASS}_ovary-${PIWI}_linkage.txt >> ${analysis}/in-trans-pingpong.stats
done
done
### in-trans-pingpong_plots.R was used to make a barchart in Figure 6 and Extended Data Figure 9.

### post-processing of the analysed data for small RNA libraries --- END ---

### processing and analysing the small RNA sequencing libraries --- END ---



### processing and analysing the paired end RNA seq library of poly-A+ RNA --- START ---

### below is the RNA library processed in this part of the script:
Deug_ovary_polyA

### RNA seq libraries
5’ AATGATACGGCGACCACCGAGATCTACAC[TCTTTCCCTACACGACGCTCTTCCGATCT]--- [R1 ->] --- [<- R2] ---[AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC] ---NNNNNN--- ATCTCGTATGCCGTCTTCTGCTTG 3’
R1: reverse complement
R2: sense

### step 1: trim adapters, filter poor-quality reads, and take only the reads that have mates
lib_pA="Deug_ovary_polyA"
mkdir -p ${analysis}/${lib_pA}/STAR
fastx_clipper -l 35 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-i <(zcat ${raw_data}/${lib_pA}_R1.fastq.gz | fastq_quality_filter -q 33 -p 40) > ${analysis}/${lib_pA}/${lib_pA}_R1.trimmed.fastq
fastx_clipper -l 35 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGA \
-i <(zcat ${raw_data}/${lib_pA}_R2.fastq.gz | fastq_quality_filter -q 33 -p 40) > ${analysis}/${lib_pA}/${lib_pA}_R2.trimmed.fastq

### pair mates using fastq-pair and take only the reads that have mates
fastq_pair ${analysis}/${lib_pA}/${lib_pA}_R1.trimmed.fastq ${analysis}/${lib_pA}/${lib_pA}_R2.trimmed.fastq
for TYPE in R1 R2; do
gzip ${analysis}/${lib_pA}/${lib_pA}_${TYPE}.trimmed.fastq.paired.fq
rm ${analysis}/${lib_pA}/${lib_pA}_${TYPE}.trimmed.fastq.single.fq
rm ${analysis}/${lib_pA}/${lib_pA}_${TYPE}.trimmed.fastq
done

### step 2: make star index
mkdir -p ${references}/Deug/indices/Deug_UStan_STAR
STAR --runThreadN 4 --runMode genomeGenerate \
--genomeDir ${references}/Deug/indices/Deug_UStan_STAR \
--genomeSAindexNbases 12 \
--genomeFastaFiles ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.fna \
--sjdbGTFfile ${references}/Deug/GCF_018153835.1_ASM1815383v1_genomic.gtf \
--sjdbOverhang 100

### step 3: STAR mapping
STAR --runThreadN 8 --genomeDir ${references}/Deug/indices/Deug_UStan_STAR \
--outSAMattributes NH HI AS nM MD \
--outFilterMismatchNmax 3 \
--outSJfilterReads All \
--quantMode GeneCounts \
--outSAMreadID Number \
--readFilesCommand zcat \
--outFileNamePrefix ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_ \
--readFilesIn ${analysis}/${lib_pA}/${lib_pA}_R1.trimmed.fastq.paired.fq.gz ${analysis}/${lib_pA}/${lib_pA}_R2.trimmed.fastq.paired.fq.gz

### step 4: generate bigwig files for the IGV browser --- START ---
### take unique mappers and normalise the read counts
READS=`(awk '{if($0~"Uniquely mapped reads number") print $NF}' ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_Log.final.out)`

### sam to bed
samtools view -@ 8 -q 10 -bS ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_Aligned.out.sam |\
bamToBed -bed12 -i - > ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR.bed

### reverse the strandness of /1 and sort bed files
awk '{if($4~"/1" && $6=="-") {$6="+"; print}
else if($4~"/1" && $6=="+") {$6="-"; print}
else print}' ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR.bed |\
sort -k1,1 -k2,2 | tr ' ' '\t' > ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR.sorted.bed

### bedgraph output per 1Mio unique mappers (plus strand)
bedtools genomecov -strand + -split -bg -i ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR.sorted.bed \
-g ${references}/Deug/indices/Deug_UStan_STAR/chrNameLength.txt |\
awk -v READS=${READS} '{print $1,$2,$3,$4/READS*1000000}' | tr ' ' '\t' > ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_plus.bg
### bedgraph output per 1Mio unique mappers (minus strand)
bedtools genomecov -strand - -split -bg -i ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR.sorted.bed \
-g ${references}/Deug/indices/Deug_UStan_STAR/chrNameLength.txt |\
awk -v READS=${READS} '{print $1,$2,$3,-$4/READS*1000000}' | tr ' ' '\t' > ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_minus.bg

### generate bigWig files
for strand in plus minus; do
bedGraphToBigWig ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_${strand}.bg \
${references}/Deug/indices/Deug_UStan_STAR/chrNameLength.txt ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_${strand}.bw
done
### step 4: generate bigwig files for the IGV browser --- END ---
### processing and analysing the paired end RNA seq library of poly-A+ RNA --- END ---


### processing and analysing the Quant seq library --- START ---

### below is the QuantSeq library processed in this part of the script:
Deug_ovary_QuantSeq

lib_Q="Deug_ovary_QuantSeq"
mkdir -p ${analysis}/${lib_Q}/STAR

### step 1: adapter trimming
### we take R1 reads from the paired end sequencing library.
### -l 50: remove reads shorter than 50nt after adaptor trimming
fastx_clipper -l 50 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-i <(zcat ${raw_data}/${lib_Q}_R1.fastq.gz | fastq_quality_filter -q 33 -p 40) |\
fastx_clipper -a AAAAAAAAAAAAAAAAAAAA | gzip > ${analysis}/${lib_Q}/${lib_Q}_R1.trimmed.fastq.gz

### step 2: STAR mapping using the trimmed reads
STAR --runThreadN 8 --genomeDir ${references}/Deug/indices/Deug_UStan_STAR \
--outSAMattributes NH HI AS nM MD \
--outFilterMismatchNmax 3 \
--outSJfilterReads All \
--quantMode GeneCounts \
--outSAMreadID Number \
--readFilesCommand zcat \
--outFileNamePrefix ${analysis}/${lib_Q}/STAR/${lib_Q}_STAR_ \
--readFilesIn ${analysis}/${lib_Q}/${lib_Q}_R1.trimmed.fastq.gz

### step 3: generate bigwig files for the IGV browser --- START ---
### take unique mappers and normalise the read counts
READS=`(awk '{if($0~"Uniquely mapped reads number") print $NF}' ${analysis}/${lib_Q}/STAR/${lib_Q}_STAR_Log.final.out)`

### using samtools to extract genome unique mappers
samtools view -@ 8 -q 10 -bS ${analysis}/${lib_Q}/STAR/${lib_Q}_STAR_Aligned.out.sam |\
bamToBed -bed12 -i - |\
sort -k1,1 -k2,2 | tr ' ' '\t' > ${analysis}/${lib_Q}/STAR/${lib_Q}_STAR.bed

### bedgraph output per 1Mio unique mappers (plus strand)
### chrNameLength.txt is a tab separated file created while making the star index, and contains chromosome names and chromosome lengths.
bedtools genomecov -strand + -split -bg -i ${analysis}/${lib_Q}/STAR/${lib_Q}_STAR.bed \
-g ${references}/Deug/indices/Deug_UStan_STAR/chrNameLength.txt |\
awk -v READS=${READS} '{print $1,$2,$3,$4/READS*1000000}' | tr ' ' '\t' > ${analysis}/${lib_Q}/STAR/${lib_Q}_STAR_plus.bg
### bedgraph output per 1Mio unique mappers (minus strand)
bedtools genomecov -strand - -split -bg -i ${analysis}/${lib_Q}/STAR/${lib_Q}_STAR.bed \
-g ${references}/Deug/indices/Deug_UStan_STAR/chrNameLength.txt |\
awk -v READS=${READS} '{print $1,$2,$3,-$4/READS*1000000}' | tr ' ' '\t' > ${analysis}/${lib_Q}/STAR/${lib_Q}_STAR_minus.bg

### generate bigWig files
for strand in plus minus; do
bedGraphToBigWig ${analysis}/${lib_Q}/STAR/${lib_Q}_STAR_${strand}.bg \
${references}/Deug/indices/Deug_UStan_STAR/chrNameLength.txt ${analysis}/${lib_Q}/STAR/${lib_Q}_STAR_${strand}.bw
done
### step 3: generate bigwig files for the IGV browser --- END ---
### processing and analysing the Quant seq library --- START ---
