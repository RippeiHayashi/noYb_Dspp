### the code was used to analyse the sequencing data for Drosophila melanogaster.

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

### below are the small RNA libraries processed in this part of the script:
w1118_ovary_oxidised_sRNA
w1118_embryo_oxidised_sRNA
tj-g4_yb-shmiR_oxidised_sRNA
tj-g4_control_oxidised_sRNA

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
### step 0-2-1: download the genomic sequence of the dm6 r6.31 assembly of Drosophila melanogaster
mkdir -p ${references}/Dmel/indices
wget --directory-prefix="${references}/Dmel" http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.31_FB2019_06/fasta/dmel-all-chromosome-r6.31.fasta.gz
### only take the major chromosomes for the analysis
zcat ${references}/Dmel/dmel-all-chromosome-r6.31.fasta.gz | fasta_formatter - -t |\
awk '{if($1=="2L" || $1=="2R" || $1=="3L" || $1=="3R" || $1=="4" || $1=="X" || $1=="Y" ) print ">chr"$1"\n"$NF}' > ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta

### step 0-2-2: make a bowtie index
bowtie-build ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta ${references}/Dmel/indices/dmel-all-chromosome-r6.31.nosmallChr

### step 0-2-3: generate the size file
fasta_formatter -i ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta -t |\
awk '{print $1,length($NF)}' > ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.sizes


### step 0-3: generating a bed file for extended genomic regions of tRNAs
### step 0-3-1: download the fasta file of annotated tRNA sequences, which also contain the genomic coordinate information
wget --directory-prefix="${references}/Dmel" http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.31_FB2019_06/fasta/dmel-all-tRNA-r6.31.fasta.gz

### step 0-3-2: make a bed file of tRNA annotations extended by 100 bp both upstream and downstream
zcat ${references}/Dmel/dmel-all-tRNA-r6.31.fasta.gz | grep ">" | grep -v "mito" | awk 'BEGIN{ FS = ";"} {print $2}' | tr -d "loc=)" | tr ':(' ' ' |\
awk '{n=split($NF,a,"."); if($0 ~ "mpement") print "chr"$1,a[1]-101,a[n]+100,". . -";
else print "chr"$1,a[1]-101,a[n]+100,". . +"}' | sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/Dmel/dmel-all-tRNA-r6.31.tRNA-extended.bed


### step 0-4: generating a bed tile for genome unique 0.5kb tiles from the dm6 genome assembly of Drosophila melanogaster
### step 0-4-1: obtain every 25mer from the genome, only sense
### convert lowercase to uppercase first
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta |\
fasta_formatter - -t | awk '{for(i=0;i<=length($NF)-25;i++) print ">"$1"@"i"_+:"substr($NF,i+1,25)"\n"substr($NF,i+1,25)
}' > ${references}/Dmel/Dmel_dm6_25mers.plus.fasta

### step 0-4-2: obtain genome unique mappers
index_genome="${references}/Dmel/indices/dmel-all-chromosome-r6.31.nosmallChr"
bowtie -f -v 0 -m 1 -S ${index_genome} ${references}/Dmel/Dmel_dm6_25mers.plus.fasta |\
samtools view -@ 16 -bS - | bamToBed -i - > ${references}/Dmel/Dmel_dm6_25mers.plus.unique-mappers.bed

### step 0-4-3: count the unique mappers in 0.5kb bins
### take if the centre of 25mer sits between the coordinate of 1 to 500 as the first 0.5kb bin
### take as the second 0.5kb bin if it sits between 501 to 1000
awk '{split($4,a,"@"); split(a[2],b,"_"); TILE[$1"_"int((b[1]+12)/500)]++
} END {for(var in TILE) print var,TILE[var]}' ${references}/Dmel/Dmel_dm6_25mers.plus.unique-mappers.bed > ${references}/Dmel/Dmel_dm6_0.5kbtiles.counts

### step 0-4-4: make a bed file that contains all kb tiles that are gt 85% mappability
awk '{split($1,a,"_"); if($2>425) print a[1],a[2]*500,(a[2]+1)*500,a[1]":"a[2]*500"-"(a[2]+1)*500"_+ . +""\n"a[1],a[2]*500,(a[2]+1)*500,a[1]":"a[2]*500"-"(a[2]+1)*500"_- . -"
}' ${references}/Dmel/Dmel_dm6_0.5kbtiles.counts |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/Dmel/Dmel_dm6_85percent_0.5kbtiles.25mers.bed


### step 0-5: preparing for the 0.2kb genomic tile analyses of the all genome mappers --- START ---
### step 0-5-1: make a bed file of 0.2kb genomic tiles
awk '{for(i=0;i<=int($2/200)-1;i++) {print $1,i*200,(i+1)*200,$1":"i*200"-"(i+1)*200":- . -""\n"$1,i*200,(i+1)*200,$1":"i*200"-"(i+1)*200":+ . +"
} {print $1,int($2/200)*200,$2,$1":"int($2/200)*200"-"$2":- . -""\n"$1,int($2/200)*200,$2,$1":"int($2/200)*200"-"$2":+ . +"}
}' ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.sizes |\
tr ' ' '\t' > ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.200nt.window.bed

### step 0-5-2: add annotations of piRNA clusters, gypsy and TE insertions, and exons
mkdir -p ${analysis}/Dmel
### step 0-5-2-1: piRNA clusters
bedtools intersect -wo -s -a ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.200nt.window.bed \
-b ${references}/Dmel/Dmel_dm6_piRNA-clusters.bed |\
awk '{print $4,$NF/200,$10}' > ${analysis}/Dmel/Dmel_0.2kb_tiles.counts
### Dmel_dm6_piRNA-clusters.bed looks as follows:
chrX    21631891        22139306        ClusterFlam     0       +
chr3L   23280548        23311746        Cluster80F      0       +
chr3L   23280548        23311746        Cluster80F      0       -

### step 0-5-2-2: Run RepeatMasker 4.1.0 by specifying "drosophila"
fastafile="${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta"
RepeatMasker -species "drosophila" -parallel 10 -xsmall ${fastafile}

### step 0-5-2-3: extract gypsy insertions --- START ---
cat ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta.out | grep -i "gypsy" |\
awk '{if($9=="+") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9; else if($9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"
}' | tr ' ' '\t' > ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta.out.gypsy.bed

### extract sense and antisense gypsy insertions, these bed files are used in the IGV browser tracks
cat ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta.out | grep -i "gypsy" |\
awk '{if($9=="+" ) print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9}' |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta.out.gypsy.plus.bed
cat ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta.out | grep -i "gypsy" |\
awk '{if($9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"}' |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta.out.gypsy.minus.bed

### annotate tiles that overlap with gypsy insertions at the sense orientation
bedtools intersect -f 0.5 -s -a ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.200nt.window.bed \
-b ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta.out.gypsy.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"gypsy_S"}' >> ${analysis}/Dmel/Dmel_0.2kb_tiles.counts

### annotate tiles that overlap with gypsy insertions at the antisense orientation
bedtools intersect -f 0.5 -S -a ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.200nt.window.bed \
-b ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta.out.gypsy.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"gypsy_AS"}' >> ${analysis}/Dmel/Dmel_0.2kb_tiles.counts
### step 0-5-2-3: extract gypsy insertions --- END ---

### step 0-5-2-4: extract TE insertions --- START ---
### remove "Simple_repeat", "Low_complexity", "rRNA" and "tRNA"
awk '{if($11 != "Simple_repeat" && $11 != "Low_complexity" && $11 != "rRNA" && $11 != "tRNA") print}' ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta.out |\
awk '{if($9=="+") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9; else if($9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"
}' | tr ' ' '\t' > ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta.out.allTE.bed

### annotate tiles that overlap with TE insertions at the sense orientation
bedtools intersect -f 0.5 -s -a ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.200nt.window.bed \
-b ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta.out.allTE.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"TE_S"}' >> ${analysis}/Dmel/Dmel_0.2kb_tiles.counts

### annotate tiles that overlap with TE insertions at the antisense orientation
bedtools intersect -f 0.5 -S -a ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.200nt.window.bed \
-b ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta.out.allTE.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"TE_AS"}' >> ${analysis}/Dmel/Dmel_0.2kb_tiles.counts
### step 0-5-2-4: extract TE insertions --- END ---

### step 0-5-2-5: add exon annotations --- START ---
### step 0-5-2-5-1: download exon coordinates from the Flybase
wget --directory-prefix="${references}/Dmel" http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.31_FB2019_06/fasta/dmel-all-exon-r6.31.fasta.gz

### step 0-5-2-5-2: make a bed file of annotated exons while merging overlapping annotations
### plus strand
zcat ${references}/Dmel/dmel-all-exon-r6.31.fasta.gz | sed 's/MD5/ MD5/g' | grep ">" | grep -v "comple" | awk '{print $3}' |
tr '(:)' ' ' | tr -d 'loc=join;complement' |\
awk '{n=split($2,a,","); if(NF==2) for(i=1;i<=n;i++) print $1" "a[i]}' | tr '.' ' ' |\
### flybase annotation is a 1-based bed format
awk '{if($1=="2L" || $1=="2R" || $1=="3L" || $1=="3R" || $1=="X" || $1=="4" || $1=="Y") print "chr"$1,$2-1,$3,"1 1 +"}' | sort -k1,1 -k2,2n | tr ' ' '\t' |\
bedtools merge -i - | awk '{print $0,"1 1 +"}' | tr ' ' '\t' > ${references}/Dmel/dmel-all-exon-r6.31_merged.bed
### minus strand
zcat ${references}/Dmel/dmel-all-exon-r6.31.fasta.gz | sed 's/MD5/ MD5/g' | grep ">" | grep "comple" | awk '{print $3}' |\
tr '(:)' ' ' | tr -d 'loc=join;complement' |\
awk '{n=split($2,a,","); if(NF==2) for(i=1;i<=n;i++) print $1" "a[i]}' | tr '.' ' ' |\
### flybase annotation is a 1-based bed format
awk '{if($1=="2L" || $1=="2R" || $1=="3L" || $1=="3R" || $1=="X" || $1=="4" || $1=="Y") print "chr"$1,$2-1,$3,"1 1 -"}' | sort -k1,1 -k2,2n | tr ' ' '\t' |\
bedtools merge -i - | awk '{print $0,"1 1 -"}' | tr ' ' '\t' >> ${references}/Dmel/dmel-all-exon-r6.31_merged.bed
### sort the bed file for downstream intersection
cat ${references}/Dmel/dmel-all-exon-r6.31_merged.bed | sort -k1,1 -k2,2n > ${references}/Dmel/dmel-all-exon-r6.31_merged_sorted.bed

### step 0-5-2-5-3: intersect exon annotations to the 0.2kb tiles
bedtools intersect -wo -s -a ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.200nt.window.bed \
-b ${references}/Dmel/dmel-all-exon-r6.31_merged_sorted.bed |\
awk '{exon[$4]+=$NF} END {for(var in exon) print var,exon[var]/200,"exon"}' >> ${analysis}/Dmel/Dmel_0.2kb_tiles.counts
### step 0-5-2-5: add exon annotations --- END ---


### step 0-6: make a bowtie index for the Drosophila melanogaster transposon sequences used in Senti G&D 2015, PMID: 26302790 (122 entries)
bowtie-build ${references}/TEs/122_dmel_TE_Senti2015.fa ${references}/TEs/indices/122_dmel_TE_Senti2015
### create the size file
cat ${references}/TEs/122_dmel_TE_Senti2015.fa | fasta_formatter - -t | awk '{print $1,length($NF)}' | tr ' ' '\t' > ${references}/TEs/122_dmel_TE_Senti2015.sizes


### step 0-7: obtain TEs that produce abundant piRNAs
### list 24 TEs that produce least abundant piRNAs in the whole ovary
lib_sRNA="w1118_ovary_oxidised_sRNA"
awk '!seen[$1" "$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_dm6_mapped_TE_3MM_mappers.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if(length(a[1])>22) TE[$1]+=b[1]} END {for(var in TE) print var,TE[var]}' |\
sort -k2,2n | head -n 24 > ${analysis}/${lib_sRNA}/${lib_sRNA}_dm6_mapped_TE_3MM_mappers.low-abundant.list
### remove those 24 TEs from downstream analyses
join <(sort -k1,1 ${references}/TEs/122_dmel_TE_Senti2015.sizes) <(sort -k1,1 ${analysis}/${lib_sRNA}/${lib_sRNA}_dm6_mapped_TE_3MM_mappers.low-abundant.list) -a 1 |\
awk '{if(NF==2) print}' > ${references}/TEs/abundant-98_dmel_TE_Senti2015.sizes


### step 0-8: generate pseudocounts of 1 at every coordinate of non-redundant TEs
mkdir -p ${references}/TEs/TE_pseudo.piRNAs/Dmel
cat ${references}/TEs/abundant-98_dmel_TE_Senti2015.sizes | while read TE length; do
if [[ ! -f ${references}/TEs/TE_pseudo.piRNAs/Dmel/${TE}_pseudo.piRNAs.txt ]]; then
### generate pseudocounts of 1 at every coordinate
for ((i=0; i<=${length}-1; i++)); do
j=$((${i}+1))
printf "${TE} ${i} ${j} AAAAAAAAAAAAAAAAAAAAAAA@1 . +""\n""${TE} ${i} ${j} AAAAAAAAAAAAAAAAAAAAAAA@1 . -""\n" |\
tr ' ' '\t' >> ${references}/TEs/TE_pseudo.piRNAs/Dmel/${TE}_pseudo.piRNAs.txt
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
index_genome="${references}/Dmel/indices/dmel-all-chromosome-r6.31.nosmallChr"
bowtie -f -v 1 -m 1 -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - |\
awk '{if($1=="chr2L" || $1=="chr2R" || $1=="chr3L" || $1=="chr3R" || $1=="chr4" || $1=="chrX") print}' |\
sort -k1,1 -k2,2n > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dm6.bed


### step 4: remove tRNA reads by intersection
### exclude reads that intersect with tRNA annotations from the genome unique mappers
bedtools intersect -v -s -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dm6.bed \
-b ${references}/Dmel/dmel-all-tRNA-r6.31.tRNA-extended.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dm6.tRNA-excluded.bed


### step 6: 0.5kb tile analysis for the genome unique mappers
### count piRNA reads (greater than 22nt) that uniquely mapped to the 0.5kb tiles
READS=`(awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]} END {print count
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dm6.tRNA-excluded.bed)`
bedtools intersect -wo -s -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dm6.tRNA-excluded.bed \
-b ${references}/Dmel/Dmel_dm6_85percent_0.5kbtiles.25mers.bed |\
awk -v READS=${READS} -v LIB=${lib_sRNA} '{split($4,a,"@"); if(length(a[1])>22) TILE[$10]+=a[2]} END {for(var in TILE) print var,TILE[var]/READS*1000000,LIB
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dm6.tRNA-excluded.gt22nt.0.5kbtiles.counts


### step 7: 0.2kb tile analysis for the genome all mappers --- START ---
### step 7-1:map reads to dm6, allowing up to 1MM, all best strata
index_genome="${references}/Dmel/indices/dmel-all-chromosome-r6.31.nosmallChr"
bowtie -f -v 1 --all --best --strata -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.bed

### step 7-2: add mapping counts to the fasta entries. Mapping counts will be used to distribute the read coverage across tiles.
awk '{READS[$4]++} END {for(var in READS) {split(var,a,"@"); print ">"var":"READS[var]"\n"a[1]}
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.mapping-counts.fa

### step 7-3: map again all best strata
index_genome="${references}/Dmel/indices/dmel-all-chromosome-r6.31.nosmallChr"
bowtie -f -v 1 --all --best --strata -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.mapping-counts.fa |\
samtools view -bS - | bamToBed -i - | sort -k1,1 -k2,2n > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.mapping-counts.bed

### step 7-4: exclude reads that intersect with tRNA annotations
bedtools intersect -v -s -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.mapping-counts.bed \
-b ${references}/Dmel/dmel-all-tRNA-r6.31.tRNA-extended.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.mapping-counts.tRNA-excluded.bed

### step 7-5: collect read counts per 200nt windows from all-best-strata
TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if($3-$2>22) count+=b[1]} END {print count}')`
awk 'BEGIN{OFS="\t"} {if($3-$2 > 22) print}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.mapping-counts.tRNA-excluded.bed |\
bedtools intersect -s -wo -F 0.5 -a ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.200nt.window.bed -b - |\
awk -v TOTAL=${TOTAL} -v LIB=${lib_sRNA} '{split($10,a,"@"); split(a[2],b,":"); TILE[$4]+=b[1]/b[2]} END {for(var in TILE) print var,TILE[var]/TOTAL*1000000,LIB"_abs"
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.mapping-counts.tRNA-excluded.gt22nt.0.2kbtiles.counts
### step 7: 0.2kb tile analysis for the genome all mappers --- END ---


### step 8: analyse TE mapping piRNAs --- START ---
### step 8-1: take all reads that mapped to the D.melanogaster genome, and use them to map against transposon sequences
awk '{READ[$4]} END {for(var in READ) {split(var,a,"@"); print ">"var"\n"a[1]}
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.mapping-counts.tRNA-excluded.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.mapping-counts.tRNA-excluded.fa

### step 8-2: bowtie mapping to TEs, --all --best --strata option, allowing up to three mismatches
TE_index="${references}/TEs/indices/122_dmel_TE_Senti2015"
bowtie -f -v 3 --all --best --strata -S ${TE_index} ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.mapping-counts.tRNA-excluded.fa |\
samtools view -b -S - | bedtools bamtobed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_dm6_mapped_TE_3MM_mappers.bed

### step 8-3: making TE mapping plots and measuring pingpong and phasing signatures --- START ---
mkdir -p ${analysis}/${lib_sRNA}/End.plots.TEs.linkage
### step 8-3-1: normalise the counts of 3' end 5' ends per 1 million total genome mappers (>22nt)
TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if(length(a[1])>22) count+=b[1]} END {print count}')`

### step 8-3-2: count 5end and 3end of plus and minus mapping reads (>22nt)
### run it per each TE
cat ${references}/TEs/abundant-98_dmel_TE_Senti2015.sizes | while read TE length; do
### use pseudocounts created in step 0-8
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_dm6_mapped_TE_3MM_mappers.bed ${references}/TEs/TE_pseudo.piRNAs/Dmel/${TE}_pseudo.piRNAs.txt |\
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


### step 9: measure in-trans ping-pong linkage --- START ---
mkdir -p ${analysis}/${lib_sRNA}/jellyfish

### step 9-1: extract g1g9, g2g10_revComp and last9 sequences from piRNA genome unique mappers
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr($4,1,9)
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dm6.tRNA-excluded.bed > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_unique_g1g9.fasta
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr($4,2,9)
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dm6.tRNA-excluded.bed | fastx_reverse_complement - > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_unique_g2g10_revComp.fasta
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr(a[1],length(a[1])-8,9)
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dm6.tRNA-excluded.bed > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_unique_last9.fasta

### step 9-2: extract g1g9, g2g10_revComp and last9 sequences from piRNA genome all mappers
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if($3-$2>22) for(i=1;i<=b[1];i++) print ">"$4":"i"\n"substr($4,1,9)}'  > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_all_g1g9.fasta
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if($3-$2>22) for(i=1;i<=b[1];i++) print ">"$4":"i"\n"substr($4,2,9)}' | fastx_reverse_complement - > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_all_g2g10_revComp.fasta
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.mapping-counts.tRNA-excluded.bed |\
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
mkdir -p ${analysis}/Dmel
bedtools intersect -wo -s -a ${references}/Dmel/Dmel_dm6_85percent_0.5kbtiles.25mers.bed \
-b ${references}/Dmel/Dmel_dm6_piRNA-clusters.bed |\
awk '{print $4,$NF/500,$10}' > ${analysis}/Dmel/dmel_0.5kb_tiles.counts

### step 10-2: add counts from small RNA libraries
for lib_sRNA in <oxidied whole ovary small RNA library> <oxidied embryonic small RNA library>; do
#for lib_sRNA in w1118_ovary_oxidised_sRNA w1118_embryo_oxidised_sRNA; do
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dm6.tRNA-excluded.gt22nt.0.5kbtiles.counts >> ${analysis}/Dmel/dmel_0.5kb_tiles.counts
done
### step 10-3: tile_analysis_plots.R was used to make a scatter plot in Figure 5.


### step 11: 0.2kb tile analysis of the genome all mappers --- START ---
### step 11-1: add counts from small RNA libraries
for lib_sRNA in <w1118 whole ovary small RNA library> <w1118 embryonic small RNA library> <tj-g4_yb-shmiR whole ovary small RNA library> <tj-g4_control whole ovary small RNA library>; do
#for lib_sRNA in w1118_ovary_oxidised_sRNA w1118_embryo_oxidised_sRNA tj-g4_yb-shmiR_oxidised_sRNA tj-g4_control_oxidised_sRNA; do
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.mapping-counts.tRNA-excluded.gt22nt.0.2kbtiles.counts >> ${analysis}/Dmel/Dmel_0.2kb_tiles.counts
done

### step 11-2: make a table using reshape2 in R and count piRNA reads from individual annotations
#R/4.0.0
library(reshape2)
table <- read.table("/g/data/lf10/rh1772/noYb_paper/analysis/Dmel/Dmel_0.2kb_tiles.counts", header=F)
table_d=dcast(table,V1~V3,value.var="V2")
table_d[is.na(table_d)] <- 0
write.table(table_d, file="/g/data/lf10/rh1772/noYb_paper/analysis/Dmel/Dmel_0.2kb_tiles.counts.table", quote=F, col.names=T, row.names = F)
colnames(table_d)

[1] "V1"                             "Cluster80F"
[3] "ClusterFlam"                    "exon"
[5] "gypsy_AS"                       "gypsy_S"
[7] "TE_AS"                          "TE_S"
[9] "tj-g4_control_oxidised_sRNA_abs"      "tj-g4_yb-shmiR_oxidised_sRNA_abs"
[11] "w1118_embryo_oxidised_sRNA_abs"     "w1118_ovary_oxidised_sRNA_abs"
#R/4.0.0

### w1118 --- START ---
### count total somatic reads
soma_total=`(awk '{if(NR>1 && $12>0) print}' ${analysis}/Dmel/Dmel_0.2kb_tiles.counts.table |\
awk '{if($11/$12 < 0.1) print}' | awk '{count+=$12} END {print count}')`

### count somatic uni-stranded cluster reads
uni_total=`(awk '{if(NR>1 && $12>0) print}' ${analysis}/Dmel/Dmel_0.2kb_tiles.counts.table |\
awk '{if($11/$12 < 0.1) print}' | awk '{if($3>0) count+=$12} END {print count}')`

### count exonic piRNA reads
### exclude tiles that overlap with uni clusters and TE annotations
exons=`(awk '{if(NR>1 && $12>0) print}' ${analysis}/Dmel/Dmel_0.2kb_tiles.counts.table |\
awk '{if($11/$12 < 0.1) print}' | awk '{if($3==0 && $7==0 && $8==0 && $4>0) count+=$12} END {print count}')`

### count somatic gypsy_AS reads outside the clusters
no_uni_gypsy_AS=`(awk '{if(NR>1 && $12>0) print}' ${analysis}/Dmel/Dmel_0.2kb_tiles.counts.table |\
awk '{if($11/$12 < 0.1) print}' | awk '{if($3==0 && $5>0 && $6==0) count+=$12} END {print count}')`

### count somatic gypsy_S reads outside the clusters
no_uni_gypsy_S=`(awk '{if(NR>1 && $12>0) print}' ${analysis}/Dmel/Dmel_0.2kb_tiles.counts.table |\
awk '{if($11/$12 < 0.1) print}' | awk '{if($3==0 && $6>0 && $5==0) count+=$12} END {print count}')`

### collect the stats
others=`echo "$soma_total - $uni_total - $exons - $no_uni_gypsy_AS - $no_uni_gypsy_S" | bc`
printf $soma_total" w1118 soma_total""\n"$uni_total" w1118 uni_total""\n"$exons" w1118 exons""\n"$no_uni_gypsy_AS" w1118 no_uni_gypsy_AS""\n"$no_uni_gypsy_S" w1118 no_uni_gypsy_S""\n"$others" w1118 others""\n" >> ${analysis}/0.2kb_tiles.counts.stats
### w1118 --- END ---

### tj-g4 x yb shmiR --- START ---
### count total somatic reads
tjg4yb_soma_total=`(awk '{if(NR>1 && $12>0) print}' ${analysis}/Dmel/Dmel_0.2kb_tiles.counts.table |\
awk '{if($11/$12 < 0.1) print}' | awk '{count+=$10} END {print count}')`

### count somatic uni-stranded cluster reads
tjg4yb_uni_total=`(awk '{if(NR>1 && $12>0) print}' ${analysis}/Dmel/Dmel_0.2kb_tiles.counts.table |\
awk '{if($11/$12 < 0.1) print}' | awk '{if($3>0) count+=$10} END {print count}')`

### count exonic piRNA reads
### exclude tiles that overlap with uni clusters and TE annotations
tjg4yb_exons=`(awk '{if(NR>1 && $12>0) print}' ${analysis}/Dmel/Dmel_0.2kb_tiles.counts.table |\
awk '{if($11/$12 < 0.1) print}' | awk '{if($3==0 && $7==0 && $8==0 && $4>0) count+=$10} END {print count}')`

### count somatic gypsy_AS reads outside the clusters
tjg4yb_no_uni_gypsy_AS=`(awk '{if(NR>1 && $12>0) print}' ${analysis}/Dmel/Dmel_0.2kb_tiles.counts.table |\
awk '{if($11/$12 < 0.1) print}' | awk '{if($3==0 && $5>0 && $6==0) count+=$10} END {print count}')`

### count somatic gypsy_S reads outside the clusters
tjg4yb_no_uni_gypsy_S=`(awk '{if(NR>1 && $12>0) print}' ${analysis}/Dmel/Dmel_0.2kb_tiles.counts.table |\
awk '{if($11/$12 < 0.1) print}' | awk '{if($3==0 && $6>0 && $5==0) count+=$10} END {print count}')`

### collect the stats
tjg4yb_others=`echo "$tjg4yb_soma_total - $tjg4yb_uni_total - $tjg4yb_exons - $tjg4yb_no_uni_gypsy_AS - $tjg4yb_no_uni_gypsy_S" | bc`
printf $tjg4yb_soma_total" tjg4yb soma_total""\n"$tjg4yb_uni_total" tjg4yb uni_total""\n"$tjg4yb_exons" tjg4yb exons""\n"$tjg4yb_no_uni_gypsy_AS" tjg4yb no_uni_gypsy_AS""\n"$tjg4yb_no_uni_gypsy_S" tjg4yb no_uni_gypsy_S""\n"$tjg4yb_others" tjg4yb others""\n" >> ${analysis}/0.2kb_tiles.counts.stats
### tj-g4 x yb shmiR --- END ---


### tj-g4 x control --- START ---
### count total somatic reads
tjg4ctl_soma_total=`(awk '{if(NR>1 && $12>0) print}' ${analysis}/Dmel/Dmel_0.2kb_tiles.counts.table |\
awk '{if($11/$12 < 0.1) print}' | awk '{count+=$9} END {print count}')`

### count somatic uni-stranded cluster reads
tjg4ctl_uni_total=`(awk '{if(NR>1 && $12>0) print}' ${analysis}/Dmel/Dmel_0.2kb_tiles.counts.table |\
awk '{if($11/$12 < 0.1) print}' | awk '{if($3>0) count+=$9} END {print count}')`

### count exonic piRNA reads
### exclude tiles that overlap with uni clusters and TE annotations
tjg4ctl_exons=`(awk '{if(NR>1 && $12>0) print}' ${analysis}/Dmel/Dmel_0.2kb_tiles.counts.table |\
awk '{if($11/$12 < 0.1) print}' | awk '{if($3==0 && $7==0 && $8==0 && $4>0) count+=$9} END {print count}')`

### count somatic gypsy_AS reads outside the clusters
tjg4ctl_no_uni_gypsy_AS=`(awk '{if(NR>1 && $12>0) print}' ${analysis}/Dmel/Dmel_0.2kb_tiles.counts.table |\
awk '{if($11/$12 < 0.1) print}' | awk '{if($3==0 && $5>0 && $6==0) count+=$9} END {print count}')`

### count somatic gypsy_S reads outside the clusters
tjg4ctl_no_uni_gypsy_S=`(awk '{if(NR>1 && $12>0) print}' ${analysis}/Dmel/Dmel_0.2kb_tiles.counts.table |\
awk '{if($11/$12 < 0.1) print}' | awk '{if($3==0 && $6>0 && $5==0) count+=$9} END {print count}')`

### collect the stats
tjg4ctl_others=`echo "$tjg4ctl_soma_total - $tjg4ctl_uni_total - $tjg4ctl_exons - $tjg4ctl_no_uni_gypsy_AS - $tjg4ctl_no_uni_gypsy_S" | bc`
printf $tjg4ctl_soma_total" tjg4ctl soma_total""\n"$tjg4ctl_uni_total" tjg4ctl uni_total""\n"$tjg4ctl_exons" tjg4ctl exons""\n"$tjg4ctl_no_uni_gypsy_AS" tjg4ctl no_uni_gypsy_AS""\n"$tjg4ctl_no_uni_gypsy_S" tjg4ctl no_uni_gypsy_S""\n"$tjg4ctl_others" tjg4ctl others""\n" >> ${analysis}/0.2kb_tiles.counts.stats
### tj-g4 x control --- END ---
### step 11: 0.2kb tile analysis of the genome all mappers --- END ---


### step 12: comparing the abundance of somatic and germline TE mapping piRNA reads (>22nt) --- START ---
for lib_sRNA in <oxidied whole ovary small RNA library> <oxidied embryonic small RNA library>; do
#for lib_sRNA in w1118_ovary_oxidised_sRNA w1118_embryo_oxidised_sRNA; do
### normalise the counts to 1 million total genome mappers (>22nt)
TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if(length(a[1])>22) count+=b[1]} END {print count}')`

### step 12-1: count piRNA mappers (>22nt) both sense and antisense
### only count the top 98TEs
awk '!seen[$1" "$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_dm6_mapped_TE_3MM_mappers.bed |\
awk -v TOTAL=${TOTAL} -v LIB=${lib_sRNA} '{split($4,a,"@"); split(a[2],b,":"); if(length(a[1])>22) TE[$1]+=b[1]} END {for(var in TE) print var,TE[var]/TOTAL*1000000,LIB
}' | sort -k1,1 | join - ${references}/TEs/abundant-98_dmel_TE_Senti2015.sizes | awk '{print $1,$2,$3}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_98TE_mappers.counts

### step 12-2: combine counts from the whole ovary and embryonic small RNA libraries
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_98TE_mappers.counts >> ${analysis}/Dmel_98TE_mappers.counts
done

### step 12-3: TE-piRNAs_analysis_plots.R was used to make a scatter plot in Extended Data Figure 6.

### step 12-4: identify TEs that are enriched more than three times in the ovary sRNAseq library than in the embryonic sRNAseq library
#R/4.0.0
library(reshape2)
table <- read.table("${analysis}/Dmel_98TE_mappers.counts", header=F)
table.d <- dcast(table,V1~V3,value.var="V2")
write.table(table.d, file="${analysis}/Dmel_98TE_mappers.counts.table", quote=F, row.names=F)
#R/4.0.0

awk '{if(NR>1 && $3>$2*3) print $1}' ${analysis}/Dmel_98TE_mappers.counts.table
412
gtwin
gypsy
gypsy10
gypsy5
Tabor
ZAM
### step 12: comparing the abundance of somatic and germline TE mapping piRNA reads (>22nt) --- END ---


### step 13: measuring the strand bias of the piRNA reads (>22nt) mapping to the abundant 98 TEs
### exclude somatic TEs
lib_sRNA="<oxidied whole ovary small RNA library>"
lib_sRNA="w1118_ovary_oxidised_sRNA"
awk '{if($1 != "412" && $1 != "gtwin" && $1 != "gypsy" && $1 != "gypsy10" && $1 != "gypsy5" && $1 != "Tabor" && $1 != "ZAM") print
}' ${references}/TEs/abundant-98_dmel_TE_Senti2015.sizes | while read TE READS; do
## some TE sequences have repeated sequences within and the same reads can map multiple times
awk '!seen[$1" "$4" "$6]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_dm6_mapped_TE_3MM_mappers.bed |\
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


### step 16: quantify miRNA reads, TE AS and S mapping piRNAs (gt22nt) from a non-oxidised sRNAseq library from Hayashi, 2016, PMID: 27851737 --- START ---
lib_sRNA="<unoxidied whole ovary small RNA library>"
### step 16-1: download the fastq file and trim the adapter sequence
wget --directory-prefix="${raw_data}" ftp.sra.ebi.ac.uk/vol1/fastq/SRR371/008/SRR3715418/SRR3715418.fastq.gz

lib_sRNA="SRR3715418"
mkdir -p ${analysis}/${lib_sRNA}
fastq_to_fasta -Q33 -i <(zcat ${raw_data}/${lib_sRNA}.fastq.gz) |\
# -c: discard non-clipped sequences, -l 18: discard sequences shorter than 18nt
# collapse reads of the same sequence to one entry while retaining the number of reads
fastx_clipper -c -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -l 18 -i - | fasta_formatter - -t |\
awk '{if(length($NF)>25 && length($NF)<49) READS[substr($NF,5,length($NF)-8)]++} END {
for(var in READS) print ">"var"@"READS[var]"\n"var}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa

### step 16-2: run bowtie to map reads to the miscRNA and take unmapped reads
misc_index="${references}/dm6/indices/dmel-all-miRNA-rRNA-snRNA-snoRNA-tRNA-r6.31_spikes"
bowtie -f -v 1 -a -S ${misc_index} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa \
--un ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-mapped.bed

### step 16-3: map reads to dm6, allowing up to 1MM, all best strata --- START ---
index_genome="${references}/Dmel/indices/dmel-all-chromosome-r6.31.nosmallChr"
bowtie -f -v 1 --all --best --strata -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.bed

### step 16-4: exclude reads that intersect with tRNA annotations
bedtools intersect -v -s -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.bed \
-b ${references}/Dmel/dmel-all-tRNA-r6.31.tRNA-extended.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.tRNA-excluded.bed

### step 16-5: count miRNA reads
miRNA=`(cat ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-mapped.bed | grep mir | awk '!seen[$4]++' |\
awk '{split($4,a,"@"); count+=a[2]} END {print count}')`

TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); {if($3-$2>22) count+=b[1]}} END {print count}')`

### step 16-6: quantify sense and antisense TE piRNAs
TE_S_all=`(bedtools intersect -s -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.tRNA-excluded.bed \
-b ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta.out.allTE.bed | awk '!seen[$4]++' |\
awk '{split($4,a,"@"); split(a[2],b,":"); {if($3-$2>22) count+=b[1]}} END {print count}')`

TE_AS_all=`(bedtools intersect -S -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.tRNA-excluded.bed \
-b ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta.out.allTE.bed | awk '!seen[$4]++' |\
awk '{split($4,a,"@"); split(a[2],b,":"); {if($3-$2>22) count+=b[1]}} END {print count}')`

TE_all=`(bedtools intersect -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.tRNA-excluded.bed \
-b ${references}/Dmel/dmel-all-chromosome-r6.31.nosmallChr.fasta.out.allTE.bed | awk '!seen[$4]++' |\
awk '{split($4,a,"@"); split(a[2],b,":"); {if($3-$2>22) count+=b[1]}} END {print count}')`

TE_S=`echo "$TE_all - $TE_AS_all" | bc`
TE_AS=`echo "$TE_all - $TE_S_all" | bc`
TE_S_AS=`echo "$TE_S_all + $TE_AS_all - $TE_all" | bc`
others=`echo "$TOTAL - $TE_all" | bc`

awk -v miRNA=${miRNA} -v TE_S=${TE_S} -v TE_AS=${TE_AS} -v TE_S_AS=${TE_S_AS} -v others=${others} -v LIB=${lib_sRNA} '{if(NR==1) print "TE_S",TE_S/miRNA,LIB"\n""TE_AS",TE_AS/miRNA,LIB"\n""TE_S_AS",TE_S_AS/miRNA,LIB"\n""others",others/miRNA,LIB
}' ${analysis}/TE_S_AS_bias.stats >> ${analysis}/miRNA_TE-mappers.stats
### TE-piRNAs_analysis_plots.R was used to make a box plot in Figure 6.
### step 16: quantify miRNA reads, TE AS and S mapping piRNAs (gt22nt) from a non-oxidised sRNAseq library from Hayashi, 2016, PMID: 27851737 --- END ---


### step 17: analysing piRNAs from the ovarian somatic cells (OSCs) for the in-trans pingpong analysis --- START ---
### step 17-1: downloading small RNA libraries of ovarian somatic cells (OSCs) from Mohn, 2015, PMID: 25977553
wget --directory-prefix="${raw_data}" ftp.sra.ebi.ac.uk/vol1/fastq/SRR174/007/SRR1746887/SRR1746887.fastq.gz

### Step 17-2. clip adapters and collapse them by sequence
lib_sRNA="SRR1746887"
mkdir -p ${analysis}/${lib_sRNA}
fastq_to_fasta -Q33 -i <(zcat ${raw_data}/${lib_sRNA}.fastq.gz) |\
# -c: discard non-clipped sequences, -l 18: discard sequences shorter than 18nt
fastx_clipper -c -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -l 18 -i - | fasta_formatter - -t |\
awk '{if(length($NF)>25 && length($NF)<49) READS[substr($NF,5,length($NF)-8)]++} END {
for(var in READS) print ">"var"@"READS[var]"\n"var}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa

### step 17-3: map to the Dmel misc RNA allowing up to 1MM --- START ---
misc_index="${references}/dm6/indices/dmel-all-miRNA-rRNA-snRNA-snoRNA-tRNA-r6.31_spikes"
bowtie -f -v 1 -a -S ${misc_index} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa \
--un ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-mapped.bed

### step 17-4: genome unique mappers
index_genome="${references}/Dmel/indices/dmel-all-chromosome-r6.31.nosmallChr"
bowtie -f -v 1 -m 1 -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - |\
awk '{if($1=="chr2L" || $1=="chr2R" || $1=="chr3L" || $1=="chr3R" || $1=="chr4" || $1=="chrX") print}' |\
sort -k1,1 -k2,2n > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dm6.bed

### step 17-5: remove tRNA reads by intersection
bedtools intersect -v -s -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dm6.bed \
-b ${references}/Dmel/dmel-all-tRNA-r6.31.tRNA-extended.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dm6.tRNA-excluded.bed

### step 17-6: map reads to dm6, allowing up to 1MM, all best strata
index_genome="${references}/Dmel/indices/dmel-all-chromosome-r6.31.nosmallChr"
bowtie -f -v 1 --all --best --strata -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.bed

### step 17-5: remove tRNA reads by intersection
bedtools intersect -v -s -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.bed \
-b ${references}/Dmel/dmel-all-tRNA-r6.31.tRNA-extended.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.tRNA-excluded.bed

### step 17-6: extract g1g9, g2g10_revComp and last9 sequences from piRNA genome unique mappers
mkdir -p ${analysis}/${lib_sRNA}/jellyfish
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr($4,1,9)
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dm6.tRNA-excluded.bed > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_unique_g1g9.fasta
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr($4,2,9)
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dm6.tRNA-excluded.bed | fastx_reverse_complement - > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_unique_g2g10_revComp.fasta
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr(a[1],length(a[1])-8,9)
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dm6.tRNA-excluded.bed > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_unique_last9.fasta

### step 17-7: extract g1g9, g2g10_revComp and last9 sequences from piRNA genome all mappers
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.tRNA-excluded.bed |\
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr($4,1,9)}'  > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_all_g1g9.fasta
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.tRNA-excluded.bed |\
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr($4,2,9)}' | fastx_reverse_complement - > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_all_g2g10_revComp.fasta
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dm6.all-best-strata.tRNA-excluded.bed |\
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr(a[1],length(a[1])-8,9)}' > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_all_last9.fasta

### step 17-8: count the occurrences of 9mers by jellyfish
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

### step 17-8-9: measure linkage using R
for CLASS in unique all; do
Rscript /scratch/lf10/rh1772/fly_piRNAs/scripts/noYb-paper/jellyfish_linkage_Dspp_mouse.R DIRECTORY="${analysis}" LIB=${lib_sRNA} CLASS=${CLASS}
done
### step 17: analysing piRNAs from the ovarian somatic cells (OSCs) for the in-trans pingpong analysis --- END ---


### step 18: collecting the in-trans pingpong linkage values across libraries
for CLASS in unique all; do
for lib_sRNA in w1118_ovary_oxidised_sRNA Dpse_ovary_oxidised_sRNA Deug_ovary_oxidised_sRNA_rep1 Deug_ovary_oxidised_sRNA_rep2 SRR1746887 SRR1104823; do
awk -v LIB=${lib_sRNA} -v CLASS=${CLASS} '{if(NR>1) print $1,$2,LIB"@"CLASS
}' ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_all_m9.linkage.txt >> ${analysis}/in-trans-pingpong.stats
done
done
### in-trans-pingpong_plots.R was used to make a barchart in Figure 6 and Extended Data Figure 9.

### post-processing of the analysed data for small RNA libraries --- END ---

### processing and analysing the small RNA sequencing libraries --- END ---
