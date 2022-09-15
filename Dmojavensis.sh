### the code was used to analyse the sequencing data for Drosophila mojavensis.

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

### define the variables
references="<directory where the genome sequence file and the gtf file are kept>"
raw_data="<directory where the fastq files are kept>"
analysis="<directory where the output files are kept>"
lib_sRNA="<name of the small RNA seq library>"
mkdir -p ${analysis}/${lib_sRNA}

### processing and analysing the small RNA sequencing libraries --- START ---

### below are the small RNA libraries processed in this part of the script:
Dbif_ovary_oxidised_sRNA
Dbif_embryo_oxidised_sRNA


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
### step 0-2-1: download the genomic sequence of the GCF_018153725.1_ASM1815372v1 assembly of the Drosophila mojavensis genome
mkdir -p ${references}/Dmoj/indices
wget --directory-prefix="${references}/Dmoj" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/153/725/GCF_018153725.1_ASM1815372v1/GCF_018153725.1_ASM1815372v1_genomic.fna.gz
wget --directory-prefix="${references}/Dmoj" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/153/725/GCF_018153725.1_ASM1815372v1/GCF_018153725.1_ASM1815372v1_genomic.gtf.gz

### step 0-2-2: make a bowtie index
gunzip ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.fna.gz
bowtie-build ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.fna ${references}/Dmoj/indices/Dmoj_UStan_RefSeq

### step 0-2-3: generate the size file
fasta_formatter -i ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.fna -t |\
awk '{print $1,length($NF)}' | tr ' ' '\t' > ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.sizes


### step 0-3: generating a bed file for extended genomic regions of tRNAs
### fetch tRNA bed coordinates from the Dmoj_UStan_RefSeq assembly
zcat ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.gtf.gz | grep 'gbkey "tRNA"'  |\
awk -F";" '{print $1}' | tr -d '"' | awk '{print $1,$4-101,$5+100,"Dmoj_"$NF,".",$7}' | awk '{if($2<0) {$2=0; print} else print}' |\
tr ' ' '\t' | sort -k1,1 -k2,2n > ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.tRNAs-extended.bed


### step 0-4: generating a bed tile for genome unique 0.5kb tiles from the Dmoj_UStan_RefSeq assembly of the D. mojavensis genome
### step 0-4-1: obtain every 25mer from the genome, only sense
### convert lowercase to uppercase first
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.fna |\
fasta_formatter -  -t | awk '{for(i=0;i<=length($NF)-25;i++) print ">"$1"@"i"_+:"substr($NF,i+1,25)"\n"substr($NF,i+1,25)
}' > ${references}/Dmoj/Dmoj_UStan_25mers.plus.fasta

### step 0-4-2: obtain genome unique mappers
index_genome="${references}/Dmoj/indices/Dmoj_UStan_RefSeq"
bowtie -p 8 -f -v 0 -m 1 -S ${index_genome} ${references}/Dmoj/Dmoj_UStan_25mers.plus.fasta |\
samtools view -@ 8 -bS - | bamToBed -i - > ${references}/Dmoj/Dmoj_UStan_25mers.plus.unique-mappers.bed

### step 0-4-3: count the unique mappers in 0.5kb bins
### take if the centre of 25mer sits between the coordinate of 1 to 500 as the first 0.5kb bin
### take as the second 0.5kb bin if it sits between 501 to 1000
awk '{split($4,a,"@"); split(a[2],b,"_"); TILE[$1"_"int((b[1]+12)/500)]++
} END {for(var in TILE) print var,TILE[var]}' ${references}/Dmoj/Dmoj_UStan_25mers.plus.unique-mappers.bed > ${references}/Dmoj/Dmoj_UStan_0.5kbtiles.25mers.counts

### step 0-4-4: make a bed file that contains all kb tiles that are gt 85% mappability
awk '{split($1,a,"_"); if($2>425) print a[1]"_"a[2],a[3]*500,(a[3]+1)*500,a[1]"_"a[2]":"a[3]*500"-"(a[3]+1)*500"_+ . +""\n"a[1]"_"a[2],a[3]*500,(a[3]+1)*500,a[1]"_"a[2]":"a[3]*500"-"(a[3]+1)*500"_- . -"
}' ${references}/Dmoj/Dmoj_UStan_0.5kbtiles.25mers.counts |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/Dmoj/Dmoj_UStan_85percent_0.5kbtiles.bed


### step 0-5: preparing for the 0.2kb genomic tile analyses of the all genome mappers --- START ---
### step 0-5-1: make a bed file of 0.2kb genomic tiles
awk '{for(i=0;i<=int($2/200)-1;i++) {print $1,i*200,(i+1)*200,$1":"i*200"-"(i+1)*200":- . -""\n"$1,i*200,(i+1)*200,$1":"i*200"-"(i+1)*200":+ . +"
} {print $1,int($2/200)*200,$2,$1":"int($2/200)*200"-"$2":- . -""\n"$1,int($2/200)*200,$2,$1":"int($2/200)*200"-"$2":+ . +"}
}' ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.sizes |\
tr ' ' '\t' > ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.200nt.window.bed


### step 0-5-2: add annotations of piRNA clusters, gypsy and TE insertions, and exons
mkdir -p ${analysis}/Dmoj
### step 0-5-2-1: piRNA clusters
bedtools intersect -wo -s -a ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.200nt.window.bed \
-b ${references}/Dmoj/Dmoj_UStan_piRNA-clusters.bed |\
awk '{print $4,$NF/200,$10}' > ${analysis}/Dmoj/Dmoj_0.2kb_tiles.counts
### ${references}/Dmoj/Dmoj_UStan_piRNA-clusters.bed looks as follows:
NW_025318667.1  24883386        24908540        uni_667 .       -

### step 0-5-2-2: Run RepeatMasker 4.1.0 by specifying "drosophila"
fastafile="${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.fna"
RepeatMasker -species "drosophila" -parallel 10 -xsmall ${fastafile}

### step 0-5-2-3: extract gypsy insertions --- START ---
cat ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.fna.out | grep -i "gypsy" |\
awk '{if($9=="+") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9; else if($9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"
}' | tr ' ' '\t' > ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.fna.out.gypsy.bed

### extract sense and antisense gypsy insertions, these bed files are used in the IGV browser tracks
cat ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.fna.out | grep -i "gypsy" |\
awk '{if($9=="+" ) print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9}' |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.fna.out.gypsy.plus.bed
cat ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.fna.out | grep -i "gypsy" |\
awk '{if($9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"}' |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.fna.out.gypsy.minus.bed

### annotate tiles that overlap with gypsy insertions at the sense orientation
bedtools intersect -f 0.5 -s -a ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.200nt.window.bed \
-b ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.fna.out.gypsy.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"gypsy_S"}' >> ${analysis}/Dmoj/Dmoj_0.2kb_tiles.counts

### annotate tiles that overlap with gypsy insertions at the antisense orientation
bedtools intersect -f 0.5 -S -a ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.200nt.window.bed \
-b ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.fna.out.gypsy.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"gypsy_AS"}' >> ${analysis}/Dmoj/Dmoj_0.2kb_tiles.counts
### step 0-5-2-3: extract gypsy insertions --- END ---

### step 0-5-2-4: extract TE insertions --- START ---
### remove "Simple_repeat", "Low_complexity", "rRNA" and "tRNA"
awk '{if($11 != "Simple_repeat" && $11 != "Low_complexity" && $11 != "rRNA" && $11 != "tRNA") print}' ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.fna.out |\
awk '{if($9=="+") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9; else if($9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"
}' | tr ' ' '\t' > ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.fna.out.allTE.bed

### annotate tiles that overlap with TE insertions at the sense orientation
bedtools intersect -f 0.5 -s -a ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.200nt.window.bed \
-b ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.fna.out.allTE.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"TE_S"}' >> ${analysis}/Dmoj/Dmoj_0.2kb_tiles.counts

### annotate tiles that overlap with TE insertions at the antisense orientation
bedtools intersect -f 0.5 -S -a ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.200nt.window.bed \
-b ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.fna.out.allTE.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"TE_AS"}' >> ${analysis}/Dmoj/Dmoj_0.2kb_tiles.counts
### step 0-5-2-4: extract TE insertions --- END ---

### step 0-5-2-5: add RNAseq reads
zcat ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.gtf.gz | tr -d '";' |\
awk 'BEGIN{OFS="\t"} {if(NR>4 && $3=="exon") print $1,$4,$5,$10":"$12,".",$7}' | sort -k1,1 -k2,2n > ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.exon.bed
bedtools intersect -wo -s -a ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.200nt.window.bed \
-b ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.exon.bed |\
awk '{exon[$4]+=$NF} END {for(var in exon) print var,exon[var]/200,"exon"}' >> ${analysis}/Dmoj/Dmoj_0.2kb_tiles.counts
### step 0-5: preparing for the 0.2kb genomic tile analyses of the all genome mappers --- END ---


### step 0-6: make a bowtie index for all Drosophila autonomous transposons downloaded from the RepBase in March 2022 (3053 entries)
mkdir -p ${references}/TEs/indices
fasta_formatter -i ${references}/TEs/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.fasta -t |\
awk '{print ">"$1"\n"$NF}' > ${references}/TEs/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.concise.fasta
bowtie-build ${references}/TEs/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.concise.fasta ${references}/TEs/indices/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.concise

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
### run bowtie to map reads to the Dmoj_UStan_RefSeq assembly, allowing up to 1MM, unique mappers
index_genome="${references}/Dmoj/indices/Dmoj_UStan_RefSeq"
bowtie -f -v 1 -m 1 -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - | sort -k1,1 -k2,2n > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Dmoj_UStan_RefSeq.bed


### step 4: remove tRNA reads by intersection
### exclude reads that intersect with tRNA annotations from the genome unique mappers
bedtools intersect -v -s -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Dmoj_UStan_RefSeq.bed \
-b ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.tRNAs-extended.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Dmoj_UStan_RefSeq.tRNA-excluded.bed


### step 5: generate the bigwig files for the genome unique mappers --- START ---
### step 5-1: make bedgraph files
### only count reads that are greater than 22nt
### for IP libraries and oxidised libraries --- START ---
### normalise read counts to one million genome unique mappers
TOTAL=`(awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]} END {print count
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Dmoj_UStan_RefSeq.tRNA-excluded.bed)`
### bedgraph output per 1Mio miRNA mappers (plus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Dmoj_UStan_RefSeq.tRNA-excluded.bed |\
bedtools genomecov -strand + -split -bg -i - -g ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.sizes |\
awk -v TOTAL=${TOTAL} '{print $1,$2,$3,$4/TOTAL*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_plus.Dmoj_UStan_RefSeq.tRNA-excluded.bg
### bedgraph output per 1Mio miRNA mappers (minus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Dmoj_UStan_RefSeq.tRNA-excluded.bed |\
bedtools genomecov -strand - -split -bg -i - -g ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.sizes |\
awk -v TOTAL=${TOTAL} '{print $1,$2,$3,-$4/TOTAL*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_minus.Dmoj_UStan_RefSeq.tRNA-excluded.bg
### for IP libraries and oxidised libraries --- END ---

### step 5-3: make bigwig files
for strand in plus minus; do
bedGraphToBigWig ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_${strand}.Dmoj_UStan_RefSeq.tRNA-excluded.bg \
${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.sizes ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_${strand}.Dmoj_UStan_RefSeq.tRNA-excluded.bw
done
### step 5: generate the bigwig files for the genome unique mappers --- END ---


### step 6: 0.5kb tile analysis for the genome unique mappers
### count piRNA reads (greater than 22nt) that uniquely mapped to the 0.5kb tiles
READS=`(awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]} END {print count
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Dmoj_UStan_RefSeq.tRNA-excluded.bed)`
bedtools intersect -wo -s -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Dmoj_UStan_RefSeq.tRNA-excluded.bed \
-b ${references}/Dmoj/Dmoj_UStan_85percent_0.5kbtiles.bed |\
awk -v READS=${READS} -v LIB=${lib_sRNA}  '{split($4,a,"@"); if(length(a[1])>22) TILE[$10]+=a[2]} END {for(var in TILE) print var,TILE[var]/READS*1000000,LIB
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Dmoj_UStan_RefSeq.tRNA-excluded.gt22nt.0.5kbtiles.counts


### step 7: 0.2kb tile analysis for the genome all mappers --- START ---
### step 7-1: map reads to the Dmoj_UStan_RefSeq assembly of the D. mojavensis genome, allowing up to 1MM, all best strata
index_genome="${references}/Dmoj/indices/Dmoj_UStan_RefSeq"
bowtie -f -v 1 --all --best --strata -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Dmoj_UStan_RefSeq.all-best-strata.bed

### step 7-2: add mapping counts to the fasta entries. Mapping counts will be used to distribute the read coverage across tiles.
awk '{READS[$4]++} END {for(var in READS) {split(var,a,"@"); print ">"var":"READS[var]"\n"a[1]}
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Dmoj_UStan_RefSeq.all-best-strata.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.mapping-counts.fa

### step 7-3: map again all best strata
index_genome="${references}/Dmoj/indices/Dmoj_UStan_RefSeq"
bowtie -f -v 1 --all --best --strata -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.mapping-counts.fa |\
samtools view -bS - | bamToBed -i - | sort -k1,1 -k2,2n > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Dmoj_UStan_RefSeq.all-best-strata.mapping-counts.bed

### step 7-4: exclude reads that intersect with tRNA annotations
bedtools intersect -v -s -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Dmoj_UStan_RefSeq.all-best-strata.mapping-counts.bed \
-b ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.tRNAs-extended.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Dmoj_UStan_RefSeq.all-best-strata.mapping-counts.tRNA-excluded.bed

### step 7-5: collect read counts per 200nt windows from all-best-strata
TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Dmoj_UStan_RefSeq.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if($3-$2>22) count+=b[1]} END {print count}')`
awk 'BEGIN{OFS="\t"} {if($3-$2 > 22) print}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Dmoj_UStan_RefSeq.all-best-strata.mapping-counts.tRNA-excluded.bed |\
bedtools intersect -s -wo -F 0.5 -a ${references}/Dmoj/GCF_018153725.1_ASM1815372v1_genomic.200nt.window.bed -b - |\
awk -v TOTAL=${TOTAL} -v LIB=${lib_sRNA} '{split($10,a,"@"); split(a[2],b,":"); TILE[$4]+=b[1]/b[2]} END {for(var in TILE) print var,TILE[var]/TOTAL*1000000,LIB"_abs"
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Dmoj_UStan_RefSeq.all-best-strata.mapping-counts.tRNA-excluded.gt22nt.0.2kbtiles.counts
### step 7: 0.2kb tile analysis for the genome all mappers --- END ---

### common analyses for all small RNA libraries --- END ---


### post-processing of the analysed data for small RNA libraries --- START ---

### step 8: 0.5kb tile analysis of the genome unique mappers
mkdir -p ${analysis}/Dmoj
### step 8-1: add annotations of piRNA cluster tiles to the tiles
bedtools intersect -wo -s -a ${references}/Dmoj/Dmoj_UStan_85percent_0.5kbtiles.bed \
-b ${references}/Dmoj/Dmoj_UStan_piRNA-clusters.bed |\
awk '{print $4,$NF/500,$10}' > ${analysis}/Dmoj/Dmoj_0.5kb_tiles.counts

### step 8-2: add counts from small RNA libraries
for lib_sRNA in <oxidied whole ovary small RNA library> <oxidied embryonic small RNA library>; do
#for lib_sRNA in Dmoj_ovary_oxidised_sRNA Dmoj_embryo_oxidised_sRNA; do
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.Dmoj_UStan_RefSeq.tRNA-excluded.gt22nt.0.5kbtiles.counts >> ${analysis}/Dmoj/Dmoj_0.5kb_tiles.counts
done
### step 8-3: tile_analysis_plots.R was used to make a scatter plot in Extended Data Figure 6.


### step 9: 0.2kb tile analysis of the genome all mappers --- START ---
### step 9-1: add counts from small RNA libraries
for lib_sRNA in <oxidied whole ovary small RNA library> <oxidied embryonic small RNA library>; do
#for lib_sRNA in Dmoj_ovary_oxidised_sRNA Dmoj_embryo_oxidised_sRNA; do
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.Dmoj_UStan_RefSeq.all-best-strata.mapping-counts.tRNA-excluded.gt22nt.0.2kbtiles.counts >> ${analysis}/Dmoj/Dmoj_0.2kb_tiles.counts
done

### step 9-2: make a table using reshape2 in R and count piRNA reads from individual annotations
#R/4.0.0
library(reshape2)
table <- read.table("${analysis}/Dmoj/Dmoj_0.2kb_tiles.counts", header=F)
table_d=dcast(table,V1~V3,value.var="V2")
table_d[is.na(table_d)] <- 0
write.table(table_d, file="${analysis}/Dmoj/Dmoj_0.2kb_tiles.counts.table", quote=F, col.names=T, row.names = F)
colnames(table_d)

[1] "V1"                       "Dmoj_embryo_oxidised_sRNA_abs"
[3] "Dmoj_ovary_oxidised_sRNA_abs" "exon"
[5] "gypsy_AS"                 "gypsy_S"
[7] "TE_AS"                    "TE_S"
[9] "uni_667"
#R/4.0.0

### count total somatic reads
soma_total=`(awk '{if(NR>1 && $3>0) print}' ${analysis}/Dmoj/Dmoj_0.2kb_tiles.counts.table |\
awk '{if($2/$3 < 0.1) print}' | awk '{count+=$3} END {print count}')`

### count somatic uni-stranded cluster reads
uni_total=`(awk '{if(NR>1 && $3>0) print}' ${analysis}/Dmoj/Dmoj_0.2kb_tiles.counts.table |\
awk '{if($2/$3 < 0.1) print}' | awk '{if($9>0) count+=$3} END {print count}')`

### count exonic piRNA reads
### exclude tiles that overlap with uni clusters and TE annotations
exons=`(awk '{if(NR>1 && $3>0) print}' ${analysis}/Dmoj/Dmoj_0.2kb_tiles.counts.table |\
awk '{if($2/$3 < 0.1) print}' | awk '{if($7==0 && $8==0 && $9==0 && $4>0) count+=$3} END {print count}')`

### count somatic gypsy_AS reads outside the clusters
no_uni_gypsy_AS=`(awk '{if(NR>1 && $3>0) print}' ${analysis}/Dmoj/Dmoj_0.2kb_tiles.counts.table |\
awk '{if($2/$3 < 0.1) print}' | awk '{if($9==0 && $5>0 && $6==0) count+=$3} END {print count}')`

### count somatic gypsy_S reads outside the clusters
no_uni_gypsy_S=`(awk '{if(NR>1 && $3>0) print}' ${analysis}/Dmoj/Dmoj_0.2kb_tiles.counts.table |\
awk '{if($2/$3 < 0.1) print}' | awk '{if($9==0 && $6>0 && $5==0) count+=$3} END {print count}')`

### collect the stats
others=`echo "$soma_total - $uni_total - $exons - $no_uni_gypsy_AS - $no_uni_gypsy_S" | bc`
printf $soma_total" Dmoj soma_total""\n"$uni_total" Dmoj uni_total""\n"$exons" Dmoj exons""\n"$no_uni_gypsy_AS" Dmoj no_uni_gypsy_AS""\n"$no_uni_gypsy_S" Dmoj no_uni_gypsy_S""\n"$others" Dmoj others""\n" >> ${analysis}/0.2kb_tiles.counts.stats
### step 9: 0.2kb tile analysis of the genome all mappers --- END ---
### post-processing of the analysed data for small RNA libraries --- END ---

### processing and analysing the small RNA sequencing libraries --- END ---
