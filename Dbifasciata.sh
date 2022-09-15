### the code was used to analyse the sequencing data for Drosophila bifasciata.

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
lib_pA="<name of the poly-A+ RNA seq library>"
lib_ChIP="<name of the ChIP seq library>"
mkdir -p ${analysis}/${lib_sRNA}
mkdir -p ${analysis}/${lib_pA}
mkdir -p ${analysis}/${lib_ChIP}


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
### step 0-2-1: download the genomic sequence of the GCA_009664405.1_UCBerk_Dbif_1.0 assembly of the Drosophila bifasciata genome
mkdir -p ${references}/Dbif/indices
wget --directory-prefix="${references}/Dbif" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/664/405/GCA_009664405.1_UCBerk_Dbif_1.0/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna.gz

### step 0-2-2: make a bowtie index
gunzip ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna.gz
bowtie-build ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna ${references}/Dbif/indices/UCBerk_Dbif_1.0

### step 0-2-3: generate the size file
fasta_formatter -i ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna -t |\
awk '{print $1,length($NF)}' | tr ' ' '\t' > ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.sizes


### step 0-3: generating a bed file for extended genomic regions of tRNAs
### bowtie mapping of the tRNA sequences from D.pseudoobscura (see Dpseudoobscura.sh) to the UCBerk_Dbif_1.0 genome, and extend the mappings by 100nt both upstream and downstream
index_genome="${references}/Dbif/indices/UCBerk_Dbif_1.0"
bowtie -f -v 1 -a -S ${index_genome} ${references}/Dpse/GCF_009870125.1_UCI_Dpse_MV25_genomic.tRNAs.fasta |\
samtools view -bS - | bamToBed -i - | sort -k1,1 -k2,2n | awk '!seen[$1" "$2" "$3" "$6]++' |\
awk '{print $1,$2-100,$3+100,$4,$5,$6}' | awk '{if($2<0) {$2=0; print} else print}' |\
tr ' ' '\t' > ${references}/Dbif/UCBerk_Dbif_1.0.tRNAs-extended.bed


### step 0-4: generating a bed tile for genome unique 0.5kb tiles from the GCF_018153835.1_ASM1815383v1 genome assembly of Drosophila eugracilis
### step 0-4-1: obtain every 25mer from the genome, only sense
### convert lowercase to uppercase first
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna |\
fasta_formatter -  -t | awk '{for(i=0;i<=length($NF)-25;i++) print ">"$1"@"i"_+:"substr($NF,i+1,25)"\n"substr($NF,i+1,25)
}' > ${references}/Dbif/Dbif_UCBerk_25mers.plus.fasta

### step 0-4-2: obtain genome unique mappers
index_genome="${references}/Dbif/indices/UCBerk_Dbif_1.0"
bowtie -p 8 -f -v 0 -m 1 -S ${index_genome} ${references}/Dbif/Dbif_UCBerk_25mers.plus.fasta |\
samtools view -@ 8 -bS - | bamToBed -i - > ${references}/Dbif/Dbif_UCBerk_25mers.plus.unique-mappers.bed

### step 0-4-3: count the unique mappers in 0.5kb bins
### take if the centre of 25mer sits between the coordinate of 1 to 500 as the first 0.5kb bin
### take as the second 0.5kb bin if it sits between 501 to 1000
awk '{split($4,a,"@"); split(a[2],b,"_"); TILE[$1"_"int((b[1]+12)/500)]++
} END {for(var in TILE) print var,TILE[var]}' ${references}/Dbif/Dbif_UCBerk_25mers.plus.unique-mappers.bed > ${references}/Dbif/Dbif_UCBerk_0.5kbtiles.25mers.counts

### step 0-4-4: make a bed file that contains all kb tiles that are gt 85% mappability
awk '{split($1,a,"_"); if($2>425) print a[1],a[2]*500,(a[2]+1)*500,a[1]":"a[2]*500"-"(a[2]+1)*500"_+ . +""\n"a[1],a[2]*500,(a[2]+1)*500,a[1]":"a[2]*500"-"(a[2]+1)*500"_- . -"
}' ${references}/Dbif/Dbif_UCBerk_0.5kbtiles.25mers.counts |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/Dbif/Dbif_UCBerk_85percent_0.5kbtiles.bed


### step 0-5: preparing for the 0.2kb genomic tile analyses of the all genome mappers --- START ---
### step 0-5-1: make a bed file of 0.2kb genomic tiles
awk '{for(i=0;i<=int($2/200)-1;i++) {print $1,i*200,(i+1)*200,$1":"i*200"-"(i+1)*200":- . -""\n"$1,i*200,(i+1)*200,$1":"i*200"-"(i+1)*200":+ . +"
} {print $1,int($2/200)*200,$2,$1":"int($2/200)*200"-"$2":- . -""\n"$1,int($2/200)*200,$2,$1":"int($2/200)*200"-"$2":+ . +"}
}' ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.sizes |\
tr ' ' '\t' > ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.200nt.window.bed

### step 0-5-2: add annotations of piRNA clusters, gypsy and TE insertions, and exons
mkdir -p ${analysis}/Dbif
### step 0-5-2-1: piRNA clusters
bedtools intersect -wo -s -a ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.200nt.window.bed \
-b ${references}/Dbif/Dbif_UCBerk_piRNA-clusters.bed |\
awk '{print $4,$NF/200,$10}' > ${analysis}/Dbif/dbif_0.2kb_tiles.counts
### ${references}/Dbif/Dbif_UCBerk_piRNA-clusters.bed looks as follows:
CM019041.1      13765399        13803237        CM1_137_uni     .       -
CM019041.1      13273939        13557426        CM1_132-135_dual        .       -
CM019041.1      13273939        13557426        CM1_132-135_dual        .       +
CM019041.1      9951953 9996737 CM1_99_uni      .       -

### step 0-5-2-2: Run RepeatMasker 4.1.0 by specifying "drosophila"
fastafile="${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna"
RepeatMasker -species "drosophila" -parallel 10 -xsmall ${fastafile}

### step 0-5-2-3: extract gypsy insertions --- START ---
cat ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna.out | grep -i "gypsy" |\
awk '{if($9=="+") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9; else if($9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"
}' | tr ' ' '\t' > ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna.out.gypsy.bed

### extract sense and antisense gypsy insertions, these bed files are used in the IGV browser tracks
cat ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna.out | grep -i "gypsy" |\
awk '{if($9=="+" ) print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9}' |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna.out.gypsy.plus.bed
cat ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna.out | grep -i "gypsy" |\
awk '{if($9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"}' |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna.out.gypsy.minus.bed

### annotate tiles that overlap with gypsy insertions at the sense orientation
bedtools intersect -f 0.5 -s -a ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.200nt.window.bed \
-b ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna.out.gypsy.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"gypsy_S"}' >> ${analysis}/Dbif/dbif_0.2kb_tiles.counts

### annotate tiles that overlap with gypsy insertions at the antisense orientation
bedtools intersect -f 0.5 -S -a ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.200nt.window.bed \
-b ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna.out.gypsy.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"gypsy_AS"}' >> ${analysis}/Dbif/dbif_0.2kb_tiles.counts
### step 0-5-2-3: extract gypsy insertions --- END ---

### step 0-5-2-4: extract TE insertions --- START ---
### remove "Simple_repeat", "Low_complexity", "rRNA" and "tRNA"
awk '{if($11 != "Simple_repeat" && $11 != "Low_complexity" && $11 != "rRNA" && $11 != "tRNA") print}' ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna.out |\
awk '{if($9=="+") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9; else if($9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"
}' | tr ' ' '\t' > ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna.out.allTE.bed

### annotate tiles that overlap with TE insertions at the sense orientation
bedtools intersect -f 0.5 -s -a ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.200nt.window.bed \
-b ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna.out.allTE.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"TE_S"}' >> ${analysis}/Dbif/dbif_0.2kb_tiles.counts

### annotate tiles that overlap with TE insertions at the antisense orientation
bedtools intersect -f 0.5 -S -a ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.200nt.window.bed \
-b ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna.out.allTE.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"TE_AS"}' >> ${analysis}/Dbif/dbif_0.2kb_tiles.counts
### step 0-5-2-4: extract TE insertions --- END ---

### step 0-5-2-5: add RNAseq reads
### see the part of processing polyA + RNAseq from line 406
lib_pA="Dbif_ovary_polyA"
READS=`(awk '{if($0~"Uniquely mapped reads number") print $NF}' ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_Log.final.out)`
bedtools intersect -s -F 0.5 -c -a ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.200nt.window.bed \
-b ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR.sorted.bed |\
awk -v READS=${READS} -v LIB=${lib_pA} '{print $4,$7/READS*1000000,LIB}' >> ${analysis}/Dbif/dbif_0.2kb_tiles.counts
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
### run bowtie to map reads to the genome, allowing up to 1MM, unique mappers
index_genome="${references}/Dbif/indices/UCBerk_Dbif_1.0"
bowtie -f -v 1 -m 1 -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - | sort -k1,1 -k2,2n > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCBerk_Dbif_1.0.bed


### step 4: remove tRNA reads by intersection
### exclude reads that intersect with tRNA annotations from the genome unique mappers
bedtools intersect -v -s -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCBerk_Dbif_1.0.bed \
-b ${references}/Dbif/UCBerk_Dbif_1.0.tRNAs-extended.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCBerk_Dbif_1.0.tRNA-excluded.bed


### step 5: generate the bigwig files for the genome unique mappers --- START ---
### step 5-1: make bedgraph files
### only count reads that are greater than 22nt
### for IP libraries and oxidised libraries --- START ---
### normalise read counts to one million genome unique mappers
TOTAL=`(awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]} END {print count
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCBerk_Dbif_1.0.tRNA-excluded.bed)`
### bedgraph output per 1Mio miRNA mappers (plus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCBerk_Dbif_1.0.tRNA-excluded.bed |\
bedtools genomecov -strand + -split -bg -i - -g ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.sizes |\
awk -v TOTAL=${TOTAL} '{print $1,$2,$3,$4/TOTAL*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_plus.UCBerk_Dbif_1.0.tRNA-excluded.bg
### bedgraph output per 1Mio miRNA mappers (minus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCBerk_Dbif_1.0.tRNA-excluded.bed |\
bedtools genomecov -strand - -split -bg -i - -g ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.sizes |\
awk -v TOTAL=${TOTAL} '{print $1,$2,$3,-$4/TOTAL*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_minus.UCBerk_Dbif_1.0.tRNA-excluded.bg
### for IP libraries and oxidised libraries --- END ---

### step 5-2: make bigwig files
for strand in plus minus; do
bedGraphToBigWig ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_${strand}.UCBerk_Dbif_1.0.tRNA-excluded.bg \
${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.sizes ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_${strand}.UCBerk_Dbif_1.0.tRNA-excluded.bw
done
### step 5: generate the bigwig files for the genome unique mappers --- END ---


### step 6: 0.5kb tile analysis for the genome unique mappers
### count piRNA reads (greater than 22nt) that uniquely mapped to the 0.5kb tiles
READS=`(awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]} END {print count
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCBerk_Dbif_1.0.tRNA-excluded.bed)`
bedtools intersect -wo -s -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCBerk_Dbif_1.0.tRNA-excluded.bed \
-b ${references}/Dbif/Dbif_UCBerk_85percent_0.5kbtiles.bed |\
awk -v READS=${READS} -v LIB=${lib_sRNA}  '{split($4,a,"@"); if(length(a[1])>22) TILE[$10]+=a[2]} END {for(var in TILE) print var,TILE[var]/READS*1000000,LIB
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCBerk_Dbif_1.0.tRNA-excluded.gt22nt.0.5kbtiles.counts



### step 7: 0.2kb tile analysis for the genome all mappers --- START ---
### step 7-1: map reads to the UCBerk_Dbif_1.0 assembly of the D. bifasciata genome, allowing up to 1MM, all best strata
index_genome="${references}/Dbif/indices/UCBerk_Dbif_1.0"
bowtie -f -v 1 --all --best --strata -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.UCBerk_Dbif_1.0.all-best-strata.bed

### step 7-2: add mapping counts to the fasta entries. Mapping counts will be used to distribute the read coverage across tiles.
awk '{READS[$4]++} END {for(var in READS) {split(var,a,"@"); print ">"var":"READS[var]"\n"a[1]}
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.UCBerk_Dbif_1.0.all-best-strata.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.mapping-counts.fa

### step 7-3: map again all best strata
index_genome="${references}/Dbif/indices/UCBerk_Dbif_1.0"
bowtie -f -v 1 --all --best --strata -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.mapping-counts.fa |\
samtools view -bS - | bamToBed -i - | sort -k1,1 -k2,2n > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.UCBerk_Dbif_1.0.all-best-strata.mapping-counts.bed

### step 7-4: exclude reads that intersect with tRNA annotations
bedtools intersect -v -s -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.UCBerk_Dbif_1.0.all-best-strata.mapping-counts.bed \
-b ${references}/Dbif/UCBerk_Dbif_1.0.tRNAs-extended.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.UCBerk_Dbif_1.0.all-best-strata.mapping-counts.tRNA-excluded.bed

### step 7-5: collect read counts per 200nt windows from all-best-strata
TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.UCBerk_Dbif_1.0.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if($3-$2>22) count+=b[1]} END {print count}')`
awk 'BEGIN{OFS="\t"} {if($3-$2 > 22) print}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.UCBerk_Dbif_1.0.all-best-strata.mapping-counts.tRNA-excluded.bed |\
bedtools intersect -s -wo -F 0.5 -a ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.200nt.window.bed -b - |\
awk -v TOTAL=${TOTAL} -v LIB=${lib_sRNA} '{split($10,a,"@"); split(a[2],b,":"); TILE[$4]+=b[1]/b[2]} END {for(var in TILE) print var,TILE[var]/TOTAL*1000000,LIB"_abs"
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.UCBerk_Dbif_1.0.all-best-strata.mapping-counts.tRNA-excluded.gt22nt.0.2kbtiles.counts
### step 7: 0.2kb tile analysis for the genome all mappers --- END ---

### common analyses for all small RNA libraries --- END ---


### post-processing of the analysed data for small RNA libraries --- START ---

### step 8: 0.5kb tile analysis of the genome unique mappers
mkdir -p ${analysis}/Dbif
### step 8-1: add annotations of piRNA cluster tiles to the tiles
bedtools intersect -wo -s -a ${references}/Dbif/Dbif_UCBerk_85percent_0.5kbtiles.bed \
-b ${references}/Dbif/Dbif_UCBerk_piRNA-clusters.bed |\
awk '{print $4,$NF/500,$10}' > ${analysis}/Dbif/dbif_0.5kb_tiles.counts

### step 8-2: add counts from small RNA libraries
for lib_sRNA in <oxidied whole ovary small RNA library> <oxidied embryonic small RNA library>; do
#for lib_sRNA in Dbif_ovary_oxidised_sRNA Dbif_embryo_oxidised_sRNA; do
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.UCBerk_Dbif_1.0.tRNA-excluded.gt22nt.0.5kbtiles.counts >> ${analysis}/Dbif/dbif_0.5kb_tiles.counts
done
### step 8-3: tile_analysis_plots.R was used to make a scatter plot in Extended Data Figure 6.


### step 9: 0.2kb tile analysis of the genome all mappers --- START ---
### step 9-1: add counts from small RNA libraries
for lib_sRNA in <oxidied whole ovary small RNA library> <oxidied embryonic small RNA library>; do
#for lib_sRNA in Dbif_ovary_oxidised_sRNA Dbif_embryo_oxidised_sRNA; do
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.UCBerk_Dbif_1.0.all-best-strata.mapping-counts.tRNA-excluded.gt22nt.0.2kbtiles.counts >> ${analysis}/Dbif/dbif_0.2kb_tiles.counts
done

### step 9-2: make a table using reshape2 in R and count piRNA reads from individual annotations
#R/4.0.0
library(reshape2)
table <- read.table("${analysis}/Dbif/dbif_0.2kb_tiles.counts", header=F)
table_d=dcast(table,V1~V3,value.var="V2")
table_d[is.na(table_d)] <- 0
write.table(table_d, file="${analysis}/Dbif/dbif_0.2kb_tiles.counts.table", quote=F, col.names=T, row.names = F)
colnames(table_d)

[1] "V1"                        "CM1_132-135_dual"
[3] "CM1_137_uni"               "CM1_99_uni"
[5] "D-bifasciata-pA"           "Dbif_embryo_oxidised_sRNA_abs"
[7] "Dbif_ovary_oxidised_sRNA_abs" "gypsy_AS"
[9] "gypsy_S"                   "TE_AS"
[11] "TE_S"
#R/4.0.0

### count total somatic reads
soma_total=`(awk '{if(NR>1 && $7>0) print}' ${analysis}/Dbif/dbif_0.2kb_tiles.counts.table |\
awk '{if($6/$7 < 0.1) print}' | awk '{count+=$7} END {print count}')`

### count somatic uni-stranded cluster reads
uni_total=`(awk '{if(NR>1 && $7>0) print}' ${analysis}/Dbif/dbif_0.2kb_tiles.counts.table |\
awk '{if($6/$7 < 0.1) print}' | awk '{if($3>0 || $4>0) count+=$7} END {print count}')`

### count exonic piRNA reads
### consider tiles of greater than 1 RPKM
### exclude tiles that overlap with uni clusters and TE annotations
READS=`(awk '{if(NR>1) pA+=$5} END {print pA}' ${analysis}/Dbif/dbif_0.2kb_tiles.counts.table)`
pA=`(awk '{if(NR>1 && $7>0) print}' ${analysis}/Dbif/dbif_0.2kb_tiles.counts.table |\
awk '{if($6/$7 < 0.1) print}' | awk -v READS=$READS '{if($3==0 && $4==0 && $10==0 && $11==0 && $5>0.2*READS/1000000) count+=$7} END {print count}')`

### count somatic gypsy_AS reads outside the clusters
no_uni_gypsy_AS=`(awk '{if(NR>1 && $7>0) print}' ${analysis}/Dbif/dbif_0.2kb_tiles.counts.table |\
awk '{if($6/$7 < 0.1) print}' | awk '{if($3==0 && $4==0 && $8>0 && $9==0) count+=$7} END {print count}')`

### count somatic gypsy_S reads outside the clusters
no_uni_gypsy_S=`(awk '{if(NR>1 && $7>0) print}' ${analysis}/Dbif/dbif_0.2kb_tiles.counts.table |\
awk '{if($6/$7 < 0.1) print}' | awk '{if($3==0 && $4==0 && $9>0 && $8==0) count+=$7} END {print count}')`

### collect the stats
others=`echo "$soma_total - $uni_total - $pA - $no_uni_gypsy_AS - $no_uni_gypsy_S" | bc`
printf $soma_total" Dbif soma_total""\n"$uni_total" Dbif uni_total""\n"$pA" Dbif exons""\n"$no_uni_gypsy_AS" Dbif no_uni_gypsy_AS""\n"$no_uni_gypsy_S" Dbif no_uni_gypsy_S""\n"$others" Dbif others""\n" >> ${analysis}/0.2kb_tiles.counts.stats
### step 9: 0.2kb tile analysis of the genome all mappers --- END ---
### post-processing of the analysed data for small RNA libraries --- END ---

### processing and analysing the small RNA sequencing libraries --- END ---



### annotating the CDS in the D.bifasciata genome based on tBlastn search of D.pseudoobscura CDS sequences --- START ---
### step 1: downloading annotated CDS sequences from the D.pseudoobscura genome
wget --directory-prefix="${references}/Dpse" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/870/125/GCF_009870125.1_UCI_Dpse_MV25/GCF_009870125.1_UCI_Dpse_MV25_protein.faa.gz
### input has 26574 protein sequences from the UCI_Dpse_MV25

### step 2: make a blast database for the UCBerk_Dbif_1.0 assembly
fastafile="${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna"
makeblastdb -in ${fastafile} -parse_seqids -dbtype nucl -title "UCBerk_Dbif_1.0" \
-out ${references}/Dbif/Dbif_UCBerk_blast

### step 3: run tblastn to predict CDS exons in the UCBerk_Dbif_1.0 assembly
### we only included proteins smaller than 4000aa (26266 entries). Searching for the homologs of larger proteins required too much memory.
### apply an evalue cutoff of 1 x e-10
input_fasta="${references}/Dpse/GCF_009870125.1_UCI_Dpse_MV25_protein.st4000aa.faa"
output="${references}/Dbif/GCF_009870125.1_UCI_Dpse_MV25.CDS.Dbif_UCBerk_tblast.st4000aa.csv"
blast_database="${references}/Dbif/Dbif_UCBerk_blast"
tblastn -db ${blast_database} -outfmt 10 -query ${input_fasta} -out ${output} -evalue 1e-10

### step 4: output a bed file and gtf file
cat ${output} | tr ',' ' ' | awk '{if($9>$10) print $2,$10,$9,$1,". -"; else if ($10>$9) print $2,$9,$10,$1,". +"}' |\
tr ' ' '\t' | sort -k1,1 -k2,2n > ${references}/Dbif/GCF_009870125.1_UCI_Dpse_MV25.CDS.Dbif_UCBerk_tblast.st4000aa.bed

### csv file looks as follows:
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
### blast hits use 1 based coordinates
awk 'BEGIN{OFS="\t"} {split($1,a,","); if(a[9]>a[10]) print a[2],"UCBerk_Dbif_1.0","exon",a[10],a[9],".","-","0","transcript_id \""a[1]"\";" ;
else if(a[10]>a[9]) print a[2],"UCBerk_Dbif_1.0","exon",a[9],a[10],".","+","0","transcript_id \""a[1]"\";";
}' ${output} | sort -k1,1 -k4,4n > ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0.CDS.tblastn_Dpse.gtf
### this gtf file is used to display the putative CDSs in the genome browser tracks
### annotating the CDS in the D.bifasciata genome based on tBlastn search of D.pseudoobscura CDS sequences --- END ---



### processing and analysing the polyA+ RNA sequencing libraries --- START ---

### below is the RNA library processed in this part of the script:
Dbif_ovary_polyA

### RNA seq libraries
5’ AATGATACGGCGACCACCGAGATCTACAC[TCTTTCCCTACACGACGCTCTTCCGATCT]--- [R1 ->] --- [<- R2] ---[AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC] ---NNNNNN--- ATCTCGTATGCCGTCTTCTGCTTG 3’
R1: reverse complement
R2: sense

### step 1: trim adapters, filter poor-quality reads, and take only the reads that have mates
lib_pA="Dbif_ovary_polyA"
mkdir -p ${analysis}/${lib_pA}
fastx_clipper -l 35 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-i <(zcat ${raw_data}/${lib_pA}_L1_1.fq.gz | fastq_quality_filter -q 33 -p 40) > ${analysis}/${lib_pA}/${lib_pA}_R1.trimmed.fastq
fastx_clipper -l 35 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGA \
-i <(zcat ${raw_data}/${lib_pA}_L1_2.fq.gz | fastq_quality_filter -q 33 -p 40) > ${analysis}/${lib_pA}/${lib_pA}_R2.trimmed.fastq

### pair mates using fastq-pair and take only the reads that have mates
fastq_pair ${analysis}/${lib_pA}/${lib_pA}_R1.trimmed.fastq ${analysis}/${lib_pA}/${lib_pA}_R2.trimmed.fastq
for TYPE in R1 R2; do
gzip ${analysis}/${lib_pA}/${lib_pA}_${TYPE}.trimmed.fastq.paired.fq
rm ${analysis}/${lib_pA}/${lib_pA}_${TYPE}.trimmed.fastq.single.fq
rm ${analysis}/${lib_pA}/${lib_pA}_${TYPE}.trimmed.fastq
done

### step 2: make star index using the CDS annotation see above lines
mkdir -p ${references}/Dbif/indices/UCBerk_Dbif_1.0_STAR
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ${references}/Dbif/indices/UCBerk_Dbif_1.0_STAR \
--genomeSAindexNbases 12 \
--genomeFastaFiles ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.fna \
--sjdbGTFfile ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0.CDS.tblastn_Dpse.gtf \
--sjdbOverhang 100

### step 3: STAR mapping
### jM option to report split reads
mkdir -p ${analysis}/${lib_pA}/STAR
STAR --runThreadN 8 --genomeDir ${references}/Dbif/indices/UCBerk_Dbif_1.0_STAR \
--outSAMattributes NH HI AS nM MD jM \
--outFilterMismatchNmax 3 \
--outSJfilterReads All \
--outSAMreadID Number \
--readFilesCommand zcat \
--outFileNamePrefix ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_ \
--readFilesIn ${analysis}/${lib_pA}/${lib_pA}_R1.trimmed.fastq.paired.fq.gz ${analysis}/${lib_pA}/${lib_pA}_R2.trimmed.fastq.paired.fq.gz

### step 4: generate bigwig files for the IGV browser --- START ---
### take unique mappers and normalise the read counts
READS=`(awk '{if($0~"Uniquely mapped reads number") print $NF}' ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_Log.final.out)`

### sam to bed
samtools view -q 10 -bS ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_Aligned.out.sam |\
bamToBed -bed12 -i - > ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR.bed

### reverse the strandness of /1 and sort bed files
awk '{if($4~"/1" && $6=="-") {$6="+"; print}
else if($4~"/1" && $6=="+") {$6="-"; print}
else print}' ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR.bed |\
sort -k1,1 -k2,2 | tr ' ' '\t' > ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR.sorted.bed

### bedgraph output per 1Mio unique mappers (plus strand)
bedtools genomecov -strand + -split -bg -i ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR.sorted.bed \
-g ${references}/Dbif/indices/UCBerk_Dbif_1.0_STAR/chrNameLength.txt |\
awk -v READS=${READS} '{print $1,$2,$3,$4/READS*1000000}' | tr ' ' '\t' > ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_plus.bg
### bedgraph output per 1Mio unique mappers (minus strand)
bedtools genomecov -strand - -split -bg -i ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR.sorted.bed \
-g ${references}/Dbif/indices/UCBerk_Dbif_1.0_STAR/chrNameLength.txt |\
awk -v READS=${READS} '{print $1,$2,$3,-$4/READS*1000000}' | tr ' ' '\t' > ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_minus.bg

### generate bigWig files
export PATH="/home/150/rh1772/bin/x86_64:$PATH"
for strand in plus minus; do
bedGraphToBigWig ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_${strand}.bg \
${references}/Dbif/indices/UCBerk_Dbif_1.0_STAR/chrNameLength.txt ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_${strand}.bw
done
### step 4: generate bigwig files for the IGV browser --- END ---

### step 5: extract reads that mapped near the 5' end of the piRNA cluster CM1_137_uni --- START ---
### intersect bam
samtools view -@ 8 -q 10 -bS ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_Aligned.out.sam |\
samtools sort -@ 8 > ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_Aligned.out.bam

### print out the coordinate of the cluster CM1_137_uni 5' end
printf "CM019041.1 13802267 13808764 CM1_137_uni_5end . -""\n" | tr ' ' '\t' |\
# -wa: Write the original entry in A for each overlap.
# this prints out the sam header in the intersected file as well.
bedtools intersect -wa -a ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_Aligned.out.bam -b - |\
samtools sort > ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_Aligned.out.CM1_137_uni_5end.bam
samtools index ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_Aligned.out.CM1_137_uni_5end.bam
### step 5: extract reads that mapped near the 5' end of the piRNA cluster CM1_137_uni --- END ---

### processing and analysing the polyA+ RNA sequencing libraries --- END ---


### processing and analysing ChIPseq libraries --- START ---

### below is the ChIPseq libraries processed in this part of the script:
Dbif_ChIP_input
Dbif_ChIP_polII-CTD

### ChIP seq libraries
5’ AATGATACGGCGACCACCGAGATCTACAC[TCTTTCCCTACACGACGCTCTTCCGATCT]--- [R1 ->] --- [<- R2] ---[AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC] ---NNNNNN--- ATCTCGTATGCCGTCTTCTGCTTG 3’

### step 1: trim adapters and filter poor-quality reads
mkdir -p ${analysis}/${lib_ChIP}
fastx_clipper -l 35 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-i <(zcat ${raw_data}/${lib_ChIP}_R1.fastq.gz | fastq_quality_filter -q 33 -p 40) | gzip > ${analysis}/${lib_ChIP}/${lib_ChIP}_R1.trimmed.fastq.gz

### step 2: map R1 reads to the D.pseudoobscura Nanopore assembly, allowing up to 3 mismatches, unique mappers
index_genome="${references}/Dbif/indices/UCBerk_Dbif_1.0"
bowtie -q -v 3 -m 1 -S ${index_genome} ${analysis}/${lib_ChIP}/${lib_ChIP}_R1.trimmed.fastq.gz |\
samtools view -bS - | bamToBed -i - | sort -k1,1 -k2,2n > ${analysis}/${lib_ChIP}/${lib_ChIP}_genome-unique-mappers.bed

### step 3: generate bigwig files for the IGV browser --- START ---
### bedgraph output per 1Mio unique mappers (plus strand)
READS=`(cat ${analysis}/${lib_ChIP}/${lib_ChIP}_genome-unique-mappers.bed | wc -l | awk '{print $NF}')`

### bedgraph output per 1Mio unique mappers
bedtools genomecov -bg -i ${analysis}/${lib_ChIP}/${lib_ChIP}_genome-unique-mappers.bed \
-g ${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.sizes |\
awk -v READS=${READS} '{print $1,$2,$3,$4/READS*1000000}' | tr ' ' '\t' > ${analysis}/${lib_ChIP}/${lib_ChIP}_genome-unique-mappers.bg

### generate bigwig
bedGraphToBigWig ${analysis}/${lib_ChIP}/${lib_ChIP}_genome-unique-mappers.bg \
${references}/Dbif/GCA_009664405.1_UCBerk_Dbif_1.0_genomic.sizes ${analysis}/${lib_ChIP}/${lib_ChIP}_genome-unique-mappers.bw
### step 3: generate bigwig files for the IGV browser --- END ---
### processing and analysing ChIPseq libraries --- END ---
