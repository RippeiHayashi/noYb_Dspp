### the code was used to analyse the sequencing data for Drosophila pseudoobscura.

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
lib_pA="<name of the poly-A+ RNA seq library>"
${lib_ChIP}="<name of the ChIP seq library>"
mkdir -p ${analysis}/${lib_sRNA}
mkdir -p ${analysis}/${lib_pA}
mkdir -p ${analysis}/${lib_ChIP}


### processing and analysing the small RNA sequencing libraries --- START ---

### below are the small RNA libraries processed in this part of the script:
Dpse_embryo_oxidised_sRNA
Dpse_ovary_oxidised_sRNA
Dpse_ovary_unoxidised_sRNA


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
### step 0-2-1: download the genomic sequence of the Nanopore assembly of Drosophila pseudoobscura from Miller, G3, 2018, PMID: 30087105
mkdir -p ${references}/Dpse/indices
wget --directory-prefix="${references}/Dpse" https://github.com/danrdanny/Drosophila15GenomesProject/raw/master/assembledGenomes/Dpse.pass.minimap2.racon.x3.pilon.x3.fasta.gz

### step 0-2-2: simplify the contig names
zcat ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.fasta.gz | fasta_formatter - -t |\
awk '{split($1,a,"_"); print ">"a[4]"\n"$NF}' > ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta

### step 0-2-3: make a bowtie index
bowtie-build ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta ${references}/Dpse/indices/dpse-nanopore

### step 0-2-4: generate the size file
fasta_formatter -i ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta -t |\
awk '{print $1,length($2)}' | tr ' ' '\t' > ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.sizes


### step 0-3: generating a bed file for extended genomic regions of tRNAs
### step 0-3-1: download the genomic sequence and the gtf file of the UCI_Dpse_MV25 genome assembly of Drosophila pseudoobscura
wget --directory-prefix="${references}/Dpse" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/870/125/GCF_009870125.1_UCI_Dpse_MV25/GCF_009870125.1_UCI_Dpse_MV25_genomic.gtf.gz
wget --directory-prefix="${references}/Dpse" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/870/125/GCF_009870125.1_UCI_Dpse_MV25/GCF_009870125.1_UCI_Dpse_MV25_genomic.fna.gz

### step 0-3-2: fetch tRNA sequences from the UCI_Dpse_MV25 assembly
gunzip ${references}/Dpse/GCF_009870125.1_UCI_Dpse_MV25_genomic.fna.gz
zcat ${references}/Dpse/GCF_009870125.1_UCI_Dpse_MV25_genomic.gtf.gz | grep 'gbkey "tRNA"'  |\
awk -F";" '{print $1}' | tr -d '"' | awk '{print $1,$4-1,$5,"Dpse_"$NF,".",$7}' | tr ' ' '\t' | sort -k1,1 -k2,2n |\
bedtools getfasta -name -s -fi ${references}/Dpse/GCF_009870125.1_UCI_Dpse_MV25_genomic.fna -bed - > ${references}/Dpse/GCF_009870125.1_UCI_Dpse_MV25_genomic.tRNAs.fasta

### step 0-3-3: bowtie mapping of tRNA sequences and extend the mappings by 100nt both upstream and downstream
index_genome="${references}/Dpse/indices/dpse-nanopore"
bowtie -f -v 1 -a -S ${index_genome} ${references}/Dpse/GCF_009870125.1_UCI_Dpse_MV25_genomic.tRNAs.fasta |\
samtools view -bS - | bamToBed -i - | sort -k1,1 -k2,2n | awk '!seen[$1" "$2" "$3" "$6]++' |\
awk '{print $1,$2-100,$3+100,$4,$5,$6}' | awk '{if($2<0) {$2=0; print} else print}' |\
tr ' ' '\t' > ${references}/Dpse/dpse-nanopore.tRNAs-extended.bed


### step 0-4: generating a bed tile for genome unique 0.5kb tiles from the D. pseudoobscura Nanopore genome assembly
### step 0-4-1: obtain every 25mer from the genome, only sense
### convert lowercase to uppercase first
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta |\
fasta_formatter - -t | awk '{for(i=0;i<=length($NF)-25;i++) print ">"$1"@"i"_+:"substr($NF,i+1,25)"\n"substr($NF,i+1,25)
}' > ${references}/Dpse/Dpse_nanopore_25mers.plus.fasta

### step 0-4-2: obtain unique mappers to the nanopore assembly
index_genome="${references}/Dpse/indices/dpse-nanopore"
bowtie -f -v 0 -m 1 -S ${index_genome} ${references}/Dpse/Dpse_nanopore_25mers.plus.fasta |\
samtools view -@ 16 -bS - | bamToBed -i - > ${references}/Dpse/Dpse_nanopore_25mers.plus.unique-mappers.bed

### step 0-4-3: count the unique mappers in 0.5kb bins
### take if the centre of 25mer sits between the coordinate of 1 to 500 as the first 0.5kb bin
### take as the second 0.5kb bin if it sits between 501 to 1000
awk '{split($4,a,"@"); split(a[2],b,"_"); TILE[$1"_"int((b[1]+12)/500)]++
} END {for(var in TILE) print var,TILE[var]}' ${references}/Dpse/Dpse_nanopore_25mers.plus.unique-mappers.bed > ${references}/Dpse/Dpse_nanopore_0.5kbtiles.counts

### step 0-4-4: make a bed file that contains all kb tiles that are gt 85% mappability
awk '{split($1,a,"_"); if($2>425) print a[1],a[2]*500,(a[2]+1)*500,a[1]":"a[2]*500"-"(a[2]+1)*500"_+ . +""\n"a[1],a[2]*500,(a[2]+1)*500,a[1]":"a[2]*500"-"(a[2]+1)*500"_- . -"
}' ${references}/Dpse/Dpse_nanopore_0.5kbtiles.counts |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/Dpse/Dpse_nanopore_85percent_0.5kbtiles.25mers.bed


### step 0-5: preparing for the 0.2kb genomic tile analyses of the all genome mappers --- START ---
### step 0-5-1: make a bed file of 0.2kb genomic tiles
awk '{for(i=0;i<=int($2/200)-1;i++) {print $1,i*200,(i+1)*200,$1":"i*200"-"(i+1)*200":- . -""\n"$1,i*200,(i+1)*200,$1":"i*200"-"(i+1)*200":+ . +"
} {print $1,int($2/200)*200,$2,$1":"int($2/200)*200"-"$2":- . -""\n"$1,int($2/200)*200,$2,$1":"int($2/200)*200"-"$2":+ . +"}
}' ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.sizes |\
tr ' ' '\t' > ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.200nt.window.bed

### step 0-5-2: add annotations of piRNA clusters, gypsy and TE insertions, and exons
mkdir -p ${analysis}/Dpse
### step 0-5-2-1: piRNA clusters
bedtools intersect -wo -s -a ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.200nt.window.bed \
-b ${references}/Dpse/Dpse_nanopore_piRNA-clusters.bed |\
awk '{print $4,$NF/200,$10}' >> ${analysis}/Dpse/dpse_0.2kb_tiles.counts
### ${references}/Dpse/Dpse_nanopore_piRNA-clusters.bed looks as follows:
utg000003l      427114  618444  cluster3l_uni   .       +
utg000040l      337118  437667  cluster40l_dual .       +
utg000040l      337118  437667  cluster40l_dual .       -

### step 0-5-2-2: Run RepeatMasker 4.1.0 by specifying "drosophila"
fastafile="${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta"
RepeatMasker -species "drosophila" -parallel 10 -xsmall ${fastafile}

### step 0-5-2-3: extract gypsy insertions --- START ---
cat ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta.out | grep -i "gypsy" |\
awk '{if($9=="+") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9; else if($9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"
}' | tr ' ' '\t' > ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta.out.gypsy.bed

### extract sense and antisense gypsy insertions, these bed files are used in the IGV browser tracks
cat ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta.out | grep -i "gypsy" |\
awk '{if($5~"utg" && $9=="+" ) print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9}' |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta.out.gypsy.plus.bed
cat ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta.out | grep -i "gypsy" |\
awk '{if($5~"utg" && $9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"}' |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta.out.gypsy.minus.bed

### annotate tiles that overlap with gypsy insertions at the sense orientation
bedtools intersect -f 0.5 -s -a ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.200nt.window.bed \
-b ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta.out.gypsy.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"gypsy_S"}' >> ${analysis}/Dpse/dpse_0.2kb_tiles.counts

### annotate tiles that overlap with gypsy insertions at the antisense orientation
bedtools intersect -f 0.5 -S -a ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.200nt.window.bed \
-b ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta.out.gypsy.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"gypsy_AS"}' >> ${analysis}/Dpse/dpse_0.2kb_tiles.counts
### step 0-5-2-3: extract gypsy insertions --- END ---

### step 0-5-2-4: extract TE insertions --- START ---
### remove "Simple_repeat", "Low_complexity", "rRNA" and "tRNA"
awk '{if($11 != "Simple_repeat" && $11 != "Low_complexity" && $11 != "rRNA" && $11 != "tRNA") print}' ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta.out |\
awk '{if($9=="+") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,$9; else if($9=="C") print $5,$6,$7,$5":"$6":"$10":"$11,$2":"$3":"$4,"-"
}' | tr ' ' '\t' > ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta.out.allTE.bed

### annotate tiles that overlap with TE insertions at the sense orientation
bedtools intersect -f 0.5 -s -a ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.200nt.window.bed \
-b ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta.out.allTE.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"TE_S"}' >> ${analysis}/Dpse/dpse_0.2kb_tiles.counts

### annotate tiles that overlap with TE insertions at the antisense orientation
bedtools intersect -f 0.5 -S -a ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.200nt.window.bed \
-b ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta.out.allTE.bed |\
awk '{TE[$4]+=$3-$2} END {for(var in TE) print var,TE[var]/200,"TE_AS"}' >> ${analysis}/Dpse/dpse_0.2kb_tiles.counts
### step 0-5-2-4: extract TE insertions --- END ---

### step 0-5-2-5: mRNA exons --- START ---
### step 0-5-2-5-1: fetch exon sequences from the UCI_Dpse_MV25 assembly
zcat ${references}/Dpse/GCF_009870125.1_UCI_Dpse_MV25_genomic.gtf.gz | tr -d '";' |\
awk 'BEGIN{OFS="\t"} {if($3=="exon") print $1,$4-1,$5,$1":"$4-1":"$5":"$7":"$10":"$12":"$14":"$16,".",$7}'  |\
bedtools getfasta -name -fi ${references}/Dpse/GCF_009870125.1_UCI_Dpse_MV25_genomic.fna -bed - > ${references}/Dpse/GCF_009870125.1_UCI_Dpse_MV25.exons.fasta

### step 0-5-2-5-2: run blastn to predict exons in the Nanopore assembly of the D.pseudoobscura genome
### make the blast database
fastafile="${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta"
makeblastdb -in ${fastafile} -parse_seqids -dbtype nucl -title "Dpse_nanopore_assembly" \
-out ${references}/Dpse/Dpse_nanopore_blast

### input has 384662 exon sequences from the UCI_Dpse_MV25 assembly
input_fasta="${references}/Dpse/GCF_009870125.1_UCI_Dpse_MV25.exons.fasta"
output="${references}/Dpse/GCF_009870125.1_UCI_Dpse_MV25.exons.blast_Dpse_nanopore.csv"
blast_database="${references}/Dpse/Dpse_nanopore_blast"
blastn -db ${blast_database} -outfmt 10 -query ${input_fasta} -out ${output}

### output csv file looks as follows:
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
### output had 189385 sequences
### out of those, 184238 entries had at least one hit with the score of greater than 97 points
### gtf format
awk '{split($1,a,","); if(a[3]>97) print a[1],$0}' ${output} | awk '!seen[$1]++' |\
awk 'BEGIN{OFS="\t"} {split($2,a,","); split(a[1],b,":"); if(a[9]>a[10]) print a[2],"nanopore","exon",a[10],a[9],".","-","0","gene_id \""b[5]"\"; gene_symbol \""b[6]"\"; transcript_id \""b[7]"\"; transcript_symbol \""b[8]"\";" ;
else if(a[10]>a[9]) print a[2],"nanopore","exon",a[9],a[10],".","+","0","gene_id \""b[5]"\"; gene_symbol \""b[6]"\"; transcript_id \""b[7]"\"; transcript_symbol \""b[8]"\";"
}' | sort -k1,1 -k4,4n > ${references}/Dpse/GCF_009870125.1_UCI_Dpse_MV25.exons.blast_Dpse_nanopore.gtf

### gtf file looks as follows:
utg000001l      nanopore        exon    16836   17178   .       -       0       gene_id "LOC117183504"; gene_symbol "XR_004468602.1"; transcript_id "GeneID"; transcript_symbol "117183504";

### convert gtf to bed
cat ${references}/Dpse/GCF_009870125.1_UCI_Dpse_MV25.exons.blast_Dpse_nanopore.gtf | tr -d '";' |\
awk 'BEGIN{OFS="\t"} {print $1,$4,$5,$10":"$12,".",$7}' | sort -k1,1 -k2,2n >  ${references}/Dpse/GCF_009870125.1_UCI_Dpse_MV25.exons.blast_Dpse_nanopore.bed

### step 0-5-2-5-3: add exon annotations
bedtools intersect -wo -s -a ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.200nt.window.bed \
-b ${references}/Dpse/GCF_009870125.1_UCI_Dpse_MV25.exons.blast_Dpse_nanopore.bed |\
awk '{exon[$4]+=$NF} END {for(var in exon) print var,exon[var]/200,"exon"}' >> ${analysis}/Dpse/dpse_0.2kb_tiles.counts
### step 0-5-2-6: mRNA exons --- START ---
### step 0-5: preparing for the 0.2kb genomic tile analyses of the all genome mappers --- END ---


### step 0-6: make a bowtie index for all Drosophila autonomous transposons downloaded from the RepBase in March 2022 (3053 entries)
mkdir -p ${references}/TEs/indices
fasta_formatter -i ${references}/TEs/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.fasta -t |\
awk '{print ">"$1"\n"$NF}' > ${references}/TEs/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.concise.fasta
bowtie-build ${references}/TEs/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.concise.fasta ${references}/TEs/indices/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.concise


### step 0-7: obtain TEs that are non-redundant and produce abundant piRNAs --- START ---
### We used the whole ovary oxidised small RNA library (Dpse_ovary_oxidised_sRNA) to obtain 150 TEs that produce most abundant piRNAs.
lib_sRNA="Dpse_ovary_oxidised_sRNA"
### step 0-7-1: take all reads that mapped to the D.pseudoobscura genome, and use them to map against transposon sequences
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); print ">"$4"\n"a[1]}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.fa

### step 0-7-2: bowtie mapping to TEs, --all --best --strata option, allowing up to three mismatches
TE_index="${references}/TEs/indices/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.concise"
bowtie -f -v 3 --all --best --strata -S ${TE_index} ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_dpse-nanopore_mapped_TE_3MM_mappers.bed

### step 0-7-3: We manually removed 41 TEs that resemble other TEs within the top 150 TEs, which left 109 non-redundant TEs.

### step 0-7-4: check the redundancy between the remaining 109 TEs by counting the piRNA reads that mapped to multiple TE sequences.
mkdir -p ${analysis}/${lib_sRNA}/TE_mapping.109TEs

### step 0-7-4-1: print out all the piRNA reads (gt22nt) mapping to individual TEs
cat ${references}/Dpse/Dpse.non-redundant.top109.TEs | while read TE READS; do
awk -v TE=${TE} '{if($1==TE && $3-$2>22) print $4}' ${analysis}/${lib_sRNA}/${lib_sRNA}_dpse-nanopore_mapped_TE_3MM_mappers.bed |\
sort | uniq > ${analysis}/${lib_sRNA}/TE_mapping.109TEs/${lib_sRNA}_${TE}.mappers.nuc
done

### step 0-7-4-2: collect all TEs mapping piRNA reads that mapped to multiple TEs
cat ${analysis}/${lib_sRNA}/TE_mapping.109TEs/*.mappers.nuc |\
awk '{READS[$1]++} END {for(var in READS) {if(READS[var]>1) print var}}' |\
sort > ${analysis}/${lib_sRNA}/TE_mapping.109TEs/${lib_sRNA}_top109.redundant-mappers.nuc

### step 0-7-4-3: per each TE, measure how many of piRNA mappers mapped to other TEs
cat ${references}/Dpse/Dpse.non-redundant.top109.TEs | while read TE READS; do
printf $TE"\n""total reads count gt22nt: "$READS"\n" >> ${analysis}/${lib_sRNA}/TE_mapping.109TEs/${lib_sRNA}_top109.stats
join -1 1 -2 1 ${analysis}/${lib_sRNA}/TE_mapping.109TEs/${lib_sRNA}_${TE}.mappers.nuc ${analysis}/${lib_sRNA}/TE_mapping.109TEs/${lib_sRNA}_top109.redundant-mappers.nuc |\
awk '{split($1,a,"@"); count+=a[2]} END {print "redundant reads count gt22nt: "count}' >> ${analysis}/${lib_sRNA}/TE_mapping.109TEs/${lib_sRNA}_top109.stats
done
### This analysis confirmed that the non-redundant 109 TEs have more than 10% of redundant piRNA mappers.
### step 0-7: obtain TEs that are non-redundant and produce abundant piRNAs --- END ---


### step 0-8: generate pseudocounts of 1 at every coordinate of non-redundant TEs
mkdir -p ${references}/TEs/TE_pseudo.piRNAs
### run it per each TE
cat ${references}/Dpse/Dpse.non-redundant.top109.TEs | while read TE READS; do
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
index_genome="${references}/Dpse/indices/dpse-nanopore"
bowtie -f -v 1 -m 1 -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - | sort -k1,1 -k2,2n > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dpse-nanopore.bed


### step 4: remove tRNA reads by intersection
### exclude reads that intersect with tRNA annotations from the genome unique mappers
bedtools intersect -v -s -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dpse-nanopore.bed \
-b ${references}/Dpse/dpse-nanopore.tRNAs-extended.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dpse-nanopore.tRNA-excluded.bed


### step 5: generate the bigwig files for the genome unique mappers --- START ---
### step 5-1: make bedgraph files
### only count reads that are greater than 22nt
### for IP libraries and oxidised libraries --- START ---
### normalise read counts to one million genome unique mappers
TOTAL=`(awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]} END {print count
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dpse-nanopore.tRNA-excluded.bed)`
### bedgraph output per 1Mio miRNA mappers (plus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dpse-nanopore.tRNA-excluded.bed |\
bedtools genomecov -strand + -split -bg -i - -g ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.sizes |\
awk -v TOTAL=${TOTAL} '{print $1,$2,$3,$4/TOTAL*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_plus.dpse-nanopore.tRNA-excluded.bg
### bedgraph output per 1Mio miRNA mappers (minus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dpse-nanopore.tRNA-excluded.bed |\
bedtools genomecov -strand - -split -bg -i - -g ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.sizes |\
awk -v TOTAL=${TOTAL} '{print $1,$2,$3,-$4/TOTAL*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_minus.dpse-nanopore.tRNA-excluded.bg
### for IP libraries and oxidised libraries --- END ---

### for non-oxidised libraries --- START ---
### normalise read counts to one million microRNA reads
MIRNA=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-mapped.bed | grep mir | tr '@' ' ' | awk '{count+=$5} END {print count}')`
### bedgraph output per 1Mio miRNA mappers (plus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dpse-nanopore.tRNA-excluded.bed |\
bedtools genomecov -strand + -split -bg -i - -g ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.sizes |\
awk -v MIRNA=${MIRNA} '{print $1,$2,$3,$4/MIRNA*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_plus.dpse-nanopore.tRNA-excluded.bg
### bedgraph output per 1Mio miRNA mappers (minus strand)
awk '{split($4,a,"@"); if($3-$2 > 22) {for(i=1;i<=a[2];i++) print}}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dpse-nanopore.tRNA-excluded.bed |\
bedtools genomecov -strand - -split -bg -i - -g ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.sizes |\
awk -v MIRNA=${MIRNA} '{print $1,$2,$3,-$4/MIRNA*1000000}' | tr ' ' '\t' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_minus.dpse-nanopore.tRNA-excluded.bg
### for non-oxidised libraries --- END ---

### step 5-2: make bigwig files
for strand in plus minus; do
bedGraphToBigWig ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_${strand}.dpse-nanopore.tRNA-excluded.bg \
${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.sizes ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers_${strand}.dpse-nanopore.tRNA-excluded.bw
done
### step 5: generate the bigwig files for the genome unique mappers --- END ---


### step 6: 0.5kb tile analysis for the genome unique mappers --- START ---
### count piRNA reads (greater than 22nt) that uniquely mapped to the 0.5kb tiles
READS=`(awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]} END {print count
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dpse-nanopore.tRNA-excluded.bed)`
bedtools intersect -wo -s -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dpse-nanopore.tRNA-excluded.bed \
-b ${references}/Dpse/Dpse_nanopore_85percent_0.5kbtiles.25mers.bed |\
awk -v READS=${READS} -v LIB=${lib_sRNA} '{split($4,a,"@"); if(length(a[1])>22) TILE[$10]+=a[2]} END {for(var in TILE) print var,TILE[var]/READS*1000000,LIB
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dpse-nanopore.tRNA-excluded.gt22nt.0.5kbtiles.counts


### step 7: 0.2kb tile analysis for the genome all mappers --- START ---
### step 7-1: map reads to the Nanopore assembly of the D. pseudoobscura genome, allowing up to 1MM, all best strata --- START ---
index_genome="${references}/Dpse/indices/dpse-nanopore"
bowtie -f -v 1 --all --best --strata -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.bed

### step 7-2: add mapping counts to the fasta entries. Mapping counts will be used to distribute the read coverage across tiles.
awk '{READS[$4]++} END {for(var in READS) {split(var,a,"@"); print ">"var":"READS[var]"\n"a[1]}
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.mapping-counts.fa

### step 7-3: map again all best strata
index_genome="${references}/Dpse/indices/dpse-nanopore"
bowtie -f -v 1 --all --best --strata -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.mapping-counts.fa |\
samtools view -bS - | bamToBed -i - | sort -k1,1 -k2,2n > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.bed

### step 7-4: exclude reads that intersect with tRNA annotations
bedtools intersect -v -s -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.bed \
-b ${references}/Dpse/dpse-nanopore.tRNAs-extended.bed > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.bed

### step 7-5: collect read counts per 200nt windows from all-best-strata
TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if($3-$2>22) count+=b[1]} END {print count}')`
awk 'BEGIN{OFS="\t"} {if($3-$2 > 22) print}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.bed |\
bedtools intersect -s -wo -F 0.5 -a ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.200nt.window.bed -b - |\
awk -v TOTAL=${TOTAL} -v LIB=${lib_sRNA} '{split($10,a,"@"); split(a[2],b,":"); TILE[$4]+=b[1]/b[2]} END {for(var in TILE) print var,TILE[var]/TOTAL*1000000,LIB"_abs"
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.gt22nt.0.2kbtiles.counts
### step 7: 0.2kb tile analysis for the genome all mappers --- END ---


### step 8: analyse TE mapping reads --- START ---
### step 8-1: take all reads that mapped to the D.pseudoobscura Nanopore assembly, and use them to map against transposon sequences
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); print ">"$4"\n"a[1]}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.fa

### step 8-2: bowtie mapping to TEs, --all --best --strata option, allowing up to three mismatches
TE_index="${references}/TEs/indices/Drosophila_all-RepBase-autonomous-Transposon-entries_Mar2022.concise"
bowtie -f -v 3 --all --best --strata -S ${TE_index} ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_dpse-nanopore_mapped_TE_3MM_mappers.bed

### step 8-3: making TE mapping plots and measuring pingpong and phasing signatures --- START ---
mkdir -p ${analysis}/${lib_sRNA}/End.plots.TEs.linkage
### step 8-3-1: normalise the counts of 3' end 5' ends per 1 million total genome mappers (>22nt)
TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if(length(a[1])>22) count+=b[1]} END {print count}')`

### step 8-3-2: count 5end and 3end of plus and minus mapping reads (>22nt)
### run it per each TE
cat ${references}/Dpse/Dpse.non-redundant.top109.TEs | while read TE READS; do
### use pseudocounts created in step 0-8
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_dpse-nanopore_mapped_TE_3MM_mappers.bed ${references}/TEs/TE_pseudo.piRNAs/${TE}_pseudo.piRNAs.txt |\
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
### step 8: analyse TE mapping piRNAs --- END ---


### step 9: measure in-trans ping-pong linkage --- START ---
mkdir -p ${analysis}/${lib_sRNA}/jellyfish

### step 9-1: extract g1g9, g2g10_revComp and last9 sequences from piRNA genome unique mappers
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr($4,1,9)
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dpse-nanopore.tRNA-excluded.bed > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_unique_g1g9.fasta
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr($4,2,9)
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dpse-nanopore.tRNA-excluded.bed | fastx_reverse_complement - > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_unique_g2g10_revComp.fasta
awk '{split($4,a,"@"); if(length(a[1])>22) for(i=1;i<=a[2];i++) print ">"$4":"i"\n"substr(a[1],length(a[1])-8,9)
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dpse-nanopore.tRNA-excluded.bed > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_unique_last9.fasta

### step 9-2: extract g1g9, g2g10_revComp and last9 sequences from piRNA genome all mappers
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if($3-$2>22) for(i=1;i<=b[1];i++) print ">"$4":"i"\n"substr($4,1,9)}'  > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_all_g1g9.fasta
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if($3-$2>22) for(i=1;i<=b[1];i++) print ">"$4":"i"\n"substr($4,2,9)}' | fastx_reverse_complement - > ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_all_g2g10_revComp.fasta
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.bed |\
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
Rscript /scratch/lf10/rh1772/fly_piRNAs/scripts/noYb-paper/jellyfish_linkage_Dspp_mouse.R DIRECTORY="${analysis}" LIB=${lib_sRNA} CLASS=${CLASS}
done
### step 9: measure in-trans ping-pong linkage --- END ---

### common analyses for all small RNA libraries --- END ---



### post-processing of the analysed data for small RNA libraries --- START ---

### step 10: 0.5kb tile analysis of the genome unique mappers
### step 10-1: add annotations of piRNA cluster tiles to the tiles
mkdir -p ${analysis}/Dpse
bedtools intersect -wo -s -a ${references}/Dpse/Dpse_nanopore_85percent_0.5kbtiles.25mers.bed \
-b ${references}/Dpse/Dpse_nanopore_piRNA-clusters.bed |\
awk '{print $4,$NF/500,$10}' > ${analysis}/Dpse/dpse_0.5kb_tiles.counts

### step 10-2: add counts from small RNA libraries
for lib_sRNA in <oxidied whole ovary small RNA library> <oxidied embryonic small RNA library>; do
#for lib_sRNA in Dpse_ovary_oxidised_sRNA_abs Dpse_embryo_oxidised_sRNA_abs; do
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.dpse-nanopore.tRNA-excluded.gt22nt.0.5kbtiles.counts >> ${analysis}/Dpse/dpse_0.5kb_tiles.counts
done
### step 10-3: tile_analysis_plots.R was used to make a scatter plot in Figure 5.


### step 11: 0.2kb tile analysis of the genome all mappers --- START ---
### step 11-1: add counts from small RNA libraries
for lib_sRNA in <oxidied whole ovary small RNA library> <oxidied embryonic small RNA library>; do
#for lib_sRNA in Dpse_ovary_oxidised_sRNA_abs Dpse_embryo_oxidised_sRNA_abs; do
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.gt22nt.0.2kbtiles.counts >> ${analysis}/Dpse/dpse_0.2kb_tiles.counts
done

### step 11-2: make a table using reshape2 in R and count piRNA reads from individual annotations
#R/4.0.0
library(reshape2)
table <- read.table("${analysis}/Dpse/dpse_0.2kb_tiles.counts", header=F)
table_d=dcast(table,V1~V3,value.var="V2")
table_d[is.na(table_d)] <- 0
write.table(table_d, file="${analysis}/Dpse/dpse_0.2kb_tiles.counts.table", quote=F, col.names=T, row.names = F)
colnames(table_d)

[1] "V1"                            "cluster3l_uni"
[3] "cluster40l_dual"               "Dpse_embryo_oxidised_sRNA_abs"
[5] "Dpse_ovary_oxidised_sRNA_abs" "exon"
[7] "gypsy_AS" "gypsy_S"
[9] "TE_AS" "TE_S"
#R/4.0.0

### count total somatic reads
soma_total=`(awk '{if(NR>1 && $5>0) print}' ${analysis}/Dpse/dpse_0.2kb_tiles.counts.table |\
awk '{if($4/$5 < 0.1) print}' | awk '{count+=$5} END {print count}')`

### count somatic uni-stranded cluster reads
uni_total=`(awk '{if(NR>1 && $5>0) print}' ${analysis}/Dpse/dpse_0.2kb_tiles.counts.table |\
awk '{if($4/$5 < 0.1) print}' | awk '{if($2>0) count+=$5} END {print count}')`

### count exonic piRNA reads
### exclude tiles that overlap with uni clusters and TE annotations
exons=`(awk '{if(NR>1 && $5>0) print}' ${analysis}/Dpse/dpse_0.2kb_tiles.counts.table |\
awk '{if($4/$5 < 0.1) print}' | awk '{if($9==0 && $10==0 && $2==0 && $6>0) count+=$5} END {print count}')`

### count somatic gypsy_AS reads outside the clusters
no_uni_gypsy_AS=`(awk '{if(NR>1 && $5>0) print}' ${analysis}/Dpse/dpse_0.2kb_tiles.counts.table |\
awk '{if($4/$5 < 0.1) print}' | awk '{if($2==0 && $7>0 && $8==0) count+=$5} END {print count}')`

### count somatic gypsy_S reads outside the clusters
no_uni_gypsy_S=`(awk '{if(NR>1 && $5>0) print}' ${analysis}/Dpse/dpse_0.2kb_tiles.counts.table |\
awk '{if($4/$5 < 0.1) print}' | awk '{if($2==0 && $8>0 && $7==0) count+=$5} END {print count}')`

### collect the stats
others=`echo "$soma_total - $uni_total - $exons - $no_uni_gypsy_AS - $no_uni_gypsy_S" | bc`
printf $soma_total" Dpse soma_total""\n"$uni_total" Dpse uni_total""\n"$exons" Dpse exons""\n"$no_uni_gypsy_AS" Dpse no_uni_gypsy_AS""\n"$no_uni_gypsy_S" Dpse no_uni_gypsy_S""\n"$others" Dpse others""\n" >> ${analysis}/0.2kb_tiles.counts.stats
### step 11: 0.2kb tile analysis of the genome all mappers --- END ---


### step 12: comparing the abundance of somatic and germline TE mapping piRNA reads (>22nt) --- START ---
for lib_sRNA in <oxidied whole ovary small RNA library> <oxidied embryonic small RNA library>; do
#for lib_sRNA in Dpse_ovary_oxidised_sRNA_abs Dpse_embryo_oxidised_sRNA_abs; do
### normalise the counts to 1 million total genome mappers (>22nt)
TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); if(length(a[1])>22) count+=b[1]} END {print count}')`

### step 12-1: count piRNA mappers (>22nt) both sense and antisense
cat ${references}/Dpse/Dpse.non-redundant.top109.TEs | while read TE READS; do
## some TE sequences have repeated sequences within and the same reads can map multiple times
awk '!seen[$1" "$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_dpse-nanopore_mapped_TE_3MM_mappers.bed |\
awk -v TE=${TE} -v TOTAL=${TOTAL} -v LIB=${lib_sRNA} '{split($4,a,"@"); split(a[2],b,":"); if($3-$2>22 && $1==TE) count+=b[1]} END {print TE,count/TOTAL*1000000,LIB
}' >> ${analysis}/${lib_sRNA}/${lib_sRNA}_109TE_mappers.counts
done

### step 12-2: combine counts from the whole ovary and embryonic small RNA libraries
cat ${analysis}/${lib_sRNA}/${lib_sRNA}_109TE_mappers.counts >> ${analysis}/Dpse_109TE_mappers.counts
done

### step 12-3: TE-piRNAs_analysis_plots.R was used to make a scatter plot in Extended Data Figure 6.

### step 12-4: identify TEs that are enriched more than three times in the ovary sRNAseq library than in the embryonic sRNAseq library
#R/4.0.0
library(reshape2)
table <- read.table("${analysis}/Dpse_109TE_mappers.counts", header=F)
table.d <- dcast(table,V1~V3,value.var="V2")
write.table(table.d, file="${analysis}/Dpse_109TE_mappers.counts.table", quote=F, row.names=F)

awk '{if(NR>1 && $3>$2*3) print $1}' ${analysis}/Dpse_109TE_mappers.counts.table
Gypsy-11_DBp-I
Gypsy-12_DAzt-I
Gypsy-26_DPse-I
Gypsy-26_DPse-LTR
Gypsy12-I_Dpse
### step 12: comparing the abundance of somatic and germline TE mapping piRNA reads (>22nt) --- END ---


### step 13: measuring the strand bias of the piRNA reads (>22nt) mapping to non-redundant 94 TEs
### we used the oxidised whole ovary small RNA library
lib_sRNA="<oxidied whole ovary small RNA library>"
lib_sRNA="Dpse_ovary_oxidised_sRNA"
### exclude somatic TEs
awk '{if($1 != "Gypsy-11_DBp-I" && $1 != "Gypsy-12_DAzt-I" && $1 != "Gypsy-26_DPse-I" && $1 != "Gypsy-26_DPse-LTR" && $1 != "Gypsy12-I_Dpse") print
}' ${references}/Dpse/Dpse.non-redundant.top109.TEs | while read TE READS; do
## some TE sequences have repeated sequences within and the same reads can map multiple times
awk '!seen[$1" "$4" "$6]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_dpse-nanopore_mapped_TE_3MM_mappers.bed |\
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


### step 15: quantify miRNA reads, TE AS and S mapping piRNAs (gt22nt) from non-oxidised sRNAseq libraries
lib_sRNA="<unoxidied whole ovary small RNA library>"
lib_sRNA="Dpse_ovary_unoxidised_sRNA"
### miRNA reads
miRNA=`(cat ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-mapped.bed | grep mir | awk '!seen[$4]++' |\
awk '{split($4,a,"@"); count+=a[2]} END {print count}')`

TOTAL=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.bed |\
awk '{split($4,a,"@"); split(a[2],b,":"); {if($3-$2>22) count+=b[1]}} END {print count}')`

### quantify sense and antisense TE piRNAs
TE_S_all=`(bedtools intersect -s -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.bed \
-b ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta.out.allTE.bed | awk '!seen[$4]++' |\
awk '{split($4,a,"@"); split(a[2],b,":"); {if($3-$2>22) count+=b[1]}} END {print count}')`

TE_AS_all=`(bedtools intersect -S -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.bed \
-b ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta.out.allTE.bed | awk '!seen[$4]++' |\
awk '{split($4,a,"@"); split(a[2],b,":"); {if($3-$2>22) count+=b[1]}} END {print count}')`

TE_all=`(bedtools intersect -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped.dpse-nanopore.all-best-strata.mapping-counts.tRNA-excluded.bed \
-b ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta.out.allTE.bed | awk '!seen[$4]++' |\
awk '{split($4,a,"@"); split(a[2],b,":"); {if($3-$2>22) count+=b[1]}} END {print count}')`

TE_S=`echo "$TE_all - $TE_AS_all" | bc`
TE_AS=`echo "$TE_all - $TE_S_all" | bc`
TE_S_AS=`echo "$TE_S_all + $TE_AS_all - $TE_all" | bc`
others=`echo "$TOTAL - $TE_all" | bc`

awk -v miRNA=${miRNA} -v TE_S=${TE_S} -v TE_AS=${TE_AS} -v TE_S_AS=${TE_S_AS} -v others=${others} -v LIB=${lib_sRNA} '{if(NR==1) print "TE_S",TE_S/miRNA,LIB"\n""TE_AS",TE_AS/miRNA,LIB"\n""TE_S_AS",TE_S_AS/miRNA,LIB"\n""others",others/miRNA,LIB
}' ${analysis}/TE_S_AS_bias.stats >> ${analysis}/miRNA_TE-mappers.stats
### TE-piRNAs_analysis_plots.R was used to make a box plot in Figure 6.


### step 16: collecting the in-trans pingpong linkage values across libraries
for CLASS in unique all; do
for lib_sRNA in w1118_ovary_oxidised_sRNA Dpse_ovary_oxidised_sRNA Deug_ovary_oxidised_sRNA_rep1 Deug_ovary_oxidised_sRNA_rep2 SRR1746887 SRR1104823; do
awk -v LIB=${lib_sRNA} -v CLASS=${CLASS} '{if(NR>1) print $1,$2,LIB"@"CLASS
}' ${analysis}/${lib_sRNA}/jellyfish/${lib_sRNA}_${CLASS}_all_m9.linkage.txt >> ${analysis}/in-trans-pingpong.stats
done
done
### in-trans-pingpong_plots.R was used to make a barchart in Figure 6 and Extended Data Figure 9.

### post-processing of the analysed data for small RNA libraries --- END ---

### processing and analysing the small RNA sequencing libraries --- END ---



### processing and analysing the paired end RNA seq library of poly-A+ RNA --- START ---

### below is the RNA library processed in this part of the script:
Dpse_ovary_polyA

### RNA seq libraries
5’ AATGATACGGCGACCACCGAGATCTACAC[TCTTTCCCTACACGACGCTCTTCCGATCT]--- [R1 ->] --- [<- R2] ---[AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC] ---NNNNNN--- ATCTCGTATGCCGTCTTCTGCTTG 3’
R1: reverse complement
R2: sense

### step 1: trim adapters, filter poor-quality reads, and take only the reads that have mates
lib_pA="Dpse_ovary_polyA"
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

### step 2: make a star index using the gtf file generated using the blastn results described in the step 0-5 of the small RNA sequencing analysis
mkdir -p ${references}/Dpse/indices/Dpse_nanopore_exon_STAR
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ${references}/Dpse/indices/Dpse_nanopore_exon_STAR \
--genomeSAindexNbases 12 \
--genomeFastaFiles ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.fasta \
--sjdbGTFfile ${references}/Dpse/GCF_009870125.1_UCI_Dpse_MV25.exons.blast_Dpse_nanopore.gtf \
--sjdbOverhang 100

### step 3: STAR mapping
mkdir -p ${analysis}/${lib_pA}/STAR
STAR --runThreadN 8 --genomeDir ${references}/Dpse/indices/Dpse_nanopore_exon_STAR \
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
-g ${references}/Dpse/indices/Dpse_nanopore_exon_STAR/chrNameLength.txt |\
awk -v READS=${READS} '{print $1,$2,$3,$4/READS*1000000}' | tr ' ' '\t' > ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_plus.bg
### bedgraph output per 1Mio unique mappers (minus strand)
bedtools genomecov -strand - -split -bg -i ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR.sorted.bed \
-g ${references}/Dpse/indices/Dpse_nanopore_exon_STAR/chrNameLength.txt |\
awk -v READS=${READS} '{print $1,$2,$3,-$4/READS*1000000}' | tr ' ' '\t' > ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_minus.bg

### generate bigWig files
for strand in plus minus; do
bedGraphToBigWig ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_${strand}.bg \
${references}/Dpse/indices/Dpse_nanopore_exon_STAR/chrNameLength.txt ${analysis}/${lib_pA}/STAR/${lib_pA}_STAR_${strand}.bw
done
### step 4: generate bigwig files for the IGV browser --- END ---
### processing and analysing the paired end RNA seq library of poly-A+ RNA --- END ---


### processing and analysing ChIPseq libraries --- START ---

### below is the ChIPseq libraries processed in this part of the script:
Dpse_ChIP_input_rep1
Dpse_ChIP_input_rep2
Dpse_ChIP_polII-CTD_rep1
Dpse_ChIP_polII-CTD_rep2

### ChIP seq libraries
5’ AATGATACGGCGACCACCGAGATCTACAC[TCTTTCCCTACACGACGCTCTTCCGATCT]--- [R1 ->] --- [<- R2] ---[AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC] ---NNNNNN--- ATCTCGTATGCCGTCTTCTGCTTG 3’

### step 1: trim adapters and filter poor-quality reads
mkdir -p ${analysis}/${lib_ChIP}
fastx_clipper -l 35 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-i <(zcat ${raw_data}/${lib_ChIP}_R1.fastq.gz | fastq_quality_filter -q 33 -p 40) | gzip > ${analysis}/${lib_ChIP}/${lib_ChIP}_R1.trimmed.fastq.gz

### step 2: map R1 reads to the D.pseudoobscura Nanopore assembly, allowing up to 3 mismatches, unique mappers
index_genome="${references}/Dpse/indices/dpse-nanopore"
bowtie -q -v 3 -m 1 -S ${index_genome} ${analysis}/${lib_ChIP}/${lib_ChIP}_R1.trimmed.fastq.gz |\
samtools view -bS - | bamToBed -i - | sort -k1,1 -k2,2n > ${analysis}/${lib_ChIP}/${lib_ChIP}_genome-unique-mappers.bed

### step 3: generate bigwig files for the IGV browser --- START ---
READS=`(cat ${analysis}/${lib_ChIP}/${lib_ChIP}_genome-unique-mappers.bed | wc -l | awk '{print $NF}')`
### bedgraph output per 1Mio unique mappers
bedtools genomecov -bg -i ${analysis}/${lib_ChIP}/${lib_ChIP}_genome-unique-mappers.bed \
-g ${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.sizes |\
awk -v READS=${READS} '{print $1,$2,$3,$4/READS*1000000}' | tr ' ' '\t' > ${analysis}/${lib_ChIP}/${lib_ChIP}_genome-unique-mappers.bg

### generate bigwig
bedGraphToBigWig ${analysis}/${lib_ChIP}/${lib_ChIP}_genome-unique-mappers.bg \
${references}/Dpse/Dpse.pass.minimap2.racon.x3.pilon.x3.trimmed.sizes ${analysis}/${lib_ChIP}/${lib_ChIP}_genome-unique-mappers.bw
### step 3: generate bigwig files for the IGV browser --- END ---
### processing and analysing ChIPseq libraries --- END ---
