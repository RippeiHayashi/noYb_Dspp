### the code is to run RepeatMasker and find insertions of env carrying gypsy elements

### fastx_toolkit
http://hannonlab.cshl.edu/fastx_toolkit/
### repeatmasker 4.1.0
http://www.repeatmasker.org/RMDownload.html
### download tandem repeat finder
http://tandem.bu.edu/trf/trf409.linux64.download.html
### ncbi-blast-2.9.0 is downloaded from NCBI through the link below.
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-src.tar.gz

DIRECTORY="<the directory where you save the genome sequences>"
ASSEMBLY="<genome assembly name>"
SPECIES="<species name>"

### unzip the genomic sequence
### fasta_formatter is from fastx_toolkit
zcat $DIRECTORY/${ASSEMBLY}_genomic.fna.gz |\
fasta_formatter - -t | awk '{print ">"$1"\n"$NF}' > $DIRECTORY/${ASSEMBLY}_genomic.trimmed.fna

### Run RepeatMasker 4.1.0 by specifying "drosophila"
fastafile="$DIRECTORY/${ASSEMBLY}_genomic.trimmed.fna"
RepeatMasker -species "drosophila" -parallel 10 -xsmall ${fastafile}


### Run tBlastn to search for gypsy ORFs, see also annotating_Yb.sh
### -outfmt 10 is to obtain outputs in a tabulated format
mkdir -p ${DIRECTORY}/tblastn_results
protein_fasta="<path to the query protein sequence file>"
protein_name="<name of the query protein>"
tblastn -db ${DIRECTORY}/blastdb/${ASSEMBLY}_NCBI_blast -outfmt 10 \
-query ${protein_fasta} \
-out ${DIRECTORY}/tblastn_results/${ASSEMBLY}_${SPECIES}_${protein_name}_tblastn_results.tab.out

### gypsy ORFs that were retrieved from the RepBase and examined in this study
Gypsy10-I_Dpse_1p
Gypsy10-I_Dpse_2p
Gypsy10-I_Dpse_3p
Gypsy17-I_Dpse_2p
Gypsy17-I_Dpse_3p
Gypsy17-I_Dpse_4p
Gypsy12-I_Dpse_1p
Gypsy12-I_Dpse_2p
Gypsy12-I_Dpse_4p
Gypsy-101_DAzt-I_1p
Gypsy-101_DAzt-I_2p
Gypsy-101_DAzt-I_3p
Gypsy-19_DAzt-I_1p
Gypsy-19_DAzt-I_2p
Gypsy-19_DAzt-I_3p
Gypsy-3_DAzt-I_1p
Gypsy-3_DAzt-I_2p
Gypsy-3_DAzt-I_3p
Gypsy-3_DAzt-I_4p
Gypsy-8_DAzt-I_1p
Gypsy-8_DAzt-I_2p
Gypsy-8_DAzt-I_3p
Gypsy_DS_1p_gag-pol
Gypsy1-I_DM_1p
Gypsy1-I_DM_2p
Gypsy1-I_DM_3p
Gypsy-6_DEu-I_1p
Gypsy-6_DEu-I_2p
Gypsy-6_DEu-I_3p
Gypsy-37_DEl-I_1p
Gypsy-37_DEl-I_2p
Gypsy-37_DEl-I_3p

### gypsy env ORFs that were newly annotated and examined in this study, see also Supplementary Table 2
>Gypsy12_Dpse_env_30541304
MIILLVALVNARITDYSHSDYVPILDGDILVWDEINYLRHSTNLTDYERMADETANLTEMFPQSHMRKLLVVDTDHIRNMLATISVHHRVARSLNILGSVLKVVAGTPDADDLEKIRINEAQLIESNNRQISINSKSQEQINRLTDSVNKLLEAAKGKQIDSAHLYETLLARNRMLASELSNLMLTIS
LAKVNVINPVILDHDDLNSIFSNQLTNVIVTNILEVSKIKVFQSNSIIHFVIQFPKIKYICKKITIFPVAHNGTVLRLDDNIVADCNDQIVTVSDCKQTTTTTFCETSTKDSCAQGLYSGGVAHCQSQPSHLSAITLVDDGIIIINDHPAAVSFDGSAALNISGTHLITFNDYAVINGSRYQNRKNVQ
SRYPGVASSPLLNVTEHKRVLSLPFLHQLSEENLNFIKEIKEEVSSRSRPIFAFCLGLGICGLVCGMAMLRLYLTKKRDARQINGLMARLSAPGTATAQGGE
>Gypsy_DS_Dbif-env
LLVVLAAVNARITDFSHANYIPVLDGEVLVFDQRSYLRHSSNISEFISMIDETEKLSDSFPQSHMRKLLDVDTDHLRTLLSVLQVHHRFARSLDFLGTALKVVAGTPDASDFLKVRVTEAQLVESNSKQIIINSETQKQINRLTDTINKIISSRKGDLVDTPHLFETLLARNRILNTEIQNLILTITL
AKANIVNPTILDHADLKSLIEQDTPIVSLLEASKIKVLQSENIIHILIAYPKVEFKCQKVSVYPVSHQQTILRLDEDTLAECERDTFAVTGCTVTTHNTFCERARRETCASSLHAGNTANCHTQPSHLNAIMPIDDGVVVINEATARVRTDDGAEVTVSGTFLITFERSAAINGTEFINLRKAPSKQP
GTVRSPLLNIIGHDPALSIPLLHRMNINNLQSILDFKEEVIAAGSPKFWFAVGAVLNVGLICSFILFMALRRKRASLRIQKALDNFNMTEDGHHSEGG


### genome assemblies examined in this study
<accession> <assembly name> <species name>
GCF_009870125.1 UCI_Dpse_MV25 D.pseudoobscura
GCA_009664405.1 UCBerk_Dbif_1.0 D.bifasciata
GCA_005876895.1 DaztRS1 D.azteca
GCF_018153835.1 ASM1815383v1 D.eugracilis
dm6 r6.31 D.melanogaster
GCF_018153725.1 ASM1815372v1 D.mojavensis

### Nanopore assemblies of D.pseudoobscura and D.eugracilis from Miller et al., G3, 2018
https://academic.oup.com/g3journal/article/8/10/3131/6027000
https://github.com/danrdanny/Drosophila15GenomesProject/tree/master/assembledGenomes
