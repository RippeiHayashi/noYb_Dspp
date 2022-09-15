### this code is to find the Yb homologs in Drosophila genome assemblies

### ncbi-blast-2.9.0 is downloaded from NCBI through the link below.
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-src.tar.gz

DIRECTORY="<the directory where you save the genome sequences>"
ASSEMBLY="<genome assembly name>"
SPECIES="<species name>"

### unzip the genomic sequence
zcat ${DIRECTORY}/${ASSEMBLY}_genomic.fna.gz > ${DIRECTORY}/${ASSEMBLY}_genomic.fna

### make the blast database
mkdir -p ${DIRECTORY}/blastdb
genome_fasta=${DIRECTORY}"/"${ASSEMBLY}"_genomic.fna"
makeblastdb -in ${genome_fasta} -parse_seqids -dbtype nucl -title "${SPECIES}_NCBI_assembly" \
-out ${DIRECTORY}/blastdb/${ASSEMBLY}_NCBI_blast

### run blastn
mkdir -p ${DIRECTORY}/tblastn_results
protein_fasta="<path to the query protein sequence file>"
protein_name="<name of the query protein>"
tblastn -db ${DIRECTORY}/blastdb/${ASSEMBLY}_NCBI_blast \
-query ${protein_fasta} \
-out ${DIRECTORY}/tblastn_results/${ASSEMBLY}_${SPECIES}_${protein_name}_tblastn_results.out
