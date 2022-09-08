#!/bin/bash

module load bioinfo-tools
module load diamond
module load augustus/3.3.3-CGP
module load BLAST+
module load blast_databases #Doing module load blast_databases sets the environment variable BLASTDB to this directory

diamond blastp -d db_test/animal_root_test_diamond.dmnd -q fasta_test/dosReisDonoghueYang_2015_superalignment.fa -o diamond_results/dosReisDonoghueYang_2015_superalignment_permissive.txt -f 6  --gapopen 6 --gapextend 2

#blastp -db swissprot -query fasta_test/dosReisDonoghueYang_2015_superalignment.fa -out swissprot_results/dosReisDonoghueYang_2015_superalignment_swissprot.txt -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'
