#!/bin/sh
thread=1
module load bioinfo-tools
module load MAFFT
module load iqtree/2.0-rc2-omp-mpi

mkdir -p /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft


while read line;
do

        mkdir -p /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test_paralogy/$line/
        mafft-linsi /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/fasta/paralogy/$line.fasta > /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test_paralogy/$line/$line.mafft
	mkdir -p /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test_paralogy/iq_tree/$line/
        iqtree -s /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test_paralogy/$line/$line.mafft -m TEST -nt $thread -bb 1000 -pre  /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test_paralogy/iq_tree/$line/TEST_$line -t PARS ;
done < /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/name_list_paralogs_new.txt




