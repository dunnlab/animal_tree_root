#!/bin/sh
thread=1
module load bioinfo-tools
module load MAFFT
module load iqtree/1.5.3-omp

mkdir -p /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft


while read line;
do

        mkdir -p /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test/$line/
        mafft-linsi /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/fasta/$line.fasta > /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test/$line/$line.mafft
        iqtree -s /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test/$line/$line.mafft -m TEST -nt $thread -bb 1000 -pre  /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test/$line/TEST_$line ;
done < /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/name_list_test.txt




