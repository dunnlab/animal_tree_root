#!/bin/sh
thread=1
module load bioinfo-tools
module load MAFFT
module load iqtree/1.5.3-omp

#mkdir mafft/test

while read line;
do

        mkdir -p /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test/$line/
        mafft-linsi /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/fasta/$line.fasta > /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test/$line/$line.mafft
        iqtree -s /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test/$line/$line.mafft -m LG+C60+F -nt $thread -bb 1000 -pre  /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test/$line/LG+C60+F_$line ;
done < /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/name_list_test.txt

