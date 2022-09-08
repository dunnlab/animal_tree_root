#!/bin/sh
thread=1
module load bioinfo-tools
module load MAFFT
module load iqtree/1.5.3-omp

mkdir -p /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft

mkdir -p /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test/component1081/
mafft-linsi /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/fasta/component1081.fasta > /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test/component1081/component1081.mafft
iqtree -s /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test/component1081/component1081.mafft -m TEST -nt $thread -bb 1000 -pre  /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test/component1081/TEST_component1081 ;




