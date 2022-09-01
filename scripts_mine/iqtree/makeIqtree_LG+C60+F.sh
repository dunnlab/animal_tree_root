#!/bin/sh

thread=1
module load bioinfo-tools
module load MAFFT
module load iqtree/1.5.3-omp

mkdir -p /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/LG+C60+F_results/
mkdir -p /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/LG+C60+F_results/${1}/
mkdir -p /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/LG+C60+F_results/iqtree/${1}
#mafft-linsi /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/fasta/all_datasets/all_components/${1}.fasta > /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/LG+C60+F_results/mafft/${1}.mafft
iqtree -s /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/LG+C60+F_results/mafft/${1}.mafft -m LG+C60+F -nt $thread -bb 1000 -pre  /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/LG+C60+F_results/iqtree/${1}/LG+C60+F_${1}
