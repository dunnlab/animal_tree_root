#!/bin/sh
name=/proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/fasta/component1
thread=1
module load bioinfo-tools
module load MAFFT
module load iqtree/1.5.3-omp

mafft-linsi $name.fasta > $name.mafft
iqtree -s $name.mafft -m LG+C60+F -nt $thread -bb 1000 -redo

