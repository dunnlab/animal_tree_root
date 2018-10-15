#!/bin/bash
#SBATCH -J iqtree
#SBATCH -p general   
#SBATCH -C avx2
#SBATCH -N 1
#SBATCH -c 20
#SBATCH -t 7-00:00:00
#SBATCH --mem 100G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuanning.li@yale.edu
#SBATCH --output=iqtree.%j.out
#SBATCH --error=iqtree.%j.err

for FILENAME in *.nex
do
iqtree -s $FILENAME -nt 20 -bb 1000 -madd GTR+C60+F+R9,LG+C60+F+R9,LG+FO*H4,GTR+FO*H4, -pre $FILENAME.iqtree 
done