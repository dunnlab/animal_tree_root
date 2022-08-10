#! /bin/bash
name=Tribe1020_Tribe1468_0
thread=16
#SBATCH -J $name
#SBATCH --mail-user=Katarzyna.Zaremba@icm.uu.se
#SBATCH --mail-type=ALL
#SBATCH -A "snic2022-22-694" 
#SBATCH -p node -n 16
#SBATCH -t 5:00:00
#SBATCH -e $name.e
#SBATCH -o $name.o
module load bioinfo-tools
module load MAFFT/7.407
module load iqtree/1.5.3-omp
mafft-linsi --thread $thread $name.fa > $name.mafft
#not available?
#trimal -in $name.mafft -out $name.gappyout.fasta -fasta -gappyout
#iqtree-omp -s $name.mafft -m LG+C60+F -nt $thread -bb 1000  

