# Animal tree root

A comparison of phylogenetic studies relevant to placing the root of the animal phylogeny.

## Repo Overview

(See each directory for more documentation)

``` text
.
├── data_processed         # Matrices and data tables used for analyses in this manuscript
│   ├── matrices           # Matrices in consistent formats with harmonized taxon names
│   └── tables             # Tabular summaries of previously published datasets and results, reconciliation results (taxon-clade maps, partition-gene maps, etc...)
├── docker                 # All files needed to create environment for reproducing analyses
└── manuscript             # Files for the manuscript, an R project and associated manuscript code
    ├── figures            # Figures for the manuscript
    └── manuscript_files   # Ancillary files for the manuscript
```

## Glossary

Throughout this repo and the manuscript itself we standardize the following terms and their definitions:

**_manuscript_**: The study from which a dataset or analysis is derived.

**_matrix_**: What we call a multiple sequence alignment that consists of one or more partitions. A manuscript can have one or more matrices.

**_gene_**: a set of homologous partitions. Globally consistent across matrices.

**_partition_**: a set of homologous sequences within a matrix, they consist of the same matrix columns

**_taxon_**: The name of a taxon as it appears in a matrix row name or tree tip. Usually but not always a species.

**_clade_**: a set of taxa, eg Cnidaria

**_sequence_**: a gene sequence in a specific partition for a specific taxon. It is a 1D character string, and a segment of a matrix row.

## git LFS

All the files in this repo that match the patterns defined in [`.gitattributes`](.gitattributes) are tracked with [git large file storage](https://git-lfs.github.com/). Install Git LFS following the instructions on the project's website to work with this repo.

Getting a fresh copy of the full repo would look something like this:

``` bash
git lfs install
git clone https://github.com/dunnlab/animal_tree_root.git
```

## Recreating full project

Because the original git repository for this project is quite large, this one pared down repo for general consumption. If you would like to recreate it, download the data from [Figshare](https://doi.org/10.6084/m9.figshare.13122881.v1). The data are split into three archived directories which can be downloaded separately or all together.

`data_raw.tar.xz` - The data from each previous study we used.
`reconciliation.tar.xz` - Files and scripts used to standardize naming and formats across the datasets used here.
`trees_new.tar.xz` - Results from the new analyses we did over the course of this study. Abbridged and summarized portiuons of these data are imported into the manuscript's R environment which is included in this repo as (`manuscript/manuscript.RData`)[manuscript/manuscript.RData]. You can examine it with the Docker environment, outlined in the (`docker`)[docker] directory.

### Download everything

``` bash
cd animal_tree_root
wget https://ndownloader.figshare.com/articles/13122881/versions/1 -O tmp.zip
unzip tmp.zip && rm tmp.zip

# this may take a bit
for t in *.tar.xz; do
echo "Expanding $t ..."
tar xf $t && rm $t
done
```

## Example Analyses

### IQ-TREE

We used [IQ-TREE 1.6.7](https://github.com/Cibiv/IQ-TREE/releases/tag/v1.6.7) for these analyses. To run the examples below as-is you will need the container image for this study set up. To do so, please see the [docker directory](docker/) and readme. Otherwise you may need to change your `iqtree` options/parameters depending on your version.

``` bash
# navigate to this repo on your computer
numthreads=4 # change this to suit your computer hardware / container settings

# you can use our docker image to run iqtree, or use your own by setting the variable iqtree=iqtree
# change /path/to/animal_tree_root to the location of this repo on your computer
iqtree="docker run --rm -it -w /animal_tree_root -v /path/to/animal_tree_root:/animal_tree_root animal_tree_root iqtree"
# iqtree=iqtree

# make a directory for output if it doesn't exist
mkdir -p examples_out

# GTR20
$iqtree -s data_processed/matrices/Philippe2009.phy -nt $numthreads -bb 1000 -m GTR20+F+G -pre examples_out/Philippe2009.GTR20 -wbt

# poisson_C60
$iqtree -s data_processed/matrices/Philippe2009.phy -nt $numthreads -bb 1000 -m Poisson+C60+F+G -pre examples_out/Philippe2009.poisson_C60 -wbt

# WAG
$iqtree -s data_processed/matrices/Philippe2009.phy -nt $numthreads -bb 1000 -m WAG+F+G -pre examples_out/Philippe2009.WAG -wbt

# Modelfinder
$iqtree -s data_processed/matrices/Philippe2009.phy -nt $numthreads -bb 1000 -mset LG,GTR20,WAG,Poisson -madd Poisson+C10+F+G,Poisson+C20+F+G,Poisson+C30+F+G,Poisson+C40+F+G,Poisson+C50+F+G,Poisson+C60+F+G,WAG+C10+F+G,WAG+C20+F+G,WAG+C30+F+G,WAG+C40+F+G,WAG+C50+F+G,WAG+C60+F+G,LG+C10+F+G,LG+C20+F+G,LG+C30+F+G,LG+C40+F+G,LG+C50+F+G,LG+C60+F+G -pre exampples_out/Philippe2009.model_test -wbt
```

### PhyloBayes MPI

We used [Phylobayes MPI](https://github.com/bayesiancook/pbmpi) compiled from commit [`01cbc7d`](https://github.com/bayesiancook/pbmpi/tree/01cbc7d9d9f192eb7be0e1dc7614169d444faa3d) in that repo on the Yale HPC cluster [Farnam](https://docs.ycrc.yale.edu/clusters-at-yale/clusters/farnam/) for these analyses. To run Phylobayes MPI you need MPI installed, and it makes most sense to run on an HPC cluster if you have one available. If you are not using slurm for job scheduling, you will need to change `srun` to `mpirun` and possibly pass some options/parameters to `mpirun`. Below are some example commands from our analyses of Philippe2009_only_choanozoa.phy

``` bash
# make a directory for output if it doesn't exist
mkdir -p examples_out

## PhyloBayes with poisson+CAT model
srun pb_mpi -cat -poisson -dgam 4 -s -d data_processed/matrices/Philippe2009_only_choanozoa.phy examples_out/Philippe2009_only_choanozoa.phy_Poisson_CAT_Chain_1
srun pb_mpi -cat -poisson -dgam 4 -s -d data_processed/matrices/Philippe2009_only_choanozoa.phy examples_out/Philippe2009_only_choanozoa.phy_Poisson_CAT_Chain_2

## PhyloBayes with GTR+CAT model
srun pb_mpi -cat -gtr -dgam 4 -s -d data_processed/matrices/Philippe2009_only_choanozoa.phy examples_out/Philippe2009_only_choanozoa.phy_GTR_CAT_Chain_1
srun pb_mpi -cat -gtr -dgam 4 -s -d data_processed/matrices/Philippe2009_only_choanozoa.phy examples_out/Philippe2009_only_choanozoa.phy_GTR_CAT_Chain_2

## PhyloBayes with poisson+nCAT60 model
srun pb_mpi -ncat 60 -poisson -dgam 4 -s -d data_processed/matrices/Philippe2009_only_choanozoa.phy examples_out/Philippe2009_only_choanozoa.phy_Poisson_CAT60_Chain_1
srun pb_mpi -ncat 60 -poisson -dgam 4 -s -d data_processed/matrices/Philippe2009_only_choanozoa.phy examples_out/Philippe2009_only_choanozoa.phy_Poisson_CAT60_Chain_2
```

Here is an example slurm submission script that would run 2 chains each across 2 cores for a maximum of 30 days.

``` bash
#!/bin/bash
#SBATCH -J phylobayes -p general
#SBATCH --ntasks=20 --mem-per-cpu 6G
#SBATCH --array=1-2
#SBATCH -t 30-00:00:00
#SBATCH --mail-type=ALL

module load PhyloBayes-MPI/20170808-foss-2016b
srun pb_mpi -cat -gtr -dgam 4 -s -d data_processed/matrices/Philippe2009_only_choanozoa.phy examples_out/Philippe2009_only_choanozoa.phy_GTR_CAT_Chain_${SLURM_ARRAY_TASK_ID}
```

## Citation

citation here
