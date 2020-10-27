# Animal tree root

A comparison of phylogenetic studies relevant to placing the root of the animal phylogeny.

## Repo Overview

``` text
.
├── data_processed         # Matrices and data tables used for analyses in this manuscript
│   ├── matrices           # Matrices in consistent formats with harmonized taxon names
│   └── tables             # Tabular summaries of previously published datasets and results
├── docker                 # All files needed to create environment for reproducing analyses
└── manuscript             # Files for the manuscript, an R project and associated manuscript code
    ├── figures            # Figures for the manuscript
    └── manuscript_files   # Ancillary files for the manuscript
```

The matrices we curated, standardized, and started from for all new phylogenetic analyses are stored in [`data_processed/matrices`](data_processed/matrices/), eponymously named for their original manuscript. Tables containing summaries, taxon-clade maps, partition-gene maps, etc. are all in [`data_processed/tables`](data_processed/tables/). We maintained a frozen set of applications for the project as a Docker container image, defined in the [`docker`](docker) directory. See the readme there for instructions for building and running an RStudio Server that is compatible with the manuscript. The manuscript is stored in the [`manuscript`](manuscript) directory as an [R Markdown file](manuscript/manuscript.rmd), and the data we use in visualization and summary are all stored in an [RData file](manuscript/manuscript.RData) to avoid needing to re-run some [time-consuming functions](manuscript/manuscript_kernel.R).

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

Because the original git repository for this project is quite large, this one is distilled to just what is needed for

- Launching new analyses based on the standardized matrices in [`data_processed/matrices`](data_processed/matrices/)
- Examinging the data used in the manuscript R analyses and figures

If you would like to recreate the full repository that includes raw output from our analyses, download the data from [Figshare](https://doi.org/10.6084/m9.figshare.13122881.v1). The data are split into three archived directories which can be downloaded separately or all together.

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

We used [IQ-TREE 1.6.7](https://github.com/Cibiv/IQ-TREE/releases/tag/v1.6.7) for these analyses. The latest versions are available from the [IQ-TREE Downloads page](http://www.iqtree.org/#download). You may need to change your `iqtree` options/parameters depending on your version - some options changed in v2.x.

``` bash
# navigate to this repo on your computer
# cd ~/repos/animal_tree_root

# make a directory for output if it doesn't exist
mkdir -p examples_out

# run modelfinder for the Philippe2009 matrix
iqtree -s data_processed/matrices/Philippe2009_only_choanozoa.phy -nt AUTO -bb 1000 -o Monosiga_ovata -mset LG,GTR20,WAG,Poisson -madd Poisson+C10+F+G,Poisson+C20+F+G,Poisson+C30+F+G,Poisson+C40+F+G,Poisson+C50+F+G,Poisson+C60+F+G,WAG+C10+F+G,WAG+C20+F+G,WAG+C30+F+G,WAG+C40+F+G,WAG+C50+F+G,WAG+C60+F+G,LG+C10+F+G,LG+C20+F+G,LG+C30+F+G,LG+C40+F+G,LG+C50+F+G,LG+C60+F+G -pre examples_out/Philippe2009.model_test -wbt
```

### PhyloBayes MPI

We used [Phylobayes MPI](https://github.com/bayesiancook/pbmpi) compiled from commit [`01cbc7d`](https://github.com/bayesiancook/pbmpi/tree/01cbc7d9d9f192eb7be0e1dc7614169d444faa3d) in that repo on the Yale HPC cluster [Farnam](https://docs.ycrc.yale.edu/clusters-at-yale/clusters/farnam/) for these analyses. To run Phylobayes MPI you need MPI installed, and it makes most sense to run on an HPC cluster if you have one available. If you are not using slurm for job scheduling, you will need to change `srun` to `mpirun` and possibly pass some options/parameters to `mpirun`. Below are some example commands from our analyses of Philippe2009_only_choanozoa.phy

Here is an example slurm submission script that would run `pbmpi` on the `Philippe2009_only_choanozoa` matrix with the GTR+CAT model. The job scheduler will run 2 chains as separate jobs, each across 20 cores for a maximum of 30 days.

``` bash
#!/bin/bash
#SBATCH -J phylobayes -p general
#SBATCH --ntasks=20 --mem-per-cpu 6G
#SBATCH --array=1-2
#SBATCH -t 30-00:00:00
#SBATCH --mail-type=ALL

module load PhyloBayes-MPI/20170808-foss-2016b
# GTR+CAT model
srun pb_mpi -cat -gtr -dgam 4 -s -d  Philippe2009_only_choanozoa.phy Philippe2009_only_choanozoa.phy_GTR_CAT_Chain_${SLURM_ARRAY_TASK_ID}
```

## Citation

### Pre-print



### Additional Datasets

Li, Yuanning; Shen, Xing-Xing; Evans, Benjamin; W. Dunn, Casey; Rokas, Antonis (2020): Rooting the animal tree of life. figshare. Dataset. https://doi.org/10.6084/m9.figshare.13122881.v1
