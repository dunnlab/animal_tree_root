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

[`data_raw.tar.xz`](https://ndownloader.figshare.com/files/25186634): The data from each previous study we used.

[`reconciliation.tar.xz`](https://ndownloader.figshare.com/files/25186628): Files and scripts used to standardize naming and formats across the datasets used here.

[`trees_new.tar.xz`](https://ndownloader.figshare.com/files/25186655): Results from the new analyses we did over the course of this study. Abbridged and summarized portiuons of these data are imported into the manuscript's R environment which is included in this repo as (`manuscript/manuscript.RData`)[manuscript/manuscript.RData]. You can examine it with the Docker environment, outlined in the (`docker`)[docker] directory.

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

## Glossary

Throughout this repo and the manuscript itself we standardize the following terms and their definitions:

**_manuscript_**: The study from which a dataset or analysis is derived.

**_matrix_**: What we call a multiple sequence alignment that consists of one or more partitions. A manuscript can have one or more matrices.

**_gene_**: a set of homologous partitions. Globally consistent across matrices.

**_partition_**: a set of homologous sequences within a matrix, they consist of the same matrix columns

**_taxon_**: The name of a taxon as it appears in a matrix row name or tree tip. Usually but not always a species.

**_clade_**: a set of taxa, eg Cnidaria

**_sequence_**: a gene sequence in a specific partition for a specific taxon. It is a 1D character string, and a segment of a matrix row.


## Citation

### Pre-print



### Additional Datasets

Li, Yuanning; Shen, Xing-Xing; Evans, Benjamin; W. Dunn, Casey; Rokas, Antonis (2020): Rooting the animal tree of life. figshare. Dataset. https://doi.org/10.6084/m9.figshare.13122881.v1
