# Animal tree root

A comparison of phylogenetic studies relevant to placing the root of the animal phylogeny.

## Repo Overview

(See each directory for more documentation)

``` text
.
├── data_processed         # Data processed for analyses in this manuscript
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

**_sequence_**: a gene sequence in a specific partition for a specific taxon. It is a 1D character string, and a segment 
of a matrix row.

## git LFS

All the files in this repo that match the patterns defined in [`.gitattributes`](.gitattributes) are tracked with [git large file storage](https://git-lfs.github.com/). Install Git LFS following the instructions on the project's website to work with this repo.

Getting a fresh copy of the full repo would look something like this:

``` bash
git lfs install
git clone https://github.com/dunnlab/animal_tree_root.git
```

## Recreating full project

Because the original repo for this project is quite large, we present a pared down repo for general consumption. If you would like to recreate it, download the data from Figshare:

```
TODO
```

## Citation

citation here
