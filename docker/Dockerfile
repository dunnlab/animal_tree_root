FROM rocker/geospatial:3.5.3

# can override these at build time
# e.g. docker build --build-arg arg_iqtreeversion=1.6.8

ARG arg_iqtreeversion=1.6.7
ARG arg_minicondaversion=4.5.12
ARG arg_gitlfsversion=2.6.0
ARG arg_hutan_ref=b5bebb4
ARG arg_bioc_version=3.8

ENV iqtreeversion=${arg_iqtreeversion}
ENV minicondaversion=${arg_minicondaversion}
ENV gitlfsversion=${arg_gitlfsversion}
ENV condaprefix=/opt/conda
ENV CONDA_ENVS_PATH=${condaprefix}/envs
ENV bioc_version=${arg_bioc_version}
ENV hutan_ref=${arg_hutan_ref}

# install git LFS
RUN wget --quiet https://github.com/git-lfs/git-lfs/releases/download/v${gitlfsversion}/git-lfs-linux-amd64-v${gitlfsversion}.tar.gz && \
    tar xvf git-lfs-linux-amd64-v${gitlfsversion}.tar.gz \
        -C /usr/local/bin/ \
        git-lfs && \
    rm git-lfs-linux-amd64-v${gitlfsversion}.tar.gz

# install iqtree
RUN wget --quiet https://github.com/Cibiv/IQ-TREE/releases/download/v${iqtreeversion}/iqtree-${iqtreeversion}-Linux.tar.gz && \
    tar xvf iqtree-${iqtreeversion}-Linux.tar.gz \
        --strip-components=2 \
        -C /usr/local/bin/ \
        iqtree-${iqtreeversion}-Linux/bin/iqtree && \
    rm iqtree-${iqtreeversion}-Linux.tar.gz

# install miniconda
RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda3-${minicondaversion}-Linux-x86_64.sh && \
    /bin/bash Miniconda3-${minicondaversion}-Linux-x86_64.sh -f -b -p ${condaprefix} && \
    rm Miniconda3-${minicondaversion}-Linux-x86_64.sh

ENV PATH=${condaprefix}/bin:$PATH

COPY conda_animal_tree_root.yml ${condaprefix}

RUN conda env create -f ${condaprefix}/conda_animal_tree_root.yml

RUN install2.r --error \
    --deps TRUE \
    ape \
    BiocManager \
    doParallel \
    ggpubr \
    igraph \
    R.utils \
    viridis

# need to get older version of quadprog
# Error: package ‘quadprog’ is not available (for R version 3.5.3)
RUN wget -O /tmp/quadprog_1.5-5.tar.gz \
         https://cran.r-project.org/src/contrib/Archive/quadprog/quadprog_1.5-5.tar.gz && \
    R CMD INSTALL /tmp/quadprog_1.5-5.tar.gz

RUN Rscript -e "library(BiocManager); BiocManager::install('Biostrings', version = $bioc_version); \
                                      BiocManager::install('ggtree', version = $bioc_version); \
                                      BiocManager::install('treeio', version = $bioc_version);"

RUN Rscript -e "library(devtools); devtools::install_github('caseywdunn/hutan@$hutan_ref')"
