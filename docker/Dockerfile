FROM rocker/tidyverse:3.5.1


RUN apt-get update \
  && apt-get install -y --no-install-recommends \
  libfftw3-dev \
  libtiff5-dev \
  libglu1-mesa-dev \
  lmodern \
  procps \
  texlive-latex-extra \
  texlive-fonts-recommended \
  vim


RUN install2.r --error \
    --deps TRUE \
    -r 'http://cran.rstudio.com' \
    -r 'http://www.bioconductor.org/packages/release/bioc' \
    ggtree \
    treeio \
    ape \
    doParallel \
    ggrepel \
    igraph \
    viridis \
  && Rscript -e "library(devtools); devtools::install_github('caseywdunn/hutan')"




# Install dependencies for python etc

# can override these at build time
# e.g. docker build --build-arg arg_iqtreeversion=1.6.8

ARG arg_iqtreeversion=1.6.7
ARG arg_minicondaversion=4.5.11
ARG arg_gitlfsversion=2.6.0

ENV iqtreeversion=${arg_iqtreeversion}
ENV minicondaversion=${arg_minicondaversion}
ENV gitlfsversion=${arg_gitlfsversion}
ENV condaprefix=/opt/conda
ENV CONDA_ENVS_PATH=${condaprefix}/envs

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
RUN apt-get -y install bzip2 && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-${minicondaversion}-Linux-x86_64.sh && \
    /bin/bash Miniconda3-${minicondaversion}-Linux-x86_64.sh -f -b -p ${condaprefix} && \
    rm Miniconda3-${minicondaversion}-Linux-x86_64.sh

ENV PATH=${condaprefix}/bin:$PATH

COPY conda_animal_root.yml ${condaprefix}

RUN conda env create -f ${condaprefix}/conda_animal_root.yml
