# How to reproduce this manuscript's application environment

## Docker

To preserve the application environment (e.g. Python, R, specific package versions) used in this project, we define it as a Docker container image. By also using this container definition (on an x86_64 CPU) you can replicate our environment.

### Setup

Install [Docker](https://docs.docker.com/install/#supported-platforms) if you have not already.

Build the container image in this directory with:

``` bash
cd animal_tree_root/docker
docker build -t animal_tree_root .
```

Docker will automatically look for a recipe file called `Dockerfile` in this directory and build the container image for you, tagging it as `animal_tree_root`.

### Running analyses

The manuscript was written in R Markdown. You can start an RStudio Server session to examine it and its environment:

``` bash
# change /path/to/animal_tree_root to the location of this repo on your computer
docker run --rm  -dP -e PASSWORD=secret123 -e USERID=$UID -v /path/to/animal_tree_root:/animal_tree_root -p 8787:8787 animal_tree_root
```

Then go to [http://localhost:8787/](http://localhost:8787/) in your browser. Log in with username: `rstudio` and password: `secret123` . You can jump to the repo's R directory with the R command:

``` R
setwd("/animal_tree_root/manuscript")
```

To start an interactive shell for running other analyses (e.g. those in `/animal_tree_root/reconciliation/scripts`):

``` bash
# change /path/to/animal_tree_root to the location of this repo on your computer
docker run --rm -it -e USERID=$UID -v /path/to/animal_tree_root:/animal_tree_root animal_tree_root bash
```
