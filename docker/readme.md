# How to reproduce this manuscript's application environment

## Docker

### Setup

Make sure you have [docker installed](https://docs.docker.com/install/#supported-platforms).

In this directory, build the container image with:

``` bash
docker build -t animal_tree_root .
```

### Running analyses

The model for executing the manuscript follows that at
https://github.com/caseywdunn/comparative_expression_2017/tree/master/docker .

To run an RStudio session:

``` bash
docker run --rm  -dP -e PASSWORD=secret123 -e USERID=$UID -v /path/to/animal_tree_root:/animal_tree_root -p 8787:8787 animal_tree_root
```

Then go to [http://localhost:8787/](http://localhost:8787/) in your browser. Log in with user/password rstudio/secret123 (Password can be anything except "rstudio"). You can jump to the repo's R directory with the R command:

``` R
setwd("/animal_tree_root/manuscript")
```

To start an interactive shell for running other analyses (e.g. those in `/animal_tree_root/reconciliation/scripts`):

``` bash
docker run --rm -it -e USERID=$UID -v /path/to/animal_tree_root:/animal_tree_root animal_tree_root bash
```
