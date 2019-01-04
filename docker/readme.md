# How to reproduce this manuscript's application environment

## Docker

### Setup

Make sure you have [docker installed](https://docs.docker.com/install/#supported-platforms).

In this directory, build the container image with:

```
docker build -t animal_root .
```


### Running analyses

The model for executing the manuscript follows that at
https://github.com/caseywdunn/comparative_expression_2017/tree/master/docker .


To run an RStudio session:

```
docker run -e PASSWORD=secret123 -v /path/to/animal_root:/animal_root --rm -p 8787:8787 animal_root
```

Then go to http://localhost:8787/ in your browser. Log in with user/password rstudio/secret123 (Password can be anything except "rstudio").

To run an interactive shell:

```
docker run -it --rm -v /path/to/animal_root:/animal_root animal_root bash
```
