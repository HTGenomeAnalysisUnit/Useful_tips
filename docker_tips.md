# Tips and tricks for Docker containers

## Avoid debian question

Include this at the beginning of your Dockerfile:
`ARG DEBIAN_FRONTEND=noninteractive`

## Apt-get

When using `apt-get` to install some packages remember to include `-y` option and clean the cache. Something like the following command:

```Dockerfile
RUN apt-get install -y your_packages \
    && apt-get -y clean all \
    && rm -rf /var/cache
```

**NB.** Avoid the last line if you mean to re-use this image with `FROM` since this may break installation of additional packages.

## Use python packages with --no-home

If you need to run python packages in the image and you plan to use the image in Singularity with `--no-home` (home directory is not mounted) then remember to set the following env vars to a writable location:

```Dockerfile
ENV NUMBA_CACHE_DIR=/tmp
ENV MPLCONFIGDIR=/tmp/matplotlib
```

## Use R packages

When you run an R environment in a container it would be better to use the image in Singularity with `--no-home --cleanenv` (home directory is not mounted and no env variables are passed from host). This is needed to ensure there are no conflicts with R packages installed in the host. It is also a good idea to set R_LIBS path to some custom location:

```Dockerfile
ENV R_LIBS /opt/R/libs
ENV R_LIBS_USER /opt/R/user_libs
```

## Reticulate

If you need to use `reticulate` in a R script in your container and you are not using conda, it is better to set the precise location of the python bin to avoid issues:

```Dockerfile
ENV RETICULATE_PYTHON=/usr/bin/python3
```

## Conda

To use conda in your container, you first need to install miniconda using the following command:

```Dockerfile
#  Install miniconda
RUN apt-get update && apt-get install -y wget
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
  /bin/bash ~/miniconda.sh -b -p /opt/conda
ENV PATH=/opt/conda/bin:${PATH}

RUN conda update -y conda
```

You can then install any package you need in the root env using something like this:

```Dockerfile
#Using a yml file
RUN conda env update -n root -f environment.yml

#Install manually
RUN conda install -n root -c bioconda your_packages
```
