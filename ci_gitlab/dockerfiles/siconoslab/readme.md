# Dockerfiles to produce 'Siconos-ready' images 

images including 
- a complete Siconos install
- the siconos-tutorials repository
- jupyter notebook install


## Files

- Dockerfile: for the CI on Gitlab

- Dockerfile.manual: a local version to be able to build images locally, without gitlab-ci

## Usage


###  build, install Siconos AND prepare a notebook with siconos-tutorials included


```bash
docker build -f ci_gitlab/dockerfiles/siconoslab/Dockerfile.manual . --target=siconoslab -t my-siconoslab
```

Requires: the siconos repository must be in . and will be used as sources for install.

```
... --build-arg OPT=VALUE  --build-arg OPT2=VALUE2 ...
```

- IMAGENAME: to change source image. Must be an image will all required deps for Siconos and with jupyterlab.
  e.g.: gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos/sources/jupyterlab (default) or gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos/sources/jupyterlab-mamba
- CONF_FILE: config file used to run cmake.
- BUILD_DIR: place to build Siconos

Then run

```bash
docker run --rm -p 8888:8888 -ti siconoslab
```
