# Continuous integration setup

This directory contains cmake/ctest and script files used in .gitlab-ci.yml
for CI on gricad-gitlab.


* siconos_conf: this directory contains some examples of setup files to be used to configure siconos.

e.g. :

cmake -DUSER_OPTIONS_FILE=<path-to-siconos>/ci_gitlab/siconos_confs/siconos_no_python.cmake


  
* dockerfiles: contains directories (with Dockerfiles and some other
required files) used to build Docker images "ready-to-use" for Siconos install.


examples of usage :

docker build -t <some_name> dockerfiles/ubuntu18.04-fclib

--> build a docker image named some_name with all required deps for Siconos

docker run -ti some_name /bin/bash

--> start a docker container based on this image.



All these dockerfiles directories are used on gricad-gitlab to build
some images used in CI process for Siconos.

The list of already built images is here:
https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/container_registry

To download and run an image, try for example:

docker run -ti gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos/ubuntu18.04 /bin/bash

## Predefined environment variables

### gitlab-ci

https://docs.gitlab.com/ee/ci/variables/predefined_variables.html

### Travis

https://docs.travis-ci.com/user/environment-variables/#default-environment-variables

