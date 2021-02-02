# Continuous integration setup
[[_TOC_]]

## Driving the continuous integration process

The behavior of the CI (which jobs are created and launched) depends on the content of the last commit message and on the targeted branch as summarized in the table below.

Job type are:

* [A] Those to **Prepare docker images**: jobs used to build and save docker images with all the dependencies required to install siconos.

    * Docker images are saved in the project registry, in 'sources',  see https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/container_registry.
    * The images are built from Dockerfiles saved in [ci_gitlab/dockerfiles](./ci_gitlab/dockerfiles)

* [B] Those to configure, build or test Siconos 

    * They are executed on one of the images saved in the registry.
    * Build and tests results are published on [siconos-dashboard](http://siconos-dashboard.univ-grenoble-alpes.fr:8080/index.php?project=siconos).

* [C] Those to generate a docker image with a fully functionnal install of Siconos.

    * Based on the results of a previous job of type [B]
    * The image is saved in the registry, https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/
    * The image is usually named siconos-<CI_COMMIT_REF_NAME>-<osname>, eg siconos-master-ubuntu20


| Commit starts with ...   |  [skip ci]  | [docker-build]                          | [all-jobs] | any other message |       
| ---                      |  ------     |----------------                         |---------------------------------------|-------------------|
| push to master           |   :x:       | :white_check_mark: all [A]<br>:white_check_mark: all [B] (debian, ubuntu ...)<br>:white_check_mark: [C] on ubuntu20.04 |useless| :white_check_mark: all [A] <br>:white_check_mark: [C] on ubuntu20.04 |
| push to any other branch |   :x:       | :white_check_mark: all [A]<br>:white_check_mark: [B] on ubuntu20.04<br>:white_check_mark: [C] on ubuntu20.04 | :white_check_mark: all [A] (debian, ubuntu ...)<br>:white_check_mark: [C] on ubuntu20.04 | :white_check_mark: [B] on ubuntu20.04<br>:white_check_mark: [C] on ubuntu20.04|


Moreover, some jobs are optional and must be started directly by clicking on the little gear next to the job name on gitlab CI interface (e.g. siconos with oce, documentation ...). 

They look like ![manual_ci.png](./manual_ci.png)

To check the complete pipeline, see https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/-/pipelines.



## CI configuration

The directory [ci_gitlab](./) contains cmake/ctest and script files used in .gitlab-ci.yml
for CI on gricad-gitlab, with:

* [ctest_siconos.sh](./ctest_siconos.sh) is the file call in CI jobs to run ctest (from configure to tests)
* [siconos_conf](./siconos_conf): contains setup files to be used to configure siconos in CI jobs,
as in

```
cmake -DUSER_OPTIONS_FILE=<path-to-siconos>/ci_gitlab/siconos_confs/siconos_no_python.cmake
```

  
* [dockerfiles](./dockerfiles): contains directories (with Dockerfiles and some other
required files) used to build Docker images "ready-to-use" for Siconos install.

Examples of usage :

```
# build a docker image named some_name with all required deps for Siconos
docker build -t <some_name> dockerfiles/ubuntu18.04-fclib
```


```
start a docker container based on this image.
docker run -ti some_name /bin/bash
```

All these dockerfiles directories are used on gricad-gitlab to build images used in CI process for Siconos.

The list of already built images is here:
https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/container_registry

To download and run an image, try for example:

```
docker login gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos
docker run -ti gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos/ubuntu18.04 /bin/bash
```

## Predefined environment variables

### gitlab-ci

https://docs.gitlab.com/ee/ci/variables/predefined_variables.html

### Travis

https://docs.travis-ci.com/user/environment-variables/#default-environment-variables

