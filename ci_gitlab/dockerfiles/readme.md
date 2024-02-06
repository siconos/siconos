# Docker

This directory contains dockerfiles and other configuration files used to build docker images with all the libraries and tools
required to configure and build siconos.

Sub-directories are named according to the standard image used as a base to prepare the image and corresponds to the targeted operating system.




For instance, to create a new image, named sico-bulls, based on Debian bookworm version, run 


```
docker build -t sico-bulls -f $CI_PROJECT_DIR/ci_gitlab/dockerfiles/debian-bookworm/Dockerfile $CI_PROJECT_DIR
```

- Image name, sico-bulls here, can be set to any value you want.
- CI_PROJECT_DIR environment variable must corresponds to the absolute path to siconos repository.


And then, to start a container based on this new image:


```
docker run -ti sico-bulls
```

If required, you can mount your siconos repository inside the new container:

```
docker run -v $CI_PROJECT_DIR:/home/siconos -ti sico-bulls
```

You will find your siconos repository in /home/siconos inside the container.
