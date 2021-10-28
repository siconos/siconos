#!/bin/bash
mkdir -p tmp
sed "s|REGISTRY_PATH|nonsmooth/siconos|g" ci_gitlab/dockerfiles/$1/Dockerfile > tmp/Dockerfile1
sed "s|ci_gitlab/dockerfiles|.|g" ./tmp/Dockerfile1 > tmp/Dockerfile
cp -R ci_gitlab/dockerfiles/*.sh tmp
docker build -t gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos/sources/$1 ./tmp/
rm -rf tmp