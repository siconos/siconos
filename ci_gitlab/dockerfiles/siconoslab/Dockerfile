ARG REGISTRY=gricad-registry.univ-grenoble-alpes.fr
ARG PROJECT=nonsmooth/siconos
ARG DEFAULT_IMAGE=sources/jupyterlab
ARG IMAGENAME=$REGISTRY/$PROJECT/$DEFAULT_IMAGE
# - Install Siconos, from a previous build (CI artifacts from another job)
# Note: for local tests, use Dockerfile.manual to build the target siconos-build that can be used as IMAGENAME.
FROM $IMAGENAME AS siconoslab
ARG BUILD_DIR=build
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python, cmake ... will not be found)
RUN cmake --install $BUILD_DIR
# install Siconos. Assume a previous built in $BUILD_DIR
RUN siconos --info
# Clone Siconos examples 
ARG EXAMPLES_BRANCH_NAME=master
RUN git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials.git
RUN cd siconos-tutorials; git switch $EXAMPLES_BRANCH_NAME



