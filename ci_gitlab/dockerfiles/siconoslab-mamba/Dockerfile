ARG REGISTRY=gricad-registry.univ-grenoble-alpes.fr
ARG PROJECT=nonsmooth/siconos
ARG DEFAULT_IMAGE=sources/mambaenv
ARG IMAGENAME=$REGISTRY/$PROJECT/$DEFAULT_IMAGE
FROM $IMAGENAME AS siconoslab
ARG BUILD_DIR=build
# ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python, cmake ... will not be found)
RUN cmake --install $BUILD_DIR
# install Siconos. Assume a previous built in $BUILD_DIR
RUN siconos --info
# Clone Siconos examples 
ARG EXAMPLES_BRANCH_NAME=master
RUN git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials.git
RUN cd siconos-tutorials; git switch $EXAMPLES_BRANCH_NAME
ENV JUPYTER_PORT=8888
ENV JUPYTER_ENABLE_LAB yes
EXPOSE $JUPYTER_PORT
#ENV SHELL /bin/bash	
#ENV PATH /opt/conda/bin:$PATH
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/opt/conda/bin/jupyter", "lab","--ip=0.0.0.0","--allow-root", "--no-browser"]
# WORKDIR /home/mambauser


