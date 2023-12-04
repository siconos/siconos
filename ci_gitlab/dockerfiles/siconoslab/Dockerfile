FROM IMAGENAME
ARG BUILD_DIR=build
USER root
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python, cmake ... will not be found)
RUN cmake --install $BUILD_DIR
USER $NB_USER
ENV SICONOS_INSTALL_DIR=/home/install-siconos PATH=/home/install-siconos/bin/:$PATH
ENV JUPYTER_ENABLE_LAB yes
RUN cd $HOME; git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials.git


