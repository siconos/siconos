FROM gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos/ubuntu18.04
RUN apt update  && apt install -y -qq \
        liboce-foundation-dev \
        liboce-modeling-dev \
        liboce-ocaf-dev \
        liboce-visualization-dev
WORKDIR /home
COPY install_oce.sh .
ENV CI_PROJECT_DIR /home
RUN sh /home/install_oce.sh
RUN apt clean && apt autoremove
