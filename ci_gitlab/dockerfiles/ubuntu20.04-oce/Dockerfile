# REGISTRY_PATH will be updated (sed) during kaniko call
FROM gricad-registry.univ-grenoble-alpes.fr/REGISTRY_PATH/sources/ubuntu20.04
RUN apt update  && apt upgrade -y && apt install -y -qq \
        liboce-foundation-dev \
        liboce-modeling-dev \
        liboce-ocaf-dev \
        liboce-visualization-dev && apt autoclean -y && apt autoremove -y&& rm -rf /var/lib/apt/lists/*
WORKDIR /home
COPY ci_gitlab/dockerfiles/install_oce.sh /home
ENV CI_PROJECT_DIR /home 
RUN sh /home/install_oce.sh
ENV PYTHONPATH /home/install/site-packages
