# REGISTRY_PATH will be updated (sed) during kaniko call
FROM gricad-registry.univ-grenoble-alpes.fr/REGISTRY_PATH/ubuntu18.04
RUN apt update  && apt install -y -qq \
        liboce-foundation-dev \
        liboce-modeling-dev \
        liboce-ocaf-dev \
        liboce-visualization-dev
WORKDIR /home
COPY ci_gitlab/dockerfiles/install_oce.sh .
ENV CI_PROJECT_DIR /home
RUN sh /home/install_oce.sh
ENV PYTHONPATH /home/install/site-packages
RUN apt clean && apt autoremove && rm -rf /var/lib/apt/lists/*
