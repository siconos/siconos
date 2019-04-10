FROM gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos/ubuntu18.04
RUN apt update  && apt install -y -qq \
    libhdf5-dev
RUN apt clean && apt autoremove
WORKDIR /home
COPY install_fclib.sh .    
RUN sh /home/install_fclib.sh