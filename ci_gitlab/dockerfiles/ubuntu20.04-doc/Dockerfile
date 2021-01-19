FROM gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos/sources/ubuntu20.04
RUN apt update  && apt upgrade -y && apt install -y -qq \
        graphviz doxygen && apt autoclean -y && apt autoremove -y&& rm -rf /var/lib/apt/lists/*
WORKDIR /home
COPY ci_gitlab/dockerfiles/requirements4doc.txt /home
RUN pip3 install -U -r /home/requirements4doc.txt
ENV LANG C.UTF-8 # Required, else doxy2swig fails!
