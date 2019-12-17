FROM gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos/ubuntu18.04-no-python
RUN apt update  && apt install -y -qq \
        swig \
        doxygen \
        python3-pip\
        python3-scipy \
        python3-pytest \
	valgrind \
        python3-lxml \
        python3-h5py
RUN apt clean && apt autoremove
