# REGISTRY_PATH will be updated (sed) during kaniko call
FROM gricad-registry.univ-grenoble-alpes.fr/REGISTRY_PATH/sources/ubuntu20.04-light
# COPY ci_gitlab/dockerfiles/install_fclib.sh /home
RUN apt update  && apt upgrade -y &&  apt install -y -qq \
        lp-solve \
        liblpsolve55-dev \
        libhdf5-dev \
        libboost-serialization-dev \
        libfreetype6-dev \
        freeglut3-dev \
        swig \
 	libxrender1 \
        libpython3-dev \
        python3 \
        python3-pip \
        python3-scipy \
        python3-pytest \
	valgrind \
        python3-lxml \
        python3-packaging \
        python3-h5py && apt autoclean -y && apt autoremove -y&& rm -rf /var/lib/apt/lists/*
# RUN sh /home/install_fclib.sh
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install -U vtk
WORKDIR /home
COPY ci_gitlab/dockerfiles/install_bullet.sh .
ENV CI_PROJECT_DIR /home
RUN sh /home/install_bullet.sh
RUN apt clean && apt autoremove && rm -rf /var/lib/apt/lists/*

