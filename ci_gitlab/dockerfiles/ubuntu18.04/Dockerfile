FROM ubuntu:18.04
RUN apt update  && apt install -y -qq \
        cmake \
        git-core \
        make \
        libboost-dev \
        libboost-serialization-dev \
        libboost-filesystem-dev \
        libgmp-dev \
        swig \
        gcc \
        gfortran \
        g++ \
        liblapack-dev \
        libatlas-base-dev \
        lp-solve \
        liblpsolve55-dev \
        libpython3-dev \
        doxygen \
        libcppunit-dev \
        libbullet-dev \
        libfreetype6-dev \
        freeglut3-dev \
        python3-pip\
        python3-scipy \
        python3-pytest \
        valgrind \
        python3-lxml \
        python3-h5py
RUN apt clean && apt autoremove
