FROM debian:buster
RUN apt update  && apt install -y -qq \
        git-core \
        make \
        libboost-dev \
        libboost-serialization-dev \
        libboost-filesystem-dev \
        libboost-timer-dev \
        libboost-chrono-dev \
        libgmp-dev \
        libhdf5-mpi-dev \
        swig \
        gcc \
        gfortran \
        g++ \
        libopenblas-dev \
        liblapacke-dev \
        lp-solve \
        liblpsolve55-dev \
        libpython3-dev \
        libsuitesparse-dev \
        doxygen \
        libcppunit-dev \
        libfreetype6-dev \
        freeglut3-dev \
        python3-pip\
        python3-scipy \
        python3-pytest \
        python3-packaging \
        valgrind \
        python3-lxml \
        python3-h5py \
        pkg-config \
	libxrender1	
        vim && apt clean && apt autoremove && rm -rf /var/lib/apt/lists/*
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install cmake
RUN python3 -m pip install vtk
WORKDIR /home
COPY ci_gitlab/dockerfiles/install_bullet.sh .
ENV CI_PROJECT_DIR /home
RUN sh /home/install_bullet.sh

