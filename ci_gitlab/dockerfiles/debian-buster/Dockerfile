FROM debian:buster
RUN apt update  && apt upgrade -y && apt install -y -qq \
        cmake \
        git-core \
        make \
        libboost-dev \
        libboost-serialization-dev \
        libboost-filesystem-dev \
        libboost-timer-dev \
        libboost-chrono-dev \
        libgmp-dev \
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
        libbullet-dev \
        libfreetype6-dev \
        freeglut3-dev \
        python3-pip\
        python3-scipy \
        python3-pytest \
        python3-packaging \
        valgrind \
        python3-lxml \
        python3-h5py \
        vim && apt autoclean -y && apt autoremove -y && rm -rf /var/lib/apt/lists/*
