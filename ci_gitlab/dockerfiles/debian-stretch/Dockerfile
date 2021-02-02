FROM debian:stretch
RUN apt update  && apt install -y -qq \
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
        libgfortran3 \
        g++ \
        libopenblas-dev \
        liblapacke-dev \
        lp-solve \
        libsuitesparse-dev \
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
        python3-packaging \
        valgrind \
        python3-lxml \
        python3-h5py && apt clean && apt autoremove && rm -rf /var/lib/apt/lists/*
