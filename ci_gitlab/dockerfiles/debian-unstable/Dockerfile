FROM debian:unstable
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
        vim

# Serialization + Generation build requires documentation packages
RUN apt install -y -qq \
        python3-breathe \
        python3-numpydoc \
        python3-sphinxcontrib.bibtex \
        python3-sphinxcontrib.youtube \
        python3-sphinxcontrib.websupport \
        python3-sphinx-rtd-theme \
        python3-sphinx-bootstrap-theme

RUN apt clean && apt autoremove && rm -rf /var/lib/apt/lists/*
