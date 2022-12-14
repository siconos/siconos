FROM ubuntu:18.04
ENV TZ=Europe/Paris
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && apt update  && apt install -y -qq \
        ntp \
        cmake \
        git-core \
        make \
        libboost-dev \
        libboost-filesystem-dev \
        libboost-timer-dev \
        libboost-chrono-dev \
        libgmp-dev \
        gcc \
        gfortran \
        libgfortran-8-dev \
        g++ \
        libopenblas-dev \
        liblapacke-dev \
        libcppunit-dev \
        lp-solve \
        liblpsolve55-dev \
        vim  \
        swig \
        doxygen \
        python3-pip\
        python3-scipy \
        python3-pytest \
	valgrind \
        python3-lxml \
        python3-h5py \
        python3-packaging \
        libopenmpi-dev \
        libparmetis-dev \
        libptscotch-dev \
        openssh-client \
        libbullet-dev \
        libfreetype6-dev \
        freeglut3-dev \
        libmumps-dev && apt autoclean -y && apt autoremove -y && rm -rf /var/lib/apt/lists/*
RUN python3 -m pip install mpi4py vtk

