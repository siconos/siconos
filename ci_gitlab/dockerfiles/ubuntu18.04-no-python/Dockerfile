FROM ubuntu:18.04
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
        gcc \
        gfortran \
        libgfortran3 \
        g++ \
        libopenblas-dev \
        liblapacke-dev \
        lp-solve \
        liblpsolve55-dev \
        python3 \
        python3-packaging \
        libpython3-dev \
        libcppunit-dev \
        libbullet-dev \
        libfreetype6-dev \
        freeglut3-dev \
        vim
RUN apt clean && apt autoremove
