FROM ubuntu:18.04
RUN apt update  && apt install -y -qq \
        cmake \
        git-core \
        make \
        libboost-dev \
        libgmp-dev \
        gcc \
        gfortran \
        g++ \
        liblapack-dev \
        libatlas-base-dev \
        lp-solve \
        liblpsolve55-dev \
        python3 \
        libpython3-dev \
        libcppunit-dev \
        libbullet-dev \
        libfreetype6-dev \
        freeglut3-dev
RUN apt clean && apt autoremove
