FROM debian:bookworm AS siconos-base
ENV DEBIAN_FRONTEND=noninteractive
RUN apt update  && apt install  --no-install-recommends -y -qq \
        git-core \
	cmake \
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
        python3-full\
        valgrind \
        pkg-config \
	libxrender1 \
	libbullet-dev \
	libbullet-extras-dev \
	libeigen3-dev \
	libeigen3-doc \
	wget \
        vim && apt clean && apt autoremove && rm -rf /var/lib/apt/lists/*
RUN useradd -ms /bin/bash siconos-user
RUN wget https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.bz2 && \
    tar -jxvf boost_1_81_0.tar.bz2 && \
    mv boost_1_81_0/boost /usr/local/include

FROM siconos-base AS siconos-base-python
USER siconos-user
WORKDIR /home/siconos-user
COPY ci_gitlab/dockerfiles/requirements.txt .
ENV VIRTUAL_ENV=/home/siconos-user/siconosenv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN pip install --no-cache-dir -U -r requirements.txt

