FROM debian:bookworm
RUN apt update  && apt install -y -qq \
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
        vim && apt clean && apt autoremove && rm -rf /var/lib/apt/lists/*
COPY ci_gitlab/dockerfiles/requirements.txt /home
ENV VIRTUAL_ENV=/home/siconosenv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN pip install -U -r /home/requirements.txt
WORKDIR /home
ENV CI_PROJECT_DIR /home

