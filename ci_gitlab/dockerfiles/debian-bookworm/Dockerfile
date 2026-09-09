FROM debian:bookworm@sha256:9344f8b8992482f80cba753f323adeaf17690076c095ccff6cc9536be98185dc AS siconos-base
ENV DEBIAN_FRONTEND=noninteractive
ARG DEBIAN_SNAPSHOT=20260731T022227Z
RUN echo "deb [check-valid-until=no] http://snapshot.debian.org/archive/debian/${DEBIAN_SNAPSHOT} bookworm main" > /etc/apt/sources.list \
        && echo "deb [check-valid-until=no] http://snapshot.debian.org/archive/debian-security/${DEBIAN_SNAPSHOT} bookworm-security main" >> /etc/apt/sources.list
RUN apt update  && apt install  --no-install-recommends -y -qq \
        git-core \
        cmake \
        make \
        #libboost-dev \
        #libboost-serialization-dev \
        #libboost-filesystem-dev \
        #libboost-timer-dev \
        #libboost-chrono-dev \
        libgmp-dev \
        libhdf5-mpi-dev \
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
        nlohmann-json3-dev \
        valgrind \
        pkg-config \
        libxrender1 \
        libbullet-dev \
        libbullet-extras-dev \
        libeigen3-dev \
        libeigen3-doc \
        wget \
        bzip2 \
        vim && apt clean && apt autoremove && rm -rf /var/lib/apt/lists/*
RUN useradd -ms /bin/bash siconos-user
RUN wget https://archives.boost.io/release/1.81.0/source/boost_1_81_0.tar.bz2 && \
        tar -jxvf boost_1_81_0.tar.bz2 && \
        cd boost_1_81_0 && ./bootstrap.sh && ./b2 install

FROM siconos-base AS siconos-base-python
USER siconos-user
WORKDIR /home/siconos-user
COPY ci_gitlab/dockerfiles/requirements-pinned-2026-07-31.txt requirements.txt
ENV VIRTUAL_ENV=/home/siconos-user/siconosenv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN pip install --no-cache-dir -r requirements.txt

