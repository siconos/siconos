FROM ubuntu:22.04 AS siconos-base
ENV TZ=Europe/Paris
ENV DEBIAN_FRONTEND=noninteractive
RUN  apt update && apt upgrade -y && apt install --no-install-recommends -y -qq \
    cmake \
    git-core \
    make \
	g++ \
    gfortran \
    libgmp-dev \
    libopenblas-dev \
    liblapacke-dev \
    libcppunit-dev \
	libhdf5-dev \
	lp-solve \
    liblpsolve55-dev \
	libbullet-dev \
	libbullet-extras-dev \
	python3-full \
	libpython3-dev \
	libxrender1 \
	valgrind \
	vim \
	bzip2 \
	wget \
	nlohmann-json3-dev \
	libgl-dev \
	libgl1 \
	libeigen3-dev \
	libeigen3-doc \
	&&  apt autoclean -y && apt autoremove -y&& rm -rf /var/lib/apt/lists/*
RUN useradd -ms /bin/bash siconos-user
RUN  wget https://archives.boost.io/release/1.81.0/source/boost_1_81_0.tar.bz2 && \
    tar -jxvf boost_1_81_0.tar.bz2 && \
    cd boost_1_81_0 && ./bootstrap.sh && ./b2 install

FROM siconos-base AS siconos-base-python
USER siconos-user
WORKDIR /home/siconos-user
COPY ci_gitlab/dockerfiles/requirements.txt .
ENV VIRTUAL_ENV=/home/siconos-user/siconosenv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN pip install --no-cache-dir -U -r requirements.txt
