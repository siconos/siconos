# 1.ubuntu 26.04 image where we add all Siconos-required ubuntu packages, as root
#
FROM ubuntu:26.04 AS siconos-base
ENV TZ=Europe/Paris
ENV DEBIAN_FRONTEND=noninteractive
RUN  apt-get update && apt-get upgrade -y && apt-get install --no-install-recommends -y -qq \
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
	libboost-all-dev\
	nlohmann-json3-dev \
	libgl-dev \
	libgl1 \
	libeigen3-dev \
	gdb \
	&&  apt-get autoclean -y && apt-get autoremove -y&& rm -rf /var/lib/apt/lists/*
RUN useradd -ms /bin/bash siconos-user

# 2. Add a python venv for siconos-user
FROM siconos-base AS siconos-base-python
USER siconos-user
WORKDIR /home/siconos-user
COPY ci_gitlab/dockerfiles/requirements.txt .
ENV VIRTUAL_ENV=/home/siconos-user/siconosenv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN pip install --no-cache-dir -U -r requirements.txt
