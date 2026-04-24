# 1.ubuntu 24.04 image where we add all Siconos-required ubuntu packages, as root
# 
FROM ubuntu:24.04 AS siconos-base
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
	libeigen3-doc \
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


# 3. Extend base image with the packages required to build the documentation
FROM siconos-base-python AS siconos-base-doc
USER root
RUN  apt-get update && apt-get install --no-install-recommends -y -qq \
	doxygen \
	patch \
	graphviz \
	texlive-latex-base &&  apt-get autoclean -y && apt-get autoremove -y&& rm -rf /var/lib/apt/lists/*
USER siconos-user
WORKDIR /home/siconos-user
COPY ci_gitlab/dockerfiles/requirements4doc.txt .
RUN pip install --no-cache-dir -U -r requirements4doc.txt






FROM siconos-base-python AS siconos-base-parallel
USER root
RUN  apt-get update && apt-get install --no-install-recommends -y -qq \
    libomp-dev \
	libopenmpi-dev \
	curl \
	&&  apt-get autoclean -y && apt-get autoremove -y&& rm -rf /var/lib/apt/lists/*


FROM siconos-base-parallel AS siconos-petsc
USER siconos-user
# Install PETSc manually (to enable 64-bit-indices and remove warnings)
COPY ci_gitlab/dockerfiles/install_petsc.sh /home/siconos-user
ENV PETSC_INSTALL_DIR="$VIRTUAL_ENV"
RUN sh /home/siconos-user/install_petsc.sh
# setup to help siconos cmake process to find petsc 
ENV petsc_ROOT="$VIRTUAL_ENV" 



