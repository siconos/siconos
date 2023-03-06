FROM debian:bullseye
RUN apt update  && apt install -y -qq \
    	cmake \
        git \
	g++ \
	gfortran \
	vim \
	python3-dev \
	python3-pip \
	wget \
	libopenblas-dev \
	liblapacke-dev \
	libcppunit-dev \
        libhdf5-mpi-dev \
        libsuitesparse-dev \
        libgmp-dev \
	pkg-config \
	libbullet-dev \
	libcppunit-dev \
	valgrind
COPY ci_gitlab/dockerfiles/install_boost.sh  /home/
ENV CI_PROJECT_DIR /home
RUN sh /home/install_boost.sh
RUN mv /home/deps/boost_1_81_0/boost /usr/local/include
COPY ci_gitlab/dockerfiles/requirements.txt /home
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install -r /home/requirements.txt


