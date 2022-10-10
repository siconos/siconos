FROM ubuntu:22.04
ENV TZ=Europe/Paris
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && apt update  && apt upgrade -y && apt install -y -qq \
        cmake \
        git-core \
        make \
	g++ \
        gfortran \
        libgmp-dev \
	libboost-dev \
        libopenblas-dev \
        liblapacke-dev \
        libcppunit-dev \
	libhdf5-dev \
	lp-solve \
        liblpsolve55-dev \
	libbullet-dev \
	libbullet-extras-dev \
	python3 \
        python3-pip \
	doxygen \
	swig \
	libxrender1 \
	graphviz \
	texlive-latex-base \
	valgrind \
	vim
WORKDIR /home
COPY ci_gitlab/dockerfiles/requirements4doc.txt /home
RUN pip3 install -U wheel numpy scipy lxml pytest matplotlib h5py pyhull vtk
RUN pip3 install -U -r /home/requirements4doc.txt
RUN apt autoclean -y && apt autoremove -y&& rm -rf /var/lib/apt/lists/*
# Something breaks the swig/docstring/doxygen doc process in bullet. Fix it
COPY ci_gitlab/dockerfiles/ubuntu22.04/bullet.patch /home
RUN patch --ignore-whitespace -uf /usr/include/bullet/BulletCollision/Gimpact/btGImpactShape.h  -i /home/bullet.patch
