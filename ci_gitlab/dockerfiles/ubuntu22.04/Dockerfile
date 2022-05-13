FROM ubuntu:22.04
ENV TZ=Europe/Paris
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && apt update  && apt install -y -qq \
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
	vim
WORKDIR /home
COPY ci_gitlab/dockerfiles/requirements4doc.txt /home
# vtk is not available in pypi for python3.10
RUN apt install -y -qq python3-vtk9
RUN pip3 install -U numpy scipy lxml pytest matplotlib h5py pyhull
RUN pip3 install -U -r /home/requirements4doc.txt
RUN apt autoclean -y && apt autoremove -y&& rm -rf /var/lib/apt/lists/*
# Something breaks the swig/docstring/doxygen doc process in bullet. Fix it
COPY ci_gitlab/dockerfiles/ubuntu22.04/bullet.patch /home
RUN patch -uf /usr/include/bullet/BulletCollision/Gimpact/btGImpactShape.h  -i /home/bullet.patch