FROM ubuntu:23.04
ENV TZ=Europe/Paris
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && apt update  && apt upgrade -y && apt install -y -qq \
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
	python3-dev \
        python3-full \
        doxygen \
	libxrender1 \
	graphviz \
	texlive-latex-base \
	valgrind \
	vim \
	bzip2 \
	wget \
	nlohmann-json3-dev
WORKDIR /home
ENV CI_PROJECT_DIR /home
#COPY ci_gitlab/dockerfiles/install_boost.sh  /home/
#RUN sh /home/install_boost.sh
#RUN mv /home/deps/boost_1_81_0/boost /usr/local/include
COPY ci_gitlab/dockerfiles/requirements.txt /home
COPY ci_gitlab/dockerfiles/requirements4doc.txt /home
RUN python3 -m venv /home/pysicoenv
ENV PATH="/home/pysicoenv/bin:$PATH"
RUN python3 -m pip install -U pip
RUN python3 -m pip install -U -r /home/requirements.txt
RUN python3 -m pip install -U -r /home/requirements4doc.txt
RUN apt autoclean -y && apt autoremove -y&& rm -rf /var/lib/apt/lists/*
# Something breaks the swig/docstring/doxygen doc process in bullet. Fix it
#COPY ci_gitlab/dockerfiles/ubuntu22.04/bullet.patch /home
#RUN patch --ignore-whitespace -uf /usr/include/bullet/BulletCollision/Gimpact/btGImpactShape.h  -i /home/bullet.patch
