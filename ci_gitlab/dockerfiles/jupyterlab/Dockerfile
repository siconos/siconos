FROM jupyter/scipy-notebook
USER root
ENV JUPYTER_ENABLE_LAB yes
RUN  apt update && apt install -y -qq \
        cmake \
        make \
        gfortran \
        g++ \
        libgfortran-10-dev \
        libboost-dev \
        libboost-filesystem-dev \
        libboost-timer-dev \
        libboost-chrono-dev \
        libgmp-dev \
	libsuitesparse \
        lp-solve \
        liblpsolve55-dev \
        libhdf5-dev \
        libboost-serialization-dev \
        libfreetype6-dev \
        freeglut3-dev \
        swig \
        doxygen \
        libcppunit-dev \
        libbullet-dev \
        libfreetype6-dev \
        freeglut3-dev && apt autoclean -y && apt autoremove -y && rm -rf /var/lib/apt/lists/*
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install pytest lxml pyhull
WORKDIR /home
COPY ci_gitlab/dockerfiles/install_bullet.sh .
ENV CI_PROJECT_DIR /home
RUN sh /home/install_bullet.sh

