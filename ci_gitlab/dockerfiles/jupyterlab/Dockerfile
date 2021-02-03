FROM jupyter/scipy-notebook
USER root
ENV JUPYTER_ENABLE_LAB yes
RUN  apt update && apt install -y -qq \
        cmake \
        make \
        gfortran \
        libgfortran-10-dev \
        libboost-dev \
        libboost-filesystem-dev \
        libboost-timer-dev \
        libboost-chrono-dev \
        libgmp-dev \
        libopenblas-dev \
        liblapacke-dev \
        lp-solve \
        liblpsolve55-dev \
        libhdf5-dev \
        libboost-serialization-dev \
        libbullet-dev \
        libfreetype6-dev \
        freeglut3-dev \
        swig \
        doxygen \
        libcppunit-dev \
        libbullet-dev \
        libfreetype6-dev \
        freeglut3-dev && apt autoclean -y && apt autoremove -y && rm -rf /var/lib/apt/lists/*
RUN pip3 install pytest
