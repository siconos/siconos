FROM jupyter/scipy-notebook
USER root
ENV JUPYTER_ENABLE_LAB yes
RUN  apt update && apt install -y -qq \
        cmake \
        make \
        gfortran \
        g++ \
        libgfortran-10-dev \
        doxygen \
        libcppunit-dev && apt autoclean -y && apt autoremove -y && rm -rf /var/lib/apt/lists/*
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install pytest lxml pyhull vtk
RUN conda install -c conda-forge boost hdf5  swig bullet -y
RUN conda install -c conda-forge lpsolve55 suitesparse -y
WORKDIR /home
# COPY ci_gitlab/dockerfiles/install_bullet.sh .
ENV CI_PROJECT_DIR /home
# RUN sh /home/install_bullet.sh
