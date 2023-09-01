# FROM jupyter/scipy-notebook:latest
# USER root
# ENV JUPYTER_ENABLE_LAB yes
# RUN  apt update && apt install -y -qq \
#         doxygen && apt autoclean -y && apt autoremove -y && rm -rf /var/lib/apt/lists/*
# RUN conda install -c conda-forge cmake  cxx-compiler gfortran  boost hdf5 cppunit swig bullet -y
# RUN conda install -c conda-forge liblapacke -y
# RUN conda install -c conda-forge lpsolve55 suitesparse -y
# RUN conda install -c conda-forge pyhull wheel numpy scipy lxml pytest matplotlib h5py pybind11 vtk packaging
# WORKDIR /home
# ENV CI_PROJECT_DIR /home

FROM mambaorg/micromamba:1.4.7
USER root
RUN apt update  && apt upgrade -y && apt install -y -qq \
        git-core \
        make \
	g++ \
        gfortran
COPY --chown=$MAMBA_USER:$MAMBA_USER ci_gitlab/sicolabenv.yml /home/sicolabenv.yml
RUN micromamba config prepend channels conda-forge && micromamba self-update && micromamba config set channel_priority strict
RUN micromamba install -y -f /home/sicolabenv.yml && \
    micromamba clean --all --yes

WORKDIR /home
ENV CI_PROJECT_DIR /home
