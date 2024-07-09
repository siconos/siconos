##  - Prepare an image including Jupyter lab, based on mamba.
# Note: the conf to start properly the image as a notebook is done
# in the last step, after build, install, to avoid all the mess
# with entrypoint and CI.
#### Mamba setup ### 
ARG MAMBAIMAGE=mambaorg/micromamba:1.5.8
FROM $MAMBAIMAGE AS mamba-base
LABEL maintainer="Siconos team <siconos-team@inria.fr>"
USER root
RUN apt update  && apt upgrade -y && apt install --no-install-recommends -y -qq \
        git-core \
        make \
	g++ \
        gfortran \
	libgl-dev  &&  apt autoclean -y && apt autoremove -y&& rm -rf /var/lib/apt/lists/*
RUN micromamba config prepend channels conda-forge && micromamba self-update && \
    micromamba clean --all --yes


FROM mamba-base AS mambaenv
USER $MAMBA_USER
RUN micromamba config prepend channels conda-forge && micromamba config set channel_priority strict
COPY --chown=$MAMBA_USER:$MAMBA_USER ci_gitlab/sicolabenv.yml /home/sicolabenv.yml
RUN micromamba install -y -f /home/sicolabenv.yml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python, cmake ... will not be found)
ENV SHELL /bin/bash
ENV PATH /opt/conda/bin:$PATH
WORKDIR /home/mambauser
