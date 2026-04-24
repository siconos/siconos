##  - Prepare an image including Jupyter lab, based on mamba.
# Note: the conf to start properly the image as a notebook is done
# in the last step, after build, install, to avoid all the mess
# with entrypoint and CI.

# 1. Mamba-based image where we add all Siconos-required ubuntu packages, as root
# 
ARG MAMBAIMAGE=mambaorg/micromamba:2.5
FROM $MAMBAIMAGE AS mamba-base
ENV DEBIAN_FRONTEND=noninteractive
LABEL maintainer="Siconos team <siconos-team@inria.fr>"
USER root
RUN apt-get update  && apt-get upgrade -y && apt-get install --no-install-recommends -y -qq \
        git-core \
        make \
	    g++ \
        gfortran \
    	libgl-dev \
		libxrender1 \
		bzip2 \
		vim \
		wget \
		libglu1-mesa-dev libgl1-mesa-dev libxmu-dev libxi-dev \
		libfreetype6-dev tk-dev rapidjson-dev libpcre2-dev \
        &&  apt-get autoclean -y && apt-get autoremove -y&& rm -rf /var/lib/apt/lists/*
RUN micromamba config prepend channels conda-forge && \
    micromamba clean --all --yes

# 2. Prepare a micromamba standard env, ready for Siconos build, for mamba user
# Based on sicolabenv.yaml.
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


# 3. Prepare another micromamba env including opencascade and pyocc, 
# ready for Siconos build with mechanisms component.
FROM mamba-base AS mamba-occenv
USER $MAMBA_USER
WORKDIR /home/mambauser
RUN micromamba config prepend channels conda-forge && micromamba config set channel_priority strict
COPY --chown=$MAMBA_USER:$MAMBA_USER ci_gitlab/siconos-opencascade.yml /home/siconos-env.yml
RUN micromamba install -n base -y -f /home/siconos-env.yml&& \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python, cmake ... will not be found)

# 4. jupyter setup for micromamba+occ
FROM mamba-occenv AS jupylab
ENV JUPYTER_PORT=8888
ENV JUPYTER_ENABLE_LAB=yes
EXPOSE $JUPYTER_PORT
ENV SHELL=/bin/bash
ENTRYPOINT ["/opt/conda/bin/jupyter", "lab","--ip=0.0.0.0", "--no-browser"]
CMD ["--port=8888"]
