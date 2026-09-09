# Build an image including opencascade
# Req: opencascade-7.4.0.tgz available where docker is launched
#### Mamba setup ###
ARG MAMBAIMAGE=mambaorg/micromamba:1.5.8
FROM $MAMBAIMAGE AS mamba-base
ENV DEBIAN_FRONTEND=noninteractive
USER root
RUN apt update  && apt upgrade -y && apt install --no-install-recommends -y -qq \
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
	&& apt autoclean -y && apt autoremove -y&& rm -rf /var/lib/apt/lists/*
RUN micromamba config prepend channels conda-forge && micromamba self-update

# Warning: works only of opencascade tgz is available in the Docker context dir
ARG MAMBA_YAML_FILE=ci_gitlab/siconos-opencascade.yml
FROM mamba-base AS mamba-occenv
USER $MAMBA_USER
WORKDIR /home/mambauser
RUN micromamba config prepend channels conda-forge && micromamba config set channel_priority strict
COPY --chown=$MAMBA_USER:$MAMBA_USER $MAMBA_YAML_FILE /home/siconos-env.yml
RUN micromamba install -n base -y -f /home/siconos-env.yml&& \
	micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python, cmake ... will not be found)


FROM mamba-occenv AS jupylab
ENV JUPYTER_PORT=8888
ENV JUPYTER_ENABLE_LAB yes
EXPOSE $JUPYTER_PORT
ENV SHELL /bin/bash
ENTRYPOINT ["/opt/conda/bin/jupyter", "lab","--ip=0.0.0.0","--allow-root", "--no-browser"]
