ARG MAMBAIMAGE=mambaorg/micromamba:1.5.8@sha256:475730daef12ff9c0733e70092aeeefdf4c373a584c952dac3f7bdb739601990
FROM $MAMBAIMAGE AS mamba-base
ENV DEBIAN_FRONTEND=noninteractive
USER root
# https://hub.docker.com/r/debian/snapshot?tag=stable-20231218
# To avoid to pin each package, we get debian snapshot, frozen at a given date
# Chosen date: creation of the mamba image
ARG DEBIAN_SNAPSHOT=20240820T174758Z
RUN echo "deb [check-valid-until=no] http://snapshot.debian.org/archive/debian/${DEBIAN_SNAPSHOT} bookworm main" > /etc/apt/sources.list \
	&& echo "deb [check-valid-until=no] http://snapshot.debian.org/archive/debian-security/${DEBIAN_SNAPSHOT} bookworm-security main" >> /etc/apt/sources.list

RUN apt update  && apt install --no-install-recommends -y -qq \
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
	doxygen \
	patch \
	graphviz \
	texlive-latex-base &&  apt autoclean -y && apt autoremove -y&& rm -rf /var/lib/apt/lists/*
RUN micromamba config prepend channels conda-forge

ARG MAMBA_YAML_FILE=ci_gitlab/siconos-freeze-2026-07.yml
ARG MAMBA_LOCKFILE=ci_gitlab/siconos-conda-lock.yml
FROM mamba-base AS mamba-occenv
USER $MAMBA_USER
WORKDIR /home/mambauser
RUN micromamba config prepend channels conda-forge && micromamba config set channel_priority strict
COPY --chown=$MAMBA_USER:$MAMBA_USER $MAMBA_LOCKFILE /home/conda-lock.yml
RUN micromamba install -n base -y -f /home/conda-lock.yml &&  micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python, cmake ... will not be found)
ENV SHELL /bin/bash
ENV PATH /opt/conda/bin:$PATH

FROM mamba-occenv AS siconos-doc
# Something breaks the swig/docstring/doxygen doc process in bullet. Fix it
COPY ci_gitlab/dockerfiles/legacy-full-2026-07/bullet.patch /home/mamba-user/
USER root
RUN patch --ignore-whitespace -uf /opt/conda/include/bullet/BulletCollision/Gimpact/btGImpactShape.h  -i /home/mamba-user/bullet.patch
USER $MAMBA_USER
