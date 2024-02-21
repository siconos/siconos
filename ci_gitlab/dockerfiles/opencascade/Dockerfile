# Build an image including opencascade
# Req: opencascade-7.4.0.tgz available where docker is launched
FROM mambaorg/micromamba:1.5.6 AS mamba-base
ENV DEBIAN_FRONTEND=noninteractive
USER root
RUN apt update  && apt upgrade -y && apt install --no-install-recommends -y -qq \
        cmake \
        git-core \
        make \
	g++ \
        gfortran \
        libgmp-dev \
	libxrender1 \
	vim \
	bzip2 \
	wget \
	libboost1.81-all-dev \
	libglu1-mesa-dev libgl1-mesa-dev libxmu-dev libxi-dev \
	libfreetype6-dev tk-dev rapidjson-dev libpcre2-dev \
	&& apt autoclean -y && apt autoremove -y&& rm -rf /var/lib/apt/lists/*
RUN micromamba config prepend channels conda-forge && micromamba self-update
#RUN wget https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.bz2 && \
#    tar -jxvf boost_1_81_0.tar.bz2 && \
#    mv boost_1_81_0/boost /usr/local/include


# Warning: works only of opencascade tgz is available in the Docker context dir
FROM mamba-base AS mamba-occenv
USER $MAMBA_USER
WORKDIR /home/mambauser
#RUN micromamba config prepend channels conda-forge && micromamba config set channel_priority strict
COPY --chown=$MAMBA_USER:$MAMBA_USER ci_gitlab/siconos-opencascade.yml /home/siconos-opencascade.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER ci_gitlab/dockerfiles/opencascade/condarc /home/mambauser/.condarc
RUN micromamba install -v -y -f /home/siconos-opencascade.yml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python, cmake ... will not be found)
ENV JUPYTER_PORT=8888
#ENV JUPYTER_ENABLE_LAB yes
#EXPOSE $JUPYTER_PORT
ENV SHELL /bin/bash
#ENTRYPOINT ["/opt/conda/bin/jupyter", "lab","--ip=0.0.0.0","--allow-root", "--no-browser"]



# USER siconos-user
# WORKDIR /home/siconos-user
# COPY ci_gitlab/dockerfiles/requirements.txt .
# ENV VIRTUAL_ENV=/home/siconos-user/siconosenv
# RUN python3 -m venv $VIRTUAL_ENV
# ENV PATH="$VIRTUAL_ENV/bin:$PATH"
# RUN pip install --no-cache-dir -U -r requirements.txt
# COPY opencascade-7.4.0.tgz .
# RUN tar -zxvf opencascade-7.4.0.tgz
# RUN cmake -S /home/opencascade-7.4.0 -B /home/cmake-build -DINSTALL_DIR=/home/install-opencascade -DBUILD_RELEASE_DISABLE_EXCEPTIONS=OFF
# RUN cmake --build /home/cmake-build -j 8
# RUN cmake --install /home/cmake-build
# RUN rm -rf /home/opencascade-7.4.0.tgz  /home/opencascade-7.4.0

# FROM siconos-base-python AS siconos-base-pyocc
