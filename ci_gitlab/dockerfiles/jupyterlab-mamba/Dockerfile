##  - Prepare an image including (and starting) Jupyter lab, based on mamba.
#### Mamba setup ### 
ARG MAMBAIMAGE=mambaorg/micromamba:1.5.5
FROM $MAMBAIMAGE AS mamba-base
LABEL maintainer="Siconos team <siconos-team@inria.fr>"
USER root
RUN apt update  && apt upgrade -y && apt install --no-install-recommends -y -qq \
        git-core \
        make \
	g++ \
        gfortran \
	libgl-dev  &&  apt autoclean -y && apt autoremove -y&& rm -rf /var/lib/apt/lists/*
RUN micromamba config prepend channels conda-forge && micromamba self-update


FROM mamba-base AS mambalab
USER $MAMBA_USER
RUN micromamba config prepend channels conda-forge && micromamba config set channel_priority strict
COPY --chown=$MAMBA_USER:$MAMBA_USER ci_gitlab/sicolabenv.yml /home/sicolabenv.yml
RUN micromamba install -y -f /home/sicolabenv.yml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python, cmake ... will not be found)
ENV JUPYTER_PORT=8888
ENV JUPYTER_ENABLE_LAB yes
EXPOSE $JUPYTER_PORT
SHELL ["/bin/bash", "-l" ,"-c"]
ENV SHELL /bin/bash
ENV PATH /opt/conda/bin:$PATH
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/opt/conda/bin/jupyter", "lab","--ip=0.0.0.0","--allow-root", "--no-browser"]
WORKDIR /home/mambauser