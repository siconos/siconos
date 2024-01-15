FROM mambaorg/micromamba:1.4.7
USER root
RUN apt update  && apt upgrade -y && apt install -y -qq \
        git-core \
        make \
	g++ \
        gfortran \
	libgl-dev
COPY --chown=$MAMBA_USER:$MAMBA_USER ci_gitlab/sicolabenv.yml /home/sicolabenv.yml
RUN micromamba config prepend channels conda-forge && micromamba self-update && micromamba config set channel_priority strict
RUN micromamba install -y -f /home/sicolabenv.yml && \
    micromamba clean --all --yes

WORKDIR /home
ENV CI_PROJECT_DIR /home
