FROM gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos/sources/ubuntu18.04
RUN apt update  && apt upgrade -y && apt install -y -qq \
    libopenmpi-dev \
    libparmetis-dev \
    libptscotch-dev \
    openssh-client \
    libmumps-dev && apt autoclean -y&& apt autoremove -y&& rm -rf /var/lib/apt/lists/*
RUN python3 -m pip install mpi4py
