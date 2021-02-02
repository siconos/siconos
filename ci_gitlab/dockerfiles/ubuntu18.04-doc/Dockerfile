# REGISTRY_PATH will be updated (sed) during kaniko call
FROM gricad-registry.univ-grenoble-alpes.fr/REGISTRY_PATH/sources/ubuntu18.04
RUN apt update  && apt install -y -qq \
        graphviz && apt autoclean -y && apt autoremove -y && rm -rf /var/lib/apt/lists/*
WORKDIR /home
COPY requirements.txt .
RUN pip3 install -U -r /home/requirements.txt
ENV LANG C.UTF-8 # Required, else doxy2swig fails!
RUN pip3 install git+https://github.com/sphinx-contrib/youtube.git
