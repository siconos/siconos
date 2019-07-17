FROM gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos/ubuntu18.04-oce
RUN apt update  && apt install -y -qq \
        graphviz
WORKDIR /home
COPY requirements.txt .
RUN pip3 install -U -r /home/requirements.txt
ENV LANG C.UTF-8 # Required, else doxy2swig fails!
RUN pip3 install git+https://github.com/sphinx-contrib/youtube.git
RUN apt clean && apt autoremove
