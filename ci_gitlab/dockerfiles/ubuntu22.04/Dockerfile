FROM ubuntu:22.04 AS siconos-base
ENV TZ=Europe/Paris
ENV DEBIAN_FRONTEND=noninteractive
RUN  apt update && apt upgrade -y && apt install --no-install-recommends -y -qq \
        cmake \
        git-core \
        make \
	g++ \
        gfortran \
        libgmp-dev \
	libboost-dev \
        libopenblas-dev \
        liblapacke-dev \
        libcppunit-dev \
	libhdf5-dev \
	lp-solve \
        liblpsolve55-dev \
	libbullet-dev \
	libbullet-extras-dev \
	python3-full \
	libpython3-dev \
	swig \
	libxrender1 \
	valgrind \
	vim \
	bzip2 \
	wget \
	nlohmann-json3-dev \
	libgl-dev \
	&&  apt autoclean -y && apt autoremove -y&& rm -rf /var/lib/apt/lists/*
RUN useradd -ms /bin/bash siconos-user


FROM siconos-base AS siconos-base-python
USER siconos-user
WORKDIR /home/siconos-user
COPY ci_gitlab/dockerfiles/requirements.txt .
ENV VIRTUAL_ENV=/home/siconos-user/siconosenv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN pip install --no-cache-dir -U -r requirements.txt

FROM siconos-base-python AS siconos-base-doc
USER root
ENV DEBIAN_FRONTEND=noninteractive
RUN  apt update && apt install --no-install-recommends -y -qq \
	doxygen \
	patch \
	graphviz \
	texlive-latex-base &&  apt autoclean -y && apt autoremove -y&& rm -rf /var/lib/apt/lists/*
USER siconos-user
WORKDIR /home/siconos-user
COPY ci_gitlab/dockerfiles/requirements4doc.txt .
RUN pip install --no-cache-dir -U -r requirements4doc.txt
# Something breaks the swig/docstring/doxygen doc process in bullet. Fix it
COPY ci_gitlab/dockerfiles/ubuntu22.04/bullet.patch .
USER root
RUN patch --ignore-whitespace -uf /usr/include/bullet/BulletCollision/Gimpact/btGImpactShape.h  -i /home/siconos-user/bullet.patch
USER siconos-user