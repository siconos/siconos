# Build an image including opencascade
# Req: opencascade-7.4.0.tgz available where docker is launched
FROM ubuntu:23.04
ENV TZ=Europe/Paris
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && apt update  && apt upgrade -y && apt install -y -qq \
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
	doxygen \
	swig \
	libxrender1 \
	graphviz \
	texlive-latex-base \
	valgrind \
	vim \
	bzip2 \
	wget \
	nlohmann-json3-dev \
	wget \
	libglu1-mesa-dev libgl1-mesa-dev libxmu-dev libxi-dev \
	libfreetype6-dev tk-dev rapidjson-dev libpcre2-dev
WORKDIR /home
ENV CI_PROJECT_DIR /home
COPY ci_gitlab/dockerfiles/requirements.txt /home
ENV VIRTUAL_ENV=/home/siconosenv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN pip install -U -r /home/requirements.txt
# RUN cd /home ;wget https://dev.opencascade.org/system/files/occt//OCC_7.4.0_release/opencascade-7.4.0.tgz
COPY opencascade-7.4.0.tgz /home
RUN cd /home; tar -zxvf opencascade-7.4.0.tgz
RUN cmake -S /home/opencascade-7.4.0 -B /home/cmake-build -DINSTALL_DIR=/home/install-opencascade -DBUILD_RELEASE_DISABLE_EXCEPTIONS=OFF
RUN cmake --build /home/cmake-build -j 8
RUN cmake --install /home/cmake-build
RUN apt install libpython3-dev -y
RUN apt autoclean -y && apt autoremove -y&& rm -rf /var/lib/apt/lists/*
RUN rm -rf /home/opencascade-7.4.0.tgz  /home/opencascade-7.4.0