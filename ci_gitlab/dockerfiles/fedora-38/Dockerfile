FROM fedora:38 AS siconos-base
RUN dnf --setopt=deltarpm=false install -y \
	make \
        cmake \
        python3 \
        python3-devel \
        lpsolve-devel \
        gcc-gfortran \
        gcc-c++ \
        boost-devel \
        boost-filesystem \
        boost-serialization \
        gmp-devel \
        cppunit-devel \
        cppunit \
        hdf5-devel \
        openblas-devel \
        lapack-devel \
        suitesparse-devel \
        swig \
	libglvnd-glx \
        git-core \
        bullet-devel \
	libXrender-devel \
	eigen3-devel \
	&& dnf clean all
RUN useradd -ms /bin/bash siconos-user

FROM siconos-base AS siconos-base-python
USER siconos-user
WORKDIR /home/siconos-user
COPY ci_gitlab/dockerfiles/requirements.txt .
ENV VIRTUAL_ENV=/home/siconos-user/siconosenv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN pip install --no-cache-dir -U -r requirements.txt

