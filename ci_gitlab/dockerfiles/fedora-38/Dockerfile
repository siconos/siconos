FROM fedora:38
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
	&& dnf clean all
WORKDIR /home
COPY ci_gitlab/dockerfiles/requirements.txt /home
RUN python3 -m pip install -U -r /home/requirements.txt
ENV CI_PROJECT_DIR /home

