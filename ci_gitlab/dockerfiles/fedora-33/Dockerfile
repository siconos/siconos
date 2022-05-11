FROM fedora:33
RUN dnf --setopt=deltarpm=false install -y \
	make \
        cmake \
        python3 \
        python3-devel \
        python3-lxml \
        python3-pytest \
        python3-scipy \
        python3-matplotlib \
        python3-packaging \
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
	python3-h5py  && dnf clean all
WORKDIR /home
RUN python3 -m pip install vtk
COPY ci_gitlab/dockerfiles/install_bullet.sh .
ENV CI_PROJECT_DIR /home
RUN sh /home/install_bullet.sh

