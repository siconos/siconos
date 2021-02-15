FROM opensuse/leap:15.0

WORKDIR /home
COPY ci_gitlab/dockerfiles/requirements.txt /home
RUN zypper install -n -y \
    cmake \
    which \
    git-core \
    gcc-c++ \
    gcc-fortran \
    libboost_headers1_66_0-devel \
    gmp-devel \
    libboost_serialization1_66_0-devel \
    libboost_timer1_66_0-devel \
    openblas-devel \
    lapacke-devel \
    liblpsolve55-0 \
    python3-devel \
    python3-pip \
    swig \
    cppunit-devel \
    libbullet-devel \
    hdf5-devel \
    doxygen \
    suitesparse-devel \
    valgrind && zypper clean --all && pip3 install --upgrade pip && pip3 install -U -r /home/requirements.txt > /dev/null


