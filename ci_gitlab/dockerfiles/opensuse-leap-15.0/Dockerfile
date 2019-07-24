FROM opensuse/leap:15.0

RUN zypper install -n -y \
    cmake \
    which \
    git-core \
    gcc-c++ \
    gcc-fortran \
    libboost_headers1_66_0-devel \
    gmp-devel \
    libboost_serialization1_66_0-devel \
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
    valgrind
RUN zypper clean --all
WORKDIR /home
RUN pip3 install --upgrade pip
COPY requirements.txt /home
RUN pip3 install -U -r /home/requirements.txt > /dev/null


