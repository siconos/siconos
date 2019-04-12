FROM opensuse/leap:15.0
RUN zypper install -y \
        cmake \
        which \
        git-core \
        gcc-c++ \
        gcc-fortran \
        boost-devel \
        gmp-devel \
        libboost_serialization-devel \
        openblas-dev \
        lapacke-devel \
        liblpsolve55-0 \
        python3-dev \
        python3-pip \
        cppunit-devel \
        libbullet-dev \
        hdf5-devel \
        doxygen \
        valgrind
RUN zypper clean --all
WORKDIR /home
RUN pip3 install --upgrade pip
COPY requirements.txt /home
RUN pip3 install -U -r /home/requirements.txt > /dev/null


