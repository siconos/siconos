FROM centos:7
RUN yum update -y && yum install -y -qq \
     https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
RUN yum update -y && yum install -y -qq \   
        wget \
        git-all \
        boost169-devel \
        boost169-serialization \
        boost169-filesystem \
        boost169-timer \
        boost169-chrono \
        gmp-devel \
        gcc-c++ \
        gcc-gfortran \
        swig3 \
        lpsolve-devel \
        cppunit-devel \
        suitesparse-devel \
        python3-devel \
        python3-pip \
        openblas-devel && yum clean all
RUN wget https://github.com/Kitware/CMake/releases/download/v3.16.0/cmake-3.16.0-Linux-x86_64.sh && \
        printf 'y\nn\n' | sh cmake-3.16.0-Linux-x86_64.sh --prefix=/usr/local && \
        pip3 install numpy scipy h5py packaging pytest lxml vtk


