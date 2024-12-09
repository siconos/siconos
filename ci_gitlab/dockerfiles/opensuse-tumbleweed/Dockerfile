FROM opensuse/tumbleweed AS siconos-base

RUN zypper install -n -y \
    cmake \
    which \
    git-core \
    gcc-c++ \
    gcc-fortran \
    libboost_headers-devel \
    gmp-devel \
    libboost_serialization-devel \
    libboost_timer-devel \
    openblas-common-devel \
    lapacke-devel \
    liblpsolve55-0 \
    python-devel \
    cppunit-devel \
    libbullet-devel \
    hdf5-devel \
    doxygen \
    suitesparse-devel \
    valgrind && zypper clean --all
RUN useradd -ms /bin/bash siconos-user   

FROM siconos-base AS siconos-base-python
USER siconos-user
WORKDIR /home/siconos-user
COPY ci_gitlab/dockerfiles/requirements.txt .
ENV VIRTUAL_ENV=/home/siconos-user/siconosenv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN pip install --no-cache-dir -U -r requirements.txt


