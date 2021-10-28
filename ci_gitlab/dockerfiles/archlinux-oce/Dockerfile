FROM archlinux/base
WORKDIR /home
ENV CI_PROJECT_DIR /home
COPY ci_gitlab/dockerfiles/install_oce.sh /home/install_oce.sh
RUN pacman -Sy && pacman --noconfirm -S \
        git \
        boost \
        boost-libs \
        gmp \
        cmake \
        make \
        gcc \
        gcc-fortran \
        swig \
        openblas \
        lapacke \
        cblas \
        suitesparse \
        python3 \
        python-pytest \
        python-scipy \
        python-h5py \
        python-pip \
        python-lxml \
        python-packaging \
        cppunit \
        bullet \
        freetype2 \
        freeglut \
        tk tcl \
        glu && pacman --noconfirm -Scc
RUN sh /home/install_oce.sh clone_oce