FROM intel/oneapi-hpckit:latest
RUN apt update  && apt upgrade -y &&  apt install -y -qq \
    	# build-essential \
        # libgmp-dev \
        # libcppunit-dev \
	# libbullet-dev \
	# libbullet-extras-dev \
	# vim \
	# bzip2 \
	# libsuitesparse-dev \
	libgl1 \
	intel-oneapi-python
RUN apt autoclean -y && apt autoremove -y&& rm -rf /var/lib/apt/lists/*
ENV PATH /opt/intel/oneapi/intelpython/latest/bin/:$PATH
ENV CI_PROJECT_DIR /home
COPY ci_gitlab/dockerfiles/sicoenv.yml /home
COPY ci_gitlab/dockerfiles/siconosintelenv.yml  /home/
RUN conda config --prepend channels conda-forge
RUN conda env create -f /home/siconosintelenv.yml
#RUN conda config set channel_priority strict
#RUN conda install -y lxml pytest matplotlib h5py pybind11 vtk 
#COPY ci_gitlab/dockerfiles/install_boost.sh  /home/
#RUN mkdir -p $CI_PROJECT_DIR/build
#RUN mkdir -p $CI_PROJECT_DIR/deps
#RUN cd $CI_PROJECT_DIR/deps; git clone https://github.com/HDFGroup/hdf5.git
##RUN cd $CI_PROJECT_DIR/build; cmake $CI_PROJECT_DIR/deps/hdf5 -DCMAKE_C_COMPILER="mpiicc;-cc=icx" -DCMAKE_CXX_COMPILER="mpiicpc;-cxx=icpx" -DCMAKE_INSTALL_PREFIX=$CI_PROJECT_DIR/install/hdf5 -DHDF5_ENABLE_PARALLEL=ON;make -j 4 ; make install
#ENV HDF5_ROOT /home/install/hdf5
#RUN sh /home/install_boost.sh
#RUN mv /home/deps/boost_1_81_0/boost /usr/local/include
ENV CXX icpx
ENV CC icx
ENV FC ifx
ENV Python_EXECUTABLE /opt/intel/oneapi/intelpython/latest/bin//python3