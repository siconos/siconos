#!bin/bash
apt update -qq && apt install -y -qq cmake git-core wget make \
			      libboost-dev libgmp-dev swig gcc gfortran g++ liblapack-dev libatlas-base-dev \
			      lp-solve liblpsolve55-dev python3-dev libpython3-dev bash swig doxygen python3-dev python3-pip graphviz

pip3 install -r ./docs/requirements.txt
pip3 install git+https://github.com/sphinx-contrib/youtube.git
mkdir build
cd build
export LANG=C.UTF-8 # Required, else doxy2swig fails!
cmake ../ -DUSER_OPTIONS_FILE=$PWD/../CI/siconos_docs.cmake -DUSE_EXHALE=ON 
make doc
