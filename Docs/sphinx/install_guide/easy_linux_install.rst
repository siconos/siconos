.. _easy_install:

Some tips to install siconos "easily" on linux-like systems
===========================================================

An easy way to get all the required libraries is to use one of the following command.
It will install the necessary packages for compiling all the components of Siconos.

Fedora
""""""

::
   yum install cmake gcc gcc-gfortran boost-devel atlas atlas-devel atlas-sse3 gmp gmp-devel cppunit cppunit-devel pcre-devel python-matplotlib scipy python-numeric python-devel swig python-py

Fedora 16, 17 or 18
"""""""""""""""""""

::
   yum install cmake gcc gcc-gfortran boost-devel atlas atlas-devel atlas-sse3 gmp gmp-devel cppunit cppunit-devel pcre-devel python-matplotlib scipy python-numeric python-devel swig pytest

Fedora 20
"""""""""

We recommand you to use openblas instead of atlas since it provides LAPACKE functions. Power users can still choose their favorite BLAS and LAPACK(E) vendor.::

  yum install cmake gcc gcc-gfortran boost-devel openblas openblas-devel gmp gmp-devel cppunit cppunit-devel python-matplotlib scipy python-numeric python-devel swig pytest

Ubuntu Precise
""""""""""""""

::
   apt-get install cmake gfortran gcc doxygen libboost-graph-dev libboost-filesystem-dev libboost-serialization-dev libatlas-base-dev libblas-dev libgmp3-dev libcppunit-dev libpcre3-dev python-matplotlib  python-numpy  python-scipy python-dev python-codespeak-lib make g++

Ubuntu Quantal
""""""""""""""
::
   apt-get install cmake gfortran gcc doxygen libboost-graph-dev libboost-filesystem-dev libboost-serialization-dev libatlas-base-dev libblas-dev libgmp3-dev libcppunit-dev libpcre3-dev python-matplotlib  python-numpy  python-scipy python-dev python-codespeak-lib make g++ python-pytest

Debian Squeeze
""""""""""""""
::
   apt-get install cmake gfortran gcc libboost-graph-dev libatlas-base-dev libgmp3-dev libcppunit-dev libpcre3-dev python-matplotlib  python-numpy  python-scipy python-dev python-codespeak-lib python-py swig2.0 doxygen libboost-filesystem-dev 

You will also need to use the backport repository to have a cmake recent enough::

  echo "deb http://backports.debian.org/debian-backports squeeze-backports main" >> /etc/apt/sources.list
  apt-get update
  apt-get -y --force-yes -t squeeze-backports install cmake

Mac Os
""""""

We recommend you to use brew to install the dependencies. Even if you have clang, you will need to install gcc because some fortran code has to be compiled in Numerics.
Here is a (probably) incomplete list of ports you have to install: wget, gcc, gmp, cppunit, doxygen, swig-python, py27-numpy, py27-py
