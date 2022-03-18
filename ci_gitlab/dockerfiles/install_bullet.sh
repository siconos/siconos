#!bin/bash

# Get number of procs
if  [ -x "$(command -v nproc)" ]; then
   export nbprocs=`nproc --all`  # linux
elif  [ -x "$(command -v sysctl)" ]; then
   export nbprocs=`sysctl -n hw.ncpu` # macos
else
   export nbprocs=2
fi

# Check if CI_PROJECT_DIR is set AND not empty
: ${CI_PROJECT_DIR:?"Please set environment variable CI_PROJECT_DIR with 'siconos' repository (absolute) path."}

mkdir -p $CI_PROJECT_DIR/build/bullet3
cd $CI_PROJECT_DIR/
git clone https://github.com/bulletphysics/bullet3.git > /dev/null
git checkout tags/3.21
cd $CI_PROJECT_DIR/build/bullet3
cmake -DBUILD_PYBULLET=ON -DBUILD_PYBULLET_NUMPY=ON -DCMAKE_BUILD_TYPE=Release -DOpenGL_GL_PREFERENCE=GLVND $CI_PROJECT_DIR/bullet3 -Wno-dev
make -j $nbprocs > /dev/null
echo "----> install bullet ..."
make install > /dev/null
# Clean up
rm -rf $CI_PROJECT_DIR/bullet3
rm -rf $CI_PROJECT_DIR/build/bullet
