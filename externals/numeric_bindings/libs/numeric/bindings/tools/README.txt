
Using the Python Generator
----
To test the python generator, you will need the lapack sources. Make sure the path to your 
lapack sources is also set in the lapack_parser.py file; look for the variable lapack_src_path.
In addition, check the same variable (lapack_src_path) in the blas_generator.py file.

On a debian-based system, acquiring the sources of BLAS and LAPACK may be achieved by 

$ cd tools
$ apt-get source blas lapack
$ apt-get install nvidia-cuda-dev

this will give you a subdirectory called 'blas-X.X' and 'lapack-X.X.X'. Running the generator may be done by

$ ./lapack_generator > output_lapack.txt
$ ./blas_generator > output_blas.txt

from within the tools directory. The generators are quite verbose, it is recommended to put their output
in a file.

To use the compile_test, you need cmake. Please do a 

$ mkdir build
$ cd build
$ cmake ..
$ make

to try to compile all test files. ]
