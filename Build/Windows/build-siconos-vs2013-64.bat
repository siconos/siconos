@echo off

REM IMPORTANT VARIABLES :: you have to set them properly
set VCVARSALL_BAT="C:\Program Files (x86)\VC\vcvarsall.bat"
set VC_ARCH=x86_amd64

call %VCVARSALL_BAT% %VC_ARCH%

REM set PATH=%PATH%;E:\OpenBLAS-64\lib
set PATH=E:\siconos-external-lib\lib;%PATH%
set PATH=E:\cppunit-1.13.2\lib;%PATH%
set PATH=E:\all_libs;%PATH%

set BOOST_INCLUDEDIR=E:/siconos-external-lib/include
set CMAKE_INSTALL_PREFIX=E:/install-vs2013-amd64
set REL_TYPE="Release"


set PATH=%CMAKE_INSTALL_PREFIX%;%PATH%

call cmake %* -DON_DASHBOARD=1 ^
-G"NMake Makefiles" ^
-DCMAKE_BUILD_TYPE=Release ^
-DMODE=Continuous ^
-DWITH_IO=0 ^
-DWITH_MECHANICS=1 ^
-DWITH_CONTROL=1 ^
-DWITH_EXAMPLES=1 ^
-DCMAKE_INSTALL_PREFIX="%CMAKE_INSTALL_PREFIX%" ^
-DINSTALL_COMMAND="nmake;/I;install" ^
-DFROM_REPO=0 ^
-DCTEST_OPTIONS="-j2;-V" ^
-DBUILD_TYPE=%REL_TYPE% ^
-Dproject_CMAKE_ARGS="-DSWIG_EXECUTABLE=C__22__/Chocolatey/bin/swig.bat;-DWITH_LAPACK=lapacke;-DCMAKE_Fortran_COMPILER=gfortran;-DCMAKE_Fortran_COMPILER_FORCED=1;__11__-DCMAKE_CXX_FLAGS=/DBOOST_ALL_NO_LIB__00__/EHsc__11__"

call nmake.exe /I
