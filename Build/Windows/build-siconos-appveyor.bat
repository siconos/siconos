@echo off

REM IMPORTANT VARIABLES :: you have to set them properly
set VCVARSALL_BAT="C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat"
set VC_ARCH=amd64

call %VCVARSALL_BAT% %VC_ARCH%

set PATH=C:\externals-siconos\bin;%PATH%
set LIB=C:\externals-siconos\lib;%LIB%
set INCLUDE=C:\externals-siconos\include;%INCLUDE%

REM We might have to switch boost version --xhub
set CMAKE_INSTALL_PREFIX=C:/install-vs2015-amd64
set REL_TYPE="Release"


set PATH=%CMAKE_INSTALL_PREFIX%;%PATH%

REM does not work if you don't have msys64 at this precise place
for /f "delims=" %%a in ('C:\msys64\usr\bin\bash -lc "cygpath -m `swig -swiglib`"') do @set SWIG_DIR=%%a

call cmake %* ^
-G"NMake Makefiles" ^
-DCMAKE_BUILD_TYPE=%REL_TYPE% ^
-DMODE=Continuous ^
-DCMAKE_INSTALL_PREFIX="%CMAKE_INSTALL_PREFIX%" ^
-DINSTALL_COMMAND="nmake;/I;install" ^
-DSWIG_DIR="%SWIG_DIR%" ^
-DBUILD_AS_CPP=1 ^
-DWITH_TESTING=1 ^
-DWITH_PYTHON_WRAPPER=0 ^
-DCMAKE_Fortran_COMPILER=gfortran ^
-DCMAKE_Fortran_COMPILER_FORCED=1;

call nmake.exe /I
call nmake.exe /I install
