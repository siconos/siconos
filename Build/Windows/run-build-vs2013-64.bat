REM @echo off

set PATH=C:\Program Files\mingw-w64\x86_64-4.8.3-posix-seh-rt_v3-rev0\mingw64\bin;%PATH%

rm -rf E:/build-vs2013-amd64
rm -rf E:/install-vs2013-amd64

mkdir E:\build-vs2013-amd64
cd /D E:\build-vs2013-amd64

call %~dp0\build-siconos-vs2013-64.bat %1\Build

cd /D E:\
call zip -r siconos-vs2013-amd64.zip install-vs2013-amd64
copy siconos-vs2013-amd64.zip %1/

REM rm -rf E:/Users/ci/install-vs2013-amd64
REM rm -rf E:/Users/ci/build-vs2013-amd64
