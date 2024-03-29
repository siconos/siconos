{ stdenv, fetchurl, cmake, pkgconfig,
  pkgs ? import <nixpkgs> {},
  gcc ? pkgs.gcc ,
  gfortran ? pkgs.gfortran,
  numerics_only ? false,
  enable_python ? true,
  pythonX ? pkgs.python36,
  blas_name ? "openblas",
  blas_implem ? pkgs.openblas.override { blas64 = false; },
  boost ? pkgs.boost,
  gmp ? pkgs.gmp,
  version ? "4.1.0",
  enable_cxx ? true,
  enable_openmp ? false,
  }:
  
with pkgs;
with stdenv.lib;

let
	pythonenv = pythonX.withPackages (ps: with ps; [pip numpy ipython h5py matplotlib lxml scipy pytest]);
        boost-dev-meta = (pkgs.boost.meta // {outputsToInstall =["out" "dev"]; });
        boost-dev = pkgs.boost // {meta = boost-dev-meta;};               
in

stdenv.mkDerivation rec {
 name = "siconos-${version}";
 
 enableParallelBuilding = true;	

 nativeBuildInputs = [
    pkgconfig
    pythonX
   ]  
  ++ optional enable_python [pythonX.pkgs.wrapPython];
 
 buildInputs = [
    gcc
    blas_implem
    gmp
    boost
    cppunit
    boost-dev
    ]
  ++ optional (numerics_only != true) [ boost ];
  
  propagatedNativeBuildInputs = with pythonX.pkgs;
    if enable_python then [ cmake swig gfortran pythonenv]
    else [cmake gfortran];

 cmakeFlags = [ "-DBLA_VENDOR=${blas_name}" ]
    ++ optional (numerics_only == true) [ "-DCOMPONENTS=externals;numerics -DWITH_CXX=OFF" ]
    ++ optional (enable_python != true) [ "-DWITH_PYTHON_WRAPPER=OFF" ]
    ++ optional (enable_openmp == true) [ "-DWITH_OPENMP=ON" ];

    
    
 hardeningDisable = [ "format" ];
 src = ./.;



 
  postFixup = ''
    echo "Create links in bin ..."
    if test -e $out/nix-support/propagated-native-build-inputs; then
        ln -s $out/nix-support/propagated-native-build-inputs $out/nix-support/propagated-user-env-packages
    fi
   '';

 
  meta = with stdenv.lib; {
    homepage = http://siconos.gforge.inria.fr/;
    description = "Nonsmooth dynamical systems simulator";
    license = licenses.asl20;
    platforms = platforms.all;
  };


}
