{ pkgs ? import <nixpkgs> {},
  gcc ? pkgs.gcc ,
  gfortran ? pkgs.gfortran,
  numerics_only ? false,
  enable_python ? true,
  python ? pkgs.python36,
  blas_name ? "openblas",
  blas_implem ? pkgs.openblas.override { blas64 = false; },
  }:
  
with pkgs;

let  version = "1.0.0" ; in
 stdenv.mkDerivation rec {

 name = "siconos-${version}";
 buildInputs = [ cmake pkgconfig gcc gfortran python blas_implem pkgs.gmp pkgs.boost]
 	     ++ stdenv.lib.optional (numerics_only == false)[pkgs.boost];	     
 cmakeFlags = [ " " ]
    ++ stdenv.lib.optional (numerics_only == true) [ "-DCOMPONENTS=externals;numerics -DWITH_CXX=OFF" ]
    ++ stdenv.lib.optional (enable_python != true) [ "-DWITH_PYTHON_WRAPPER=OFF" ]
    ++ stdenv.lib.optional (enable_python != true) [ "-DWITH_BLAS=${blas_name}" ];
    
    

hardeningDisable = [ "format" ];
 src = ./.;
    }