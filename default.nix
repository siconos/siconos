{ pkgs ? import <nixpkgs> {} }:
with pkgs;

{
   siconos-numerics-only = callPackage ./siconos.nix {
   			 numerics_only=true;
			 enable_python=false;
			 blas_implem=pkgs.openblas.override { blas64 = false; };};
   siconos-full = callPackage ./siconos.nix {
			 enable_python=false;
			 blas_implem=pkgs.openblas.override { blas64 = false; };};
}
