{ pkgs ? import <nixpkgs> {} }:
with pkgs;

{
   siconos-numerics-only = callPackage ./siconos.nix {
   			 numerics_only=true;
			 enable_python=false;
			 blas_implem=pkgs.openblas.override { blas64 = false; };};
   siconos-numerics-python3 = callPackage ./siconos.nix {
   			 numerics_only=true;
			 enable_python=true;
			 blas_implem=pkgs.openblas.override { blas64 = false; };
			 };
   siconos-numerics-python2 = callPackage ./siconos.nix {
   			 numerics_only=true;
			 enable_python=true;
			 blas_implem=pkgs.openblas.override { blas64 = false; };
			 pythonX=python2;};
   siconos-full = callPackage ./siconos.nix {
			 enable_python=false;
			 blas_implem=pkgs.openblas.override { blas64 = false; };};
   siconos-full-python = callPackage ./siconos.nix {
			 enable_python=true;
			 blas_implem=pkgs.openblas.override { blas64 = false; };};
   siconos-test = callPackage ./sico_light.nix {
			 enable_python=true;
			 blas_implem=pkgs.openblas.override { blas64 = false; };};


   sico-py = pkgs.python35.withPackages (ps: with ps; [ ./siconos-numerics-python3 pip numpy ipython h5py matplotlib pytest lxml scipy]);

}
