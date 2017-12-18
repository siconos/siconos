Some words about nix and siconos
================================

Using nix to build/install siconos

Prereq: source /applis/site/nix.sh

* List which derivations (i.e. which siconos config) are available from current sources

IN SOURCE_DIR:

nix-env -qaP -f .

siconos-full           siconos-1.0.0
siconos-numerics-only  siconos-1.0.0

* build siconos, install it in the store

nix-build -A attribute-name

==> build attribute-name (attribute-name being for example siconos-numerics-only)

Creates a link (@result) in the current dir to siconos in nixstore.

* build and install siconos in your nix environment:

nix-env -f . -iA fclib-gcc6

check for example ~/.nix-profile/lib/ for libsiconos...so



About nix :

* https://nixos.org/nixos/manual/
* http://nix-cookbook.readthedocs.io/en/latest/nix-pills.html




