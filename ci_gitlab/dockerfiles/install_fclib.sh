#!/usr/bin/env bash

# -- Fclib download --
cd "$HOME"||exit 1
rm -rf "$HOME/fclib"
git clone https://github.com/FrictionalContactLibrary/fclib.git

# Creates a directory to build libs.
mkdir -p "$HOME/build"
cd "$HOME/build"||exit 1

# Configure and build fclib
cmake "$HOME/fclib" -DFCLIB_HEADER_ONLY=OFF
make
make install
rm -rf "$HOME/fclib" "$HOME/build"
