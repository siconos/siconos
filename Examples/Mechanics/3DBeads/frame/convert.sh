#!/bin/sh

# Convert ppm files into jpg
for f in *ppm ; do convert -quality 100 $f `basename $f ppm`jpg; done

