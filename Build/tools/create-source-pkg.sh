#!/bin/sh
case $1 in -x) set -x; shift ;; esac
#set -e


# Usage :: cd /path/to/siconos; sh Build/tools/create-source-pkg.sh;
# archives are created in DESTDIR = BASE_DESTDIR-VERSION
# You can tune the BASE_DESTDIR variable
# You have to checkout at the right git branch/commit/tag

if [ "$#" -ne 1 ]; then
	echo "You have to specify one (1) argument: the version of siconos we are packaging"
	exit 1
fi

VERSION=$1

BASE_DESTDIR=/tmp/`whoami`/siconos-source

MODULES="Numerics Kernel Mechanics Control IO Front-End Examples"

DESTDIR=${BASE_DESTDIR}-${VERSION}
if [ -e ${DESTDIR} ]; then
	rm -rf ${DESTDIR}
fi

git rev-parse --is-inside-work-tree >/dev/null 2>&1 || { echo "Not in a git repo, exiting"; exit 1; }

TOPLEVEL=$(git rev-parse --show-toplevel)

for MODULE in ${MODULES}
do
	cd ${TOPLEVEL}/${MODULE}
	MODULE_DESTDIR=${DESTDIR}/siconos-$(echo ${MODULE} | tr '[:upper:]' '[:lower:]')-${VERSION}-Source
	for i in $(git ls-files)
	do
		mkdir -p ${MODULE_DESTDIR}/$(dirname $i)
		if test -L $i; then
			slink=$(readlink $i);
			slink_dir=$(dirname $i)
			FILES_TO_COPY=$(git ls-files ${slink_dir}/$slink)
			# dirname of a symlink, even if the latter is pointing to a directory,
			# is the directory of the symlink
			if test -d $i; then
				mkdir -p ${MODULE_DESTDIR}/$i
			fi
			cp -a ${FILES_TO_COPY} ${MODULE_DESTDIR}/$i
		else
			cp -Pa $i ${MODULE_DESTDIR}/$i
		fi
	done
	MODULE_DESTDIR_PARDIR=$(dirname ${MODULE_DESTDIR})
	MODULE_DESTDIR_NAME=$(basename ${MODULE_DESTDIR})
	cd ${MODULE_DESTDIR_PARDIR}
        tar zcvf ${MODULE_DESTDIR_NAME}.tar.gz ${MODULE_DESTDIR_NAME}
	zip -9 ${MODULE_DESTDIR_NAME}.zip -r ${MODULE_DESTDIR_NAME}
done
