#!/bin/sh
case $1 in -x) set -x; shift ;; esac

#set -e
#set -x

# Usage :: cd /path/to/siconos; sh Build/tools/create-source-pkg.sh;
# archives are created in DESTDIR = BASE_DESTDIR-VERSION
# You can tune the BASE_DESTDIR variable
# You have to checkout at the right git branch/commit/tag

if [ "$#" -ne 1 ]; then
	echo "You have to specify one (1) argument: the version of siconos we are packaging"
	exit 1
fi

VERSION=$1

BASE_DESTDIR=/tmp/`whoami`/siconos-all-in-one-source

DESTDIR=${BASE_DESTDIR}-${VERSION}
if [ -e ${DESTDIR} ]; then
	rm -rf ${DESTDIR}
fi

git rev-parse --is-inside-work-tree >/dev/null 2>&1 || { echo "Not in a git repo, exiting"; exit 1; }

TOPLEVEL=$(git rev-parse --show-toplevel)

cd ${TOPLEVEL}
for i in $(git ls-files)
do
	mkdir -p ${DESTDIR}/$(dirname $i)
	if test -L $i; then
		slink=$(readlink $i);
		slink_dir=$(dirname $i)
		FILES_TO_COPY=$(git ls-files ${slink_dir}/$slink)
		# dirname of a symlink, even if the latter is pointing to a directory,
		# is the directory of the symlink
		if test -d $i; then
			mkdir -p ${DESTDIR}/$i
		fi
		cp -a ${FILES_TO_COPY} ${DESTDIR}/$i
	else
		cp -Pa $i ${DESTDIR}/$i
	fi
done
DESTDIR_PARDIR=$(dirname ${DESTDIR})
DESTDIR_NAME=$(basename ${DESTDIR})
cd ${DESTDIR_PARDIR}
tar zcvf ${DESTDIR_NAME}.tar.gz ${DESTDIR_NAME}
zip -9 ${DESTDIR_NAME}.zip -r ${DESTDIR_NAME}
