

SYMLINKS=$(git ls-files -s | gawk '/120000/{print $4}')
for i in ${SYMLINKS}
do
  BASE_DIR=$(dirname $i)
  SLINK=$(cat $i)
  TARGET=${BASE_DIR}/${SLINK}
  if test -d ${TARGET}; then
    rm -rf $i
	cp -ra ${TARGET} $i
  else
    rm -f $i
	cp -a ${TARGET} $i
  fi
done