

# Generate this list with
# git ls-files -s | gawk '/120000/{print $4}'
SYMLINKS="Build/Pipol/rc.fedora-core14
Build/Pipol/rc.fedora-core17
control/src/tests/PID.ref
control/src/tests/TestMain.cpp
examples/Mechanics/JointsTests/NE_1DS_1Knee_MLCP.cmake
examples/Mechanics/JointsTests/NE_1DS_1Knee_MLCP_MoreauJeanCombinedProjection.cmake
examples/Mechanics/JointsTests/NE_3DS_3Knee_1Prism_GMP.cmake
examples/Mechanics/JointsTests/NE_3DS_3Knee_1Prism_MLCP.cmake
examples/Mechanics/JointsTests/NE_3DS_3Knee_1Prism_MLCP_MoreauJeanCombinedProjection.cmake
examples/Mechanics/JointsTests/NE_BouncingBeam.cmake
examples/Mechanics/JointsTests/PrismaticTest.cmake
examples/cmake"


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
