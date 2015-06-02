#view Param
mbtb.MBTB_MAX_BODIES_NUMBER 
mbtb.MBTB_MAX_CONTACTS_NUMBER
mbtb.MBTB_MAX_JOINTS_NUMBER

bodyDraw=array.array('I',[1 for i in range(mbtb.MBTB_MAX_BODIES_NUMBER)])
bodyTrans=array.array('d',[2.5 for i in range(mbtb.MBTB_MAX_BODIES_NUMBER)])
initVel=numpy.array([(0,0,0,0,0,0) for i in range(mbtb.MBTB_MAX_BODIES_NUMBER)])
               
#view parameters
ContactArtefactLength=10.0
ArtefactThershold=1e-7
NominalForce=0
contactDraw1=array.array('I',[ 1 for i in range(mbtb.MBTB_MAX_CONTACTS_NUMBER)])
contactDraw2=array.array('I',[ 1 for i in range(mbtb.MBTB_MAX_CONTACTS_NUMBER)])

contactTrans1=array.array('d',[2.7  for i in range(mbtb.MBTB_MAX_CONTACTS_NUMBER)])
contactTrans2=array.array('d',[2.7  for i in range(mbtb.MBTB_MAX_CONTACTS_NUMBER)])

NBARTEFACTS=0
contactNormalFromFace1=array.array('I',[ 1 for i  in range(mbtb.MBTB_MAX_CONTACTS_NUMBER)])
contactOffsetP1=array.array('I',[1 for i in range(mbtb.MBTB_MAX_CONTACTS_NUMBER)])


#3D parameters
with3D=1
freqUpdate=1
freqOutput=1
dumpGraphic=0
drawMode=mbtb.MBTB_ARTEFACT_REACTION+mbtb.MBTB_ARTEFACT_NORMAL+mbtb.MBTB_ARTEFACT_P1P2
#drawMode=mbtb.MBTB_ARTEFACT_P1P2
#Simulation parameters
stepSize=1e-4
stepNumber=60000
TSdeactivateYPosThreshold=1e-4
TSdeactivateYVelThreshold=0.0
TSactivateYPosThreshold=0.0
TSactivateYVelThreshold=100

TSProjectionMaxIteration=20
TSConstraintTol=1e-06
TSConstraintTolUnilateral=1e-08
TSLevelOfProjection=1
#solver parameters
withProj=0
withReduced=2

gotoPos=0
solverTol=1e-4
solverIt=100000

et=0.0

# link with dylib
apple=0
