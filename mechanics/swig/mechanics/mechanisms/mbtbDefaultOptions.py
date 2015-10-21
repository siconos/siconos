""""Default values for bodies and contact parameters.
Should be change in bodydef.py if needed.
"""
from siconos.mechanics.mechanisms import mbtb
import numpy as np
import array

mbtb.MBTB_MAX_BODIES_NUMBER
mbtb.MBTB_MAX_CONTACTS_NUMBER
mbtb.MBTB_MAX_JOINTS_NUMBER

initVel = np.array([(0, 0, 0, 0, 0, 0)
                    for i in range(mbtb.MBTB_MAX_BODIES_NUMBER)])
fctfext = np.array(['' for i in range(mbtb.MBTB_MAX_BODIES_NUMBER)])
fctmext = np.array(['' for i in range(mbtb.MBTB_MAX_BODIES_NUMBER)])
fctfint = np.array(['' for i in range(mbtb.MBTB_MAX_BODIES_NUMBER)])
fctmint = np.array(['' for i in range(mbtb.MBTB_MAX_BODIES_NUMBER)])
fctfintjacq = np.array(['' for i in range(mbtb.MBTB_MAX_BODIES_NUMBER)])
fctmintjacq = np.array(['' for i in range(mbtb.MBTB_MAX_BODIES_NUMBER)])
fctfintjacv = np.array(['' for i in range(mbtb.MBTB_MAX_BODIES_NUMBER)])
fctmintjacv = np.array(['' for i in range(mbtb.MBTB_MAX_BODIES_NUMBER)])

boundaryCondition = np.array(['' for i in range(mbtb.MBTB_MAX_BODIES_NUMBER)])
boundaryConditionIndex = np.array(
    [np.array([], dtype=int) for i in range(mbtb.MBTB_MAX_BODIES_NUMBER)])

NBARTEFACTS = 0

contactNormalFromFace1 = array.array(
    'I', [1 for i in range(mbtb.MBTB_MAX_CONTACTS_NUMBER)])
contactOffsetP1 = array.array(
    'I', [1 for i in range(mbtb.MBTB_MAX_CONTACTS_NUMBER)])

# View and rendering parameters
with3D = 1                   # rendering window

ContactArtefactLength = 10.0
ArtefactThreshold = 1e-7
NominalForce = 0

bodyDraw = array.array('I', [1 for i in range(mbtb.MBTB_MAX_BODIES_NUMBER)])
bodyTrans = array.array('d', [2.5 for i in range(mbtb.MBTB_MAX_BODIES_NUMBER)])

contactDraw1 = array.array(
    'I', [1 for i in range(mbtb.MBTB_MAX_CONTACTS_NUMBER)])
contactDraw2 = array.array(
    'I', [1 for i in range(mbtb.MBTB_MAX_CONTACTS_NUMBER)])

contactTrans1 = array.array(
    'd', [2.7 for i in range(mbtb.MBTB_MAX_CONTACTS_NUMBER)])
contactTrans2 = array.array(
    'd', [2.7 for i in range(mbtb.MBTB_MAX_CONTACTS_NUMBER)])

freqUpdate = 1    # frequency update of the view windows
freqOutput = 1    # frequency update of data output
dumpGraphic = 0   # if 1 dump graphic file

drawMode = mbtb.MBTB_ARTEFACT_REACTION + mbtb.MBTB_ARTEFACT_NORMAL \
    + mbtb.MBTB_ARTEFACT_P1P2
#drawMode = mbtb.MBTB_ARTEFACT_P1P2

#Default simulation parameters

stepSize = 1e-4
stepNumber = 60000

TSTheta = 0.5
TSGamma = 0.5

TSNewtonTolerance = 1e-10
TSNewtonMaxIteration = 15

TSdeactivateYPosThreshold = 1e-4
TSdeactivateYVelThreshold = 0.0
TSactivateYPosThreshold = 0.0
TSactivateYVelThreshold = 100

TSProjectionMaxIteration = 20
TSConstraintTol = 1e-06
TSConstraintTolUnilateral = 1e-08
TSLevelOfProjection = 1
#solver parameters
withProj = 0
withReduced = 2

gotoPos = 0
solverTol = 1e-4
solverIt = 100000

et = 0.0

# link with dylib
apple = 0
