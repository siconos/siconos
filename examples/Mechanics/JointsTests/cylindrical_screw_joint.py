
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.kernel as Kernel

# A demonstration of how to couple the two free axes of a
# CylindricalJointR in order to construct a screw relation. (Coupled
# rotational and translational motion.)

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    # Bouncy contact with the ground
    io.add_Newton_impact_friction_nsl('contact', mu=0.3, e=0.6)

    # Definition of a bar
    io.add_primitive_shape('Bar', 'Box', (0.2, 0.2, 1))
    io.add_object('bar', [Contactor('Bar')], [0,0,1], mass=1)

    # Definition of the ground
    io.add_primitive_shape('Ground', 'Box', (2, 3, 0.1))
    io.add_object('ground', [Contactor('Ground')], [0,0,-0.05])

    # Add a cylindrical joint with a coupling between its two degrees
    # of freedom with a ratio of 5.0 (rotation of 5 radians for every
    # translation of 1.0 units)
    io.add_joint('joint1', 'bar', None, [[0,0,0]], [[0,0,1]], 'CylindricalJointR',
                coupled=[(0,1,5.0)], absolute=True)

# Load and run the simulation
with Hdf5(mode='r+') as io:
    io.run(t0=0,
           T=3,
           h=0.001,
           theta=0.5,
           Newton_max_iter=1,
           itermax=1000,
           tolerance=1e-12,
           projection_itermax=3,
           projection_tolerance=1e-5,
           projection_tolerance_unilateral=1e-5,
           time_stepping=Kernel.TimeSteppingDirectProjection,
           osi=Kernel.MoreauJeanDirectProjectionOSI,
    )
