
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.kernel as Kernel
import numpy as np

# A demonstration of how to couple the two free axes of a
# CylindricalJointR in order to construct a screw relation. (Coupled
# rotational and translational motion.)

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    # Bouncy contact with the ground
    io.add_Newton_impact_friction_nsl('contact', mu=0.3, e=0.9)

    c = np.cos(np.pi/4)
    s = np.sin(np.pi/4)

    # Definition of a long and short bar
    io.add_primitive_shape('LongBar', 'Box', (0.2, 0.2, 1))
    io.add_primitive_shape('ShortBar', 'Box', (0.2, 0.2, 0.5))

    ## Rotating gear ratio demonstration

    # Put two bars at equal 45' angles and let them drop
    io.add_object('bar1', [Contactor('LongBar')], [0,-1,0.5], mass=1)
    io.add_object('bar2', [Contactor('LongBar')], [c*0.4,-0.8,1+s*0.3], mass=1,
                 orientation=[(0,1,0),np.pi/4])
    io.add_object('bar3', [Contactor('ShortBar')], [-c*0.15,-0.6,1+s*0.05], mass=1,
                 orientation=[(0,1,0),-np.pi/4])

    # Definition of the ground
    io.add_primitive_shape('Ground', 'Box', (2, 3, 0.1))
    io.add_object('ground', [Contactor('Ground')], [0,0,-0.05])

    # Fix bar1 in place
    io.add_joint('jnt1', 'bar1', None, None, None, 'FixedJointR')

    # Add a pivot joint from bar2 to bar1 and another from bar1 to bar3.
    io.add_joint('jnt2','bar2','bar1',[[0,0,0.9]],[[0,1,0]],'PivotJointR',absolute=True)
    io.add_joint('jnt3','bar1','bar3',[[0,0,0.9]],[[0,1,0]],'PivotJointR',absolute=True)

    # We specified a gear ratio of 2.0 between dof 0 of jnt2 and dof
    # 0 of jnt3.  The short bar must maintain an angle
    # theta3=2.0*theta2, and it therefore spins twice as fast.
    io.add_joint('jnt4','bar2','bar3',None,None,'CouplerJointR',
                coupled=[[0, 0, 2.0]], references=['jnt2','jnt3','bar1'])

    ## Rotation-prismatic gear ratio demonstration

    # Put two bars at equal 45' angles and let them drop
    io.add_object('bar4', [Contactor('LongBar')], [0,0.6,0.5], mass=1)
    io.add_object('bar5', [Contactor('LongBar')], [0,0.8,1], mass=1)
    io.add_object('bar6', [Contactor('ShortBar')], [0,1.0,1], mass=1
                 , velocity=[0,0,0,10,10,10])

    # Fix bar1 in place
    io.add_joint('jnt5', 'bar4', None, None, None, 'FixedJointR')

    # Add a slider joint from bar5 to bar4
    io.add_joint('jnt6','bar4','bar5',None,[[0,0,1]],'PrismaticJointR',absolute=True)

    # Add a pivot joint from bar6 to bar5
    io.add_joint('jnt7','bar5','bar6',[[0,0.6,1]],[[0,1,0]],'PivotJointR',absolute=True)

    # We specified a gear ratio of 0.2 between dof 0 of jnt6 and dof
    # 0 of jnt7.  Since jnt6 must maintain an angle
    # theta6=2.0*dist7, bar6 therefore spins as bar5 drops.
    io.add_joint('jnt8','bar4','bar6',None,None,'CouplerJointR',
                coupled=[[0, 0, 3.0]], references=['jnt6','jnt7','bar5'])

# Load and run the simulation
with Hdf5(mode='r+') as io:
    io.run(t0=0,
           T=5,
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
