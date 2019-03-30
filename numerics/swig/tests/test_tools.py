import siconos.numerics as sn
import numpy as np

def test_projection():
    r = np.array([1.0,1.0,1.0])
    status = sn.projectionOnCone(r, 0.5)
    rr = np.array([1.0,1.0,1.0,1.0,1.0])
    status_r = sn.projectionOnRollingCone(rr, 0.5, 0.5)
    return status, status_r
