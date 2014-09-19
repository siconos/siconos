from math import cos, sin
from numpy.linalg import norm

class Contactor:
    def __init__(self,
                 shape_name,
                 collision_group=0,
                 relative_position=[0, 0, 0],
                 relative_orientation=[0, 1, 0, 0]):

        if len(relative_orientation) == 2:
            # axis + angle
            axis = relative_orientation[0]
            assert len(axis) == 3
            angle = relative_orientation[1]
            assert type(angle) is float
            n = sin(angle / 2) / norm(axis)

            ori = [cos(angle / 2.), axis[0] * n, axis[1] * n, axis[2] * n]
        else:
            # a given quaternion
            assert len(relative_orientation) == 4
            ori = relative_orientation

        self.name = shape_name
        self.group = collision_group
        self.position = relative_position
        self.orientation = ori

