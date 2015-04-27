from math import cos, sin
from numpy.linalg import norm

class Shape(object):
    """
    A shape is associated to some body with relatives position and orientation.
    """
    def __init__(self,
                 shape_name,
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
        self.position = relative_position
        self.orientation = ori


class Avatar(Shape):
    """
    An Avatar shape is a shape associated to a body for visualisation
    purpose only.
    """

    def __init__(self,
                 shape_name,
                 parameters=None,
                 relative_position=[0, 0, 0],
                 relative_orientation=[0, 1, 0, 0]):

        self.parameters = parameters
        super(Avatar, self).__init__(shape_name,
                                     relative_position,
                                     relative_orientation)


class Contactor(Shape):
    """
    A Contactor shape belongs to a group and is associated to a body.
    Contact laws must then be defined between groups.
    """

    def __init__(self,
                 shape_name,
                 collision_group=0,
                 relative_position=[0, 0, 0],
                 relative_orientation=[0, 1, 0, 0]):

        self.group = collision_group

        super(Contactor, self).__init__(shape_name,
                                        relative_position,
                                        relative_orientation)
