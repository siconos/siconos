from math import cos, sin
from numpy.linalg import norm


class Shape(object):
    """
    A shape with translation and orientation. The given position is meant to be
    relative to some body reference frame.
    """
    def __init__(self,
                 shape_name,
                 relative_translation=[0, 0, 0],
                 relative_orientation=[1, 0, 0, 0]):

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
        self.translation = relative_translation
        self.orientation = ori


class Avatar(Shape):
    """
    An Avatar shape is a shape that may be associated to some  body
    for visualisation purpose only.
    """

    def __init__(self,
                 shape_name,
                 parameters=None,
                 relative_translation=[0, 0, 0],
                 relative_orientation=[1, 0, 0, 0]):

        self.parameters = parameters
        super(Avatar, self).__init__(shape_name,
                                     relative_translation,
                                     relative_orientation)


class Contactor(Shape):
    """
    A Contactor is associated to a shape. It belongs to a
    collision group. Depending on the geometrical engine used,
    some informations may be added to the contactor, such as
    the kind of contact and the concerned part of the shape.
    Note that contact laws must then be defined between collision
    groups.
    """

    def __init__(self,
                 shape_name,
                 collision_group=0,
                 contact_type=None,
                 contact_index=None,
                 relative_translation=[0, 0, 0],
                 relative_orientation=[1, 0, 0, 0]):

        self.group = collision_group
        self.contact_type = contact_type
        self.contact_index = contact_index

        super(Contactor, self).__init__(shape_name,
                                        relative_translation,
                                        relative_orientation)
