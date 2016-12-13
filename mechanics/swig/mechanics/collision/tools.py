from math import cos, sin
from numpy.linalg import norm


class MovedShape(object):
    """A shape with a reference to some shape data and that is moved by a
    given relative translation and relative orientation. The
    orientation may be given as a quaternion or in the form (axis,
    angle) where len(axis) == 3 and angle is a scalar.

    The default value for the translation is
    [0, 0, 0] and [1, 0, 0, 0] for the orientation.
    """
    def __init__(self,
                 shape_data,
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

        self.data = shape_data
        self.translation = relative_translation
        self.orientation = ori


class Shape(MovedShape):
    """
    A MovedShape with optional instance name and parameters.
    """

    def __init__(self,
                 shape_data,
                 instance_name=None,
                 parameters=None,
                 relative_translation=[0, 0, 0],
                 relative_orientation=[1, 0, 0, 0]):

        self.instance_name = instance_name
        self.parameters = parameters
        super(Shape, self).__init__(shape_data,
                                    relative_translation,
                                    relative_orientation)


class Contactor(Shape):
    """A Contactor is a Shape that belongs to a collision
    group. Depending on the geometrical engine used, some informations
    may be added to the contactor, such as the kind of contact and the
    concerned part of the shape.  Note that contact laws must then be
    defined between collision groups.
    """

    def __init__(self,
                 shape_data,
                 instance_name=None,
                 parameters=None,
                 collision_group=0,
                 contact_type=None,
                 contact_index=None,
                 relative_translation=[0, 0, 0],
                 relative_orientation=[1, 0, 0, 0]):

        self.group = collision_group
        self.contact_type = contact_type
        self.contact_index = contact_index

        super(Contactor, self).__init__(shape_data,
                                        instance_name,
                                        parameters,
                                        relative_translation,
                                        relative_orientation)
