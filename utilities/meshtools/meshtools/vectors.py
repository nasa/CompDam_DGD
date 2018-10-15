from abaqus import *
from abaqusConstants import *

import numpy as np

from points import *

class Vector(object):
    '''
    A vector object
    '''

    def __init__(self, startPoint, endPoint):
        '''
        Expects a start Point and an end Point
        '''

        self.startPoint = startPoint
        self.endPoint = endPoint
        self.components = endPoint.coordinates - startPoint.coordinates

    @classmethod
    def fromComponents(cls, i, j, k):
        return cls(startPoint=Point(), endPoint=Point(x=i, y=j, z=k))

    @classmethod
    def fromNodes(cls, startNode, endNode):
        return cls(startPoint=Point.fromMeshNode(startNode), endPoint=Point.fromMeshNode(endNode))

    def __str__(self):
        return "(" + str(self.i) + "i, " + str(self.j) + "j, " + str(self.k) + "k)"

    @property
    def unit(self):
        return self.components/self.magnitude

    @property
    def i(self):
        return self.components[0]
    
    @property
    def j(self):
        return self.components[1]

    @property
    def k(self):
        return self.components[2]

    @property
    def magnitude(self):
        return Point.distance(self.startPoint,self.endPoint)

    @staticmethod
    def cross(a, b):
        """Returns the corss product of two vectors a and b"""
        return np.cross(a.components, b.components)

    @staticmethod
    def dot(a,b):
        """Returns the dot product of two vectors a and b"""
        return np.dot(a.components, b.components)

    @staticmethod
    def angle(a,b):
        """Returns the angle between two vectors a and b in radians"""
        if a.magnitude == 0.0 or b.magnitude == 0.0:
            return np.pi
        angle = acos(float(Vector.dot(a,b)/(a.magnitude*b.magnitude)))
        if np.isnan(angle):
            if (a == b).all():
                return 0.0
            else:
                return np.pi
        return angle

    def principalDirections(a):
        """Returns a VectorArray with the two principal directions (ie 2D) for a given vector, a"""
        vecs = list()
        vecs.append(Vector.fromComponents(1,0,0) if a.i > 0 else Vector.fromComponents(-1,0,0))
        vecs.append(Vector.fromComponents(0,1,0) if a.j > 0 else Vector.fromComponents(0,-1,0))
        return VectorArray(vecs)


class VectorArray(object):
    '''
    A container for vectors
    '''

    def __init__(self, vectors):
        '''
        Expects a list of vectors
        '''
        if len(vectors) < 1:
            raise ValueError("Expect at least one vector")
        self.vectors = vectors

    def __str__(self):
        out = ""
        for v in self.vectors:
            out = out + str(v) + "\n"
        return out

    def __getitem__(self, index):
        return self.vectors.__getitem__(index)

    def __setitem__(self, index):
        return self.vectors.__setitem__(index)

    def __delitem__(self, index):
        return self.vectors.__delitem__(index)

    @classmethod
    def fromSetOfPoints(cls, startPoint, endPoints):
        '''
        Creates a vector from the startPoint to each endPoint where endPoints is a PointArray
        '''
        vectors = list()
        for endPoint in endPoints.points:
            vectors.append(Vector(startPoint=startPoint, endPoint=endPoint))
        return cls(vectors=vectors)
