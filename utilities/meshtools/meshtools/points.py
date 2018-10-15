from abaqus import *
from abaqusConstants import *

import numpy as np

class Point(object):
    '''
    Representation for a point in 3D cartesian coordinates
    '''

    def __init__(self, x=0.0, y=0.0, z=0.0, label=0):
        self.coordinates = np.array([x, y, z], float)
        self.label = label

    @classmethod
    def fromNumpyArray(cls, numpyArray):
        if numpyArray.size == 3:
            return cls(x=numpyArray[0], y=numpyArray[1], z=numpyArray[2])
        else:
            raise TypeError("Expected numpy array with size 3")

    @classmethod
    def fromMeshNode(cls, meshNodeObject):
        coords = meshNodeObject.coordinates
        return cls(x=coords[0], y=coords[1], z=coords[2], label=meshNodeObject.label)

    @classmethod
    def fromVertex(cls, vertexObject):
        coords = vertexObject.pointOn[0]
        return cls(x=coords[0], y=coords[1], z=coords[2])

    @classmethod
    def fromMeshElementAtCentroid(cls, meshElement):
        '''Generates a point at the element centroid'''
        cPt = PointArray.fromMeshNodeArray(meshElement.getNodes()).center
        return cls(x=cPt.x, y=cPt.y, z=cPt.z, label=meshElement.label)

    def __str__(self):
        return "(" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ", Label: " + str(self.label) + ")"

    def __repr__(self):
        return "Point (" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ", Label: " + str(self.label) + ")"

    def __eq__(self, other):
        return all(self.coordinates == other.coordinates)

    def __add__(self, other):
        return self.fromNumpyArray(self.coordinates + other.coordinates)

    def __sub__(self, other):
        return self.fromNumpyArray(self.coordinates - other.coordinates)

    @property
    def x(self):
        return self.coordinates[0]

    @property
    def y(self):
        return self.coordinates[1]
    
    @property
    def z(self):
        return self.coordinates[2]

    @property
    def r(self):
        return sqrt((self.coordinates[0])**2 + (self.coordinates[1])**2 + (self.coordinates[2])**2)

    @staticmethod
    def distance(point1, point2):
        """Returns the distance between two points"""
        return sqrt((point1.x-point2.x)**2 + (point1.y-point2.y)**2 + (point1.z-point2.z)**2)


class PointArray(object):
    '''
    A sequence of points
    '''

    def __init__(self, points):
        '''
        Expects a list of points
        '''
        if len(points) < 1:
            raise ValueError("Expect at least one point")
        self.points = points

    def __str__(self):
        out = ""
        for p in self.points:
            out = out + str(p) + ", "
        return out

    def __iter__(self):
        return iter(self.points)

    @classmethod
    def fromMeshNodeArray(cls, meshNodeArray):
        points = list()
        for n in meshNodeArray:
            points.append(Point.fromMeshNode(n))
        return cls(points)

    @classmethod
    def fromMeshNodeLabels(cls, part, meshNodeLabels):
        return cls.fromMeshNodeArray(part.nodes.sequenceFromLabels(meshNodeLabels))

    @classmethod
    def fromVertexArray(cls, vertexArray):
        raise NotImplementedError()

    @classmethod
    def fromMeshElementArrayAtCentroids(cls, meshElementArray):
        '''Generates a point array of mesh element centroids'''
        points = list()
        for e in meshElementArray:
            cPt = PointArray.fromMeshNodeArray(e.getNodes()).center
            cPt.label = e.label
            points.append(cPt)
        return cls(points)

    @property
    def center(self):
        """Returns the point that is at the center of all the points in array"""

        numPoints = len(self.points)
        sumPoints = Point()

        for p in self.points:
            sumPoints = sumPoints + p

        return Point.fromNumpyArray(sumPoints.coordinates/numPoints)

    @property
    def labels(self):
        return tuple([p.label for p in self.points])

