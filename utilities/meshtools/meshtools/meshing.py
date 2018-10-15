from abaqus import *
from abaqusConstants import *
import mesh

import itertools
import numpy as np
import operator
import abc

from points import *


class Mesh(object):
    '''
    Container for items related to meshing
    '''


    @staticmethod
    def sortEdgesForBiasSeed(part, edges, center):
        '''
        Arguments:
        part: reference to the part object
        edges: a tuple of edges 
        center: a Point instance indicating the center of the fine mesh region

        This function returns the tuples e1 and e2 which corresponds to the arguments needed in seedEdgeByBias()
        where end1Edges = e1 and end2Edges = e2. In the seedEdgeByBias(), the smallest elements are 
        positioned near the first vertice on the edge for edges passed in to end1Edges and visa versa.
        '''

        # Check that e is a tuple of edges TODO
        # Check that c is a tuple of three floats TODO

        e1 = list()
        e2 = list()

        for e in edges:
            (v1, v2) = e.getVertices()
            v1Distance = Point.distance(Point.fromVertex(part.vertices[v1]), center)
            v2Distance = Point.distance(Point.fromVertex(part.vertices[v2]), center)

            if v1Distance < v2Distance:
                e1.append(e)
            else:
                e2.append(e)

        return (tuple(e1), tuple(e2))


    @staticmethod
    def elementsFromNodeArray(nodeArray):
        '''
        Given a <meshNodeArray> this function returns a list of labels of the elements 
        that are connected to the nodes in <meshNodeArray>
        '''
        # TODO: Check that the argument is a node Array

        # Get list of tuples of elements connected to each node (contains duplicates)
        allEls = [node.getElements() for node in nodeArray]

        # Flatten into a single list
        flattened = list(itertools.chain(*allEls))

        # Remove duplicates and return a list of the labels
        return list({x.label for x in flattened})

