# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 16:20:00 2018

@author: diego
"""


import numpy as np
from numpy.linalg import inv


from math import sin, cos, sqrt, pow, radians, degrees, atan2, asin, tan


class TransformationRT(object):
    def __init__(self, ellipsoid, origin):
        self._ellipsoid = ellipsoid
        self._origin = origin

    @property
    def ellipsoid(self):
        return self._ellipsoid

    @property
    def origin(self):
        return self._origin

    @property
    def a(self):
        return np.array([[1, 0, 0],
                    [0, sin(radians(self.origin[0])), cos(radians(self.origin[0]))],
                    [0, -cos(radians(self.origin[0])), sin(radians(self.origin[0]))]])
                    
    @property
    def b(self):
        return np.array([[-sin(radians(self.origin[1])), cos(radians(self.origin[1])), 0],
                    [-cos(radians(self.origin[1])), -sin(radians(self.origin[1])), 0],
                    [0, 0, 1]])
    
    def cartesianOrigin(self):
        return self.geodesic2Cartesian(self.origin)


    def geodesic2Cartesian(self, point):
        return self.ellipsoid.geodesic2Cartesian(point[0], point[1], point[2])

    def cartesian2geodesic(self, point):
        return self.ellipsoid.cartesian2Geodesic(point[0], point[1], point[2])


    def geo2local(self, point):
        x0,y0,z0 = self.cartesianOrigin()
        x,y,z = self.geodesic2Cartesian(point)         
                   
        c = np.array([[x-x0],
                    [y-y0],
                    [z-z0]])

        return self.a.dot(self.b).dot(c).flatten().tolist()

    def local2geo(self, point):
        x0,y0,z0 = self.cartesianOrigin()
        
        cartesian  = np.array([
                            [point[0]], 
                            [point[1]], 
                            [point[2]]])

        e = inv(self.a.dot(self.b)).dot(cartesian)
        f = np.array([[x0],[y0],[z0]])
        g = e+f
        return self.cartesian2geodesic(g.flatten().tolist())