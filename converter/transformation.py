# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 16:20:00 2018

@author: diego
"""


import numpy as np
from numpy.linalg import inv


from math import sin, cos, radians, tan, sqrt




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


class TransformationNBR14166(object):
    def __init__(self, ellipsoid, origin):
        self._ellipsoid = ellipsoid
        self._origin = origin

    @property
    def ellipsoid(self):
        return self._ellipsoid

    @property
    def origin(self):
        return self._origin


    def geo2local(self, point):
        _1seg = (1./3600)

        y0,x0,z0 = [radians(i) for i in self.origin]
        y1,x1,z1 = [radians(i) for i in point]
        Ht = 150

        N0 = self.ellipsoid.normalRadiusOfCurvature(y0)
        N1 = self.ellipsoid.normalRadiusOfCurvature(y1)


        M0 = (self.ellipsoid.a * ( 1. - self.ellipsoid.es) ) / ((1. - self.ellipsoid.es * sin(y0)**2.)**(3./2.))

        R0 = sqrt(M0 * N0)

        c = (R0 + Ht) / R0

        E = (1. + 3. * tan(y0)) / (6. * N0**2)

        D = (3. * self.ellipsoid.es * sin(y0) * cos(y0) * (_1seg)) / (2. * (1. - self.ellipsoid.es * sin(y0)**2))

        C = tan(y0) / (2 * M0 * N0 * (_1seg))

        B = 1. / (M0*(_1seg))

        deltaphi2 = (y1 - y0) * 3600
        deltalam2 = (x1 - x0) * 3600

        deltaphi1 = deltaphi2 * (1. - 3.9173 * 10**-12* (deltaphi2**2))

        deltalam1 = deltalam2 * (1. - 3.9173 * 10**-12* (deltalam2**2))

        xp = (deltalam1) *  cos(y1) * N1 * _1seg * c


        var1 = (deltaphi1 + C * xp**2 + D * (deltaphi1)**2 + E * xp**2* (deltaphi1) + E * C * xp**4)

        yp = 1/B * ((var1) * c)

        return [yp,xp,0]

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