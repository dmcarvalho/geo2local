# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 17:51:53 2018

@author: diego
"""

from math import sin, cos, sqrt, pow, asin, tan, radians, degrees, atan2


ellipsoids = {"airy": { "a": 6377563.396, "b": 6356256.91, "description": "Airy 1830"},
                "mod_airy": { "a": 6377340.189, "b": 6356034.446, "description": "Modified Airy"},
                "clrk66": { "a": 6378206.4, "b": 6356583.8, "description": "Clarke 1866"},
                "new_intl": { "a": 6378157.5, "b": 6356772.2, "description": "New International 1967"},
                "plessis": { "a": 6376523, "b": 6355863, "description": "Plessis 1817 (France)"},
                "SEasia": { "a": 6378155, "b": 6356773.3205, "description": "Southeast Asia"},
                "walbeck": { "a": 6376896, "b": 6355834.8467, "description": "Walbeck"},
                "sphere": { "a": 6370997, "b": 6370997, "description": "Normal Sphere (r6370997) "},
                "MERIT": { "a": 6378137, "b": 6356752.29821597, "description": "MERIT 1983"},
                "SGS85": { "a": 6378136, "b": 6356751.30156878, "description": "Soviet Geodetic System 85"},
                "GRS80": { "a": 6378137, "b": 6356752.31414036, "description": "GRS 1980(IUGG 1980)"},
                "IAU76": { "a": 6378140, "b": 6356755.28815753, "description": "IAU 1976"},
                "APL4.9": { "a": 6378137, "b": 6356751.79631182, "description": "Appl. Physics. 1965"},
                "NWL9D": { "a": 6378145, "b": 6356759.76948868, "description": "Naval Weapons Lab. 1965"},
                "andrae": { "a": 6377104.43, "b": 6355847.41523333, "description": "Andrae 1876 (Den. Iclnd.)"},
                "danish": { "a": 6377019.2563, "b": 6355762.52544567, "description": "Andrae 1876 (Denmark Iceland)"},
                "aust_SA": { "a": 6378160, "b": 6356774.71919531, "description": "Australian Natl & S. Amer. 1969"},
                "GRS67": { "a": 6378160, "b": 6356774.51609072, "description": "GRS 67(IUGG 1967)"},
                "GSK2011": { "a": 6378136.5, "b": 6356751.7579556, "description": "GSK-2011"},
                "bessel": { "a": 6377397.155, "b": 6356078.96281819, "description": "Bessel 1841"},
                "bess_nam": { "a": 6377483.865, "b": 6356165.38296633, "description": "Bessel 1841 (Namibia)"},
                "clrk80": { "a": 6378249.145, "b": 6356514.96582849, "description": "Clarke 1880 mod."},
                "clrk80ign": { "a": 6378249.2, "b": 6356515, "description": "Clarke 1880 (IGN)."},
                "CPM": { "a": 6375738.7, "b": 6356666.22191211, "description": "Comm. des Poids et Mesures 1799"},
                "delmbr": { "a": 6376428, "b": 6355957.92616372, "description": "Delambre 1810 (Belgium)"},
                "engelis": { "a": 6378136.05, "b": 6356751.32272154, "description": "Engelis 1985"},
                "evrst30": { "a": 6377276.345, "b": 6356075.41314024, "description": "Everest 1830"},
                "evrst48": { "a": 6377304.063, "b": 6356103.03899315, "description": "Everest 1948"},
                "evrst56": { "a": 6377301.243, "b": 6356100.2283681, "description": "Everest 1956"},
                "evrst69": { "a": 6377295.664, "b": 6356094.6679152, "description": "Everest 1969"},
                "evrstSS": { "a": 6377298.556, "b": 6356097.5503009, "description": "Everest (Sabah & Sarawak)"},
                "fschr60": { "a": 6378166, "b": 6356784.28360711, "description": "Fischer (Mercury Datum) 1960"},
                "fschr60m": { "a": 6378155, "b": 6356773.32048274, "description": "Modified Fischer 1960"},
                "fschr68": { "a": 6378150, "b": 6356768.33724439, "description": "Fischer 1968"},
                "helmert": { "a": 6378200, "b": 6356818.16962789, "description": "Helmert 1906"},
                "hough": { "a": 6378270, "b": 6356794.34343434, "description": "Hough"},
                "intl": { "a": 6378388, "b": 6356911.94612795, "description": "International 1909 (Hayford)"},
                "krass": { "a": 6378245, "b": 6356863.01877305, "description": "Krassovsky 1942"},
                "kaula": { "a": 6378163, "b": 6356776.99208691, "description": "Kaula 1961"},
                "lerch": { "a": 6378139, "b": 6356754.29151034, "description": "Lerch 1979"},
                "mprts": { "a": 6397300, "b": 6363806.28272251, "description": "Maupertius 1738"},
                "PZ90": { "a": 6378136, "b": 6356751.36179569, "description": "PZ-90"},
                "WGS60": { "a": 6378165, "b": 6356783.28695944, "description": "WGS 60"},
                "WGS66": { "a": 6378145, "b": 6356759.76948868, "description": "WGS 66"},
                "WGS72": { "a": 6378135, "b": 6356750.52001609, "description": "WGS 72"},
                "WGS84": { "a": 6378137, "b": 6356752.31424518, "description": "WGS 84"}}

class Ellipsoid(object):
    def __init__(self, *args, **kwargs):
        def get_ab(kwargs):
            if kwargs.has_key("a") and kwargs.has_key("b"):
                self._a = float(kwargs.get("a"))
                self._b = float(kwargs.get("b"))
            
        if kwargs.has_key("name"):
            self._name = kwargs.get("name")
            if ellipsoids.has_key(self._name):
                el = ellipsoids[self._name]
                self._a = el["a"]
                self._b = el["b"]
            else:
                pass        
        
            
    @property
    def name(self):
        return self._name
        
    '''
    The linear parameters
    '''    
    @property
    def a(self):
        return self._a
   
    @property
    def b(self):
        return self._b

    @property
    def ra(self):
        return 1./self._a
   
    @property
    def rb(self):
        return 1./self._b
        
        
    '''
    The eccentricities
    '''
    @property
    def alpha(self):
        return asin(self.e)
   
    @property
    def e(self):
        return sqrt( (pow(self.a, 2) - pow(self.b, 2)) / pow(self.a, 2) )
        
    @property
    def es(self):
        return pow(self.e, 2)
   
    @property
    def e2(self):
        return tan(self.alpha)
        
    @property
    def e2s(self):
        return pow(self.e2, 2)
        
    @property
    def e3(self):
        if (0 != self.alpha):
            return sin(self.alpha) / sqrt(2 - sin(self.alpha)*sin(self.alpha))
        else:
            return 0
                
    @property
    def e3s(self):
        return pow(self.e3, 2)
        
        
    @property
    def one_es(self):
        return 1. - self.es
        
        
    @property
    def rone_es(self):
        return 1./self.one_es
        
        
    '''
    The flattenings
    '''
    @property
    def f(self):
        return 1. - cos(self.alpha)
   
    @property
    def f2(self):
        return sqrt( (pow(self.a, 2) - pow(self.b, 2)) / pow(self.a, 2) )
        
    @property
    def n(self):
        return pow (tan(self.alpha/2), 2)
   
    @property
    def rf(self):
        return 1./self.f
        
    @property
    def rf2(self):
        return 1./self.f2
        
    @property
    def rn(self):
        return 1./self.n
        
    
    '''
    Transformations
    '''
    def normalRadiusOfCurvature(self, phi):
        return  self.a / sqrt(1 - pow(self.e,2) * pow(sin(phi),2))
    
    def hypot(self, x, y):
        return sqrt(pow(x,2) + pow(y,2)) 
        
    def geocentric_radius(self, phi):
         return self.hypot(pow(self.a,2)*cos(phi), pow(self.b,2)*sin(phi)) \
         / self.hypot(self.a*cos(phi), self.b*sin(phi));
    
    
    def geodesic2Cartesian(self, phi, _lambda, h):
        n = self.normalRadiusOfCurvature(radians(phi))
        x = (n+h) * cos(radians(phi)) * cos(radians(_lambda))
        y = (n+h) * cos(radians(phi)) * sin(radians(_lambda))
        z = (n * (1 - pow(self.e,2)) + h) * sin(radians(phi))
        return x,y,z
    
    
    def cartesian2Geodesic(self, x, y, z):
        p = self.hypot(x,y)        
        theta = atan2(z * self.a, p * self.b)
        c = cos(theta)
        s = sin(theta)
        phi = atan2 (z + (self.e2s * self.b * pow(s,3)), p - (self.es *  self.a * pow(c,3)));
        lam =  atan2 (y, x)
        s = sin(phi)    
        N = self.normalRadiusOfCurvature(phi)
        c = cos(phi)
        if abs(c) < 1e-6:
            r = self.geocentric_radius(phi)
            h = abs(z) - r
        else:
            h = (p / c) - N
        return  degrees(phi), degrees(lam), h