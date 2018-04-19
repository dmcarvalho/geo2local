# -*- coding: utf-8 -*-



semi_major = 6378137
semi_minor = 6356752.3141

phi0 = -29.68511910
lambda0 = -53.80338140
h0 = 135.788

phi1 = -29.65460913
lambda1 = -53.83412740
h1 = 450.118


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
                    [0, sin(radians(self.origin['phi')), cos(radians(self.origin['phi')))],
                    [0, -cos(radians(self.origin['phi']))), sin(radians(self.origin['phi')))]])
                    
    @property
    def b(self):
        return np.array([[-sin(radians(self.origin['lambda'])), cos(radians(self.origin['lambda'])), 0],
                    [-cos(radians(self.origin['lambda'])), -sin(radians(self.origin['lambda'])), 0],
                    [0, 0, 1]])
    

    @property
    def cartesianOrigin(self):
        return self.geodesic2Cartesian(self.origin)


    def geodesic2Cartesian(self, point):
        return self.ellipsoid
                    .convertGeodesic2Cartesian(point['phi'), point['lambda'], point['h'))



    def geo2local(self, point):
        x0,y0,z0 = self.cartesianOrigin()
        x,y,z = self.geodesic2Cartesian(point)         
                   
        c = np.array([[x-x0],
                    [y-y0],
                    [z-z0]])

        return self.a.dot(self.b).dot(c)

    def local2geo(self, point):
        x0,y0,z0 = self.cartesianOrigin()
        
        cartesian  = np.array([
                            [point['phi')], 
                            [point['lambda']], 
                            [point['h')]])

        e = inv(a.dot(b)).dot(cartesian)
        f = np.array([[x0],[y0],[z0]])
        return e+f
