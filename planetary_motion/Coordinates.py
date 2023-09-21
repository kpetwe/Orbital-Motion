from math import *
from . import Trig

"""
Series of methods that allow us to convert between various coordinate systems when
calculating / presenting orbit data.
Conversions mentioned include:
    Rectangular - Spherical coordinates
    Equatorial - Ecliptic coordinates

This class also allows us to derive the actual Ecliptic and Spherical coordinates
from the orbital elements.
"""
"""
    Derivations of the Ecliptical and Spherical coordinates using the orbital elements.
"""
def EclipticCoords(r, N, v, w, i):
    xeclip = r * ( Trig.cosd(N) * Trig.cosd(v+w) - Trig.sind(N) * Trig.sind(v+w) * Trig.cosd(i) )
    yeclip = r * ( Trig.sind(N) * Trig.cosd(v+w) + Trig.cosd(N) * Trig.sind(v+w) * Trig.cosd(i) )
    zeclip = r * Trig.sind(v+w) * Trig.sind(i)
    return xeclip, yeclip, zeclip

"""
    Coordinate conversions between rectangular coordinates and spherical cooinates.
    
    Conversions can also apply for ecliptic and horizontal coordinates. 
    Replace RA, Decl with:
        Long, Lat (Ecliptic)
        Azimuth, Altitude (Horizontal)
"""
# Spherical -> Rectangular Coordinates
def rect(RA, Decl, r = 1):
    x = r * Trig.cosd(RA) * Trig.cosd(Decl)
    y = r * Trig.sind(RA) * Trig.cosd(Decl)
    z = r * Trig.sind(Decl)
    return x, y, z

# Rectangular -> Spherical Coordinates
def sphere(x, y, z):
    r = sqrt(x * x + y * y + z * z)
    RA = Trig.atan2d(y, x)
    Decl = Trig.atan2d(z, sqrt(x * x + y * y))
    return RA, Decl, r


"""
    The following methods are for converting between equatorial and ecliptic coordinate systems.

    Rotation of the coordinate systems is done along the x-axis (which points to
    the Vernal Point-- the common point of origin in ecliptic and equatorial coordinates)
    and is done through the angle of the obliquity of the ecliptic (~23.4 degrees).
"""
# Equatorial -> Ecliptic Conversion
def eqEc(xequat, yequat, zequat, oblecl = 23.4):
    xeclip = xequat
    yeclip = yequat * Trig.cosd(-oblecl) - zequat * Trig.sind(-oblecl)
    zeclip = yequat * Trig.sind(-oblecl) + zequat * Trig.cosd(-oblecl)
    return xeclip, yeclip, zeclip
    

# Ecliptic -> Equatorial Conversion
def ecEq(xeclip, yeclip, zeclip, oblecl = 23.4):
    xequat = xeclip
    yequat = yeclip * Trig.cosd(oblecl) - zeclip * Trig.sind(oblecl)
    zequat = yeclip * Trig.sind(oblecl) + zeclip * Trig.cosd(oblecl)
    return xequat, yequat, zequat

# TODO: Test
"""
    TESTING -- Azimuthal coordinates conversion
    - converts RA, Decl -> az, alt
    - correct for alt
"""
def azimuthal(RA, Decl, LST, lat):
    # -180 to +180 degrees-- 0 = object directly south,
    # negative = east of south, positive - west of south
    # outside of interval-- rev
    HA = LST - RA
    print('HA = ', HA)
    x,y,z = rect(HA, Decl)
    xhor = x * Trig.sind(lat) - z * Trig.cosd(lat)
    yhor = y
    zhor = x * Trig.cosd(lat) + z * Trig.sind(lat)

    az, alt, r = sphere(xhor, yhor, zhor)
    az += 180

    return az, alt 
