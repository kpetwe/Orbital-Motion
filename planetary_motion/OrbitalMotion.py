from math import *
from . import Coordinates
from . import OrbitalElements
from . import Perturbations
from . import Trig

"""
Calculates planetary positions of planets with an error of at most 1-2 arc minutes.
All formulae used in this file are courtesy of https://www.stjarnhimlen.se/comp/tutorial.html#4.
Contains methods for calculating between coordinate systems (elliptical and equatorial)
and for standard heliocentric and geocentric models.
"""

"""
    Computes the time value "d" using a date in the year.
"""
# d = 0.0 at 2000 Jan 0.0 or 1999 Dec 31.0 0:00 UT
# D = date, y = year in 4 digits, M = month (1-12)
def date(D, m, y, UT = 0.0):
    d = 367 * y - 7 * (y + (m + 9) // 12) // 4 - 3 * ((y + (m - 9) // 7 ) // 100 + 1) // 4 + 275 * m // 9 + D - 730515
    return d + UT / 24.0

# TODO: test
"""
    TESTING -- Sidereal Time
    long = local longitude
    East Longitude positive, West negative
"""
def sidereal(UT, Ls, long):
    GMST0 = Ls + 180
    GMST = GMST0 + UT
    LST = GMST + long
    return LST


"""
    Longitudinal Precession: since all calculations are relative to the 
    celestial equator / ecliptic center at the moment, this function converts
    to a standard epoch.
"""
def precession(d, epoch=2000.0):
    return 3.82394E-5 * ( 365.2422 * ( epoch - 2000.0 ) - d )


"""
    Calculates the Sun's position on day d
    Returns the ecliptic and equatorial coordinates in relation to earth
    of the sun, along with the mean longitude and mean anomaly
"""
def sunPos(d):
    w, a, e, M = OrbitalElements.sunElm(d)
    # Obliquity of the ecliptic (deg)
    oblecl = 23.4393 - 3.563E-7 * d
    # Sun's mean Longitude
    L = Trig.rev(w + M)
    # Auxiliary Angle (Eccentric Anomaly)
    E = M + (180/pi) * e * Trig.sind(M) * (1 + e * Trig.cosd(M))
    
    r, v = distAnomaly(a, E, e)
    # Compute Sun's ecliptic rectangular coords
    x, y, z = Coordinates.EclipticCoords(r, 0, v, w, 0)
    # Rotate coordinates using calculated oblecl
    xeq, yeq, zeq = Coordinates.ecEq(x, y, z, oblecl)

    # Convert to RA and Decl
    RA, Decl, r = Coordinates.sphere(xeq, yeq, zeq)
    return x, y, z, xeq, yeq, zeq, RA, Decl, r, L, M, oblecl


"""
    The following methods are used for calculating the moon and planets'
    positions on day d. It's important to note that the unit measurements described here
    for the moon are in earth radii, while the unit measurements for the other planets are in
    astronomical units (au).
"""
# Returns the ecliptic and spherical coords
# Precondition: n is an integer from 0 to 5 where:
# 0 - Moon, 1 - Mercury, 2 - Venus, 3 - Mars, 4 - Jupiter, 5 - Saturn
def getPlanPos(n, d):
    N, i, w, a, e, M = OrbitalElements.getPlan(n,d)

    E = eccentricAnomaly(M, e)
    r, v = distAnomaly(a, E, e)

    # Ecliptic Coords
    xeclip, yeclip, zeclip = Coordinates.EclipticCoords(r, N, v, w, i)

    # Convert to longitude, latitude, distance
    long, lat, r = Coordinates.SphericalCoords(xeclip, yeclip, zeclip)

    # Perturbations
    long, lat, r = Perturbations.perturPlanets(n, d, M, N, w, long, lat, r)
    xeclip, yeclip, zeclip = Coordinates.rect(long, lat, r)

    return xeclip, yeclip, zeclip, long, lat, r

"""
    Converts to geocentric position using the sun's rectangular coordinates
"""
def togeo(sunCoords, planCoords):
    xgeoc = sunCoords[0] + planCoords[0]
    ygeoc = sunCoords[1] + planCoords[1]
    zgeoc = sunCoords[2] + planCoords[2]
    return xgeoc, ygeoc, zgeoc

"""
    Given a date value and a planet, returns the geocentric 
    (equatorial) coordinates of that planet.
"""
def getPlanPosGeo(n, d):
    sun_arr = OrbitalElements.sunPos(d)
    plan_arr = OrbitalElements.getPlan(n, d)
    plan_geo = togeo(sun_arr, plan_arr)
    plan_geo_eq = Coordinates.ecEq(plan_geo[0], plan_geo[1], plan_geo[2], sun_arr[11])
    return plan_geo_eq


"""
    TODO: ADD description for this method
"""
def eccentricAnomaly(M, e):
    # Eccentric Anomaly Calculation
    E0 = M + (180/pi) * e * Trig.sind(M) * (1 + e * Trig.cosd(M))

    # Since Eccentricity of the Moon's orbit is larger than of the Earth's,
    # Our initial approximation will have a greater error than 0.005 degrees
    E1 = E0 - (E0 - (180/pi) * e * Trig.sind(E0) - M) / (1 - e * Trig.cosd(E0))
    while (abs(E0 - E1) > 0.0005):
        E0 = E1
        E1 = E0 - (E0 - (180/pi) * e * Trig.sind(E0) - M) / (1 - e * Trig.cosd(E0))
    return E1


"""
    TODO: Add description for this method
"""
def distAnomaly(a, E, e):
    # Compute rectangular coords in orbital plane
    x = a * (Trig.cosd(E) - e)
    y = a * sqrt(1 - e * e) * Trig.sind(E)
    
    # Convert to Distance and true anomaly
    r = sqrt(x*x + y*y) 
    v = Trig.rev(Trig.atan2d(y, x))
    return r, v
