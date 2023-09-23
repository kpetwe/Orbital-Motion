from math import *
from . import Coordinates
from . import OrbitalElements
from . import Perturbations
from . import Trig

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
    Returns required data for sun
    d = day
    m = month
    y = year (4 digit number)
    ut = universal time (as a decimal)
    ct = coordinate type (either "rect" or "sphere")
    csys = coordinate system (either "eclip" or "equat")
        Note: spherical ecliptical = Longitude / Latitude / magnitude of distance
        and spherical equatorial = RA / Decl / magnitude of distance
"""
def getSunData(d = 31, m = 12, y = 1999, ut = 0.0, ct = "rect", csys = "eclip"):
    d = date(d, m, y, ut)
    return sunPos(d, ct, csys)

def sunPos(d = 0.0, ct = "rect", csys = "eclip"):    
    w, a, e, M = OrbitalElements.sunElm(d)
    # Obliquity of the ecliptic (deg)
    oblecl = 23.4393 - 3.563E-7 * d
    # Sun's mean Longitude
    L = Trig.rev(w + M)

    if (csys == "pertur"):
        return L, M, 0, oblecl

    # Auxiliary Angle (Eccentric Anomaly)
    E = M + (180/pi) * e * Trig.sind(M) * (1 + e * Trig.cosd(M))
    r, v = distAnomaly(a, E, e)

    # Compute Sun's ecliptic rectangular coords
    sun_arr = correctCoords(Coordinates.EclipticCoords(r, 0, v, w, 0), oblecl, ct, csys)
    return sun_arr[0], sun_arr[1], sun_arr[2], oblecl


"""
    Returns required data for a given planet (or the moon)
    d = day
    m = month
    y = year (4 digit number)
    ut = universal time (as a decimal)
    p = planet ("Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn") -- case insensitive
    ct = coordinate type (either "rect" or "sphere")
    csys = coordinate system (either "eclip" or "equat")
        (only applies to geocentric system-- heliocentric only uses the ecliptic coordinate system)
        Note for geocentric: spherical ecliptical = Longitude / Latitude / r
        and spherical equatorial = RA / Decl / r
    geoc = bool for whether you want geocentric or heliocentric coordinates**
    ** the moon's coordinates are only provided using the geocentric model
"""
def getPlanData(d = 31, m = 12, y = 1999, ut = 0.0, p = "mercury", ct = "rect", csys = "equat", geoc = True):
    d = date(d, m, y, ut)
    n = planetNum(p)
    if (n == -1):
        print("Please enter a valid planet name.")
        return 0, 0, 0
    if (n == 0):
        moon_arr = correctCoords(getPlanPos(n, d), sunPos(d, csys="pertur")[3], ct, csys)
        return moon_arr
    if (geoc):
        return getPlanPosGeo(n, d, ct, csys)

    return getPlanPos(n, d, ct)

"""
    Returns the correct version of the coordinates as requested by the user
"""
def correctCoords(arr, oblecl, ct, csys):
    if (csys == "equat"):
        arr = Coordinates.ecEq(arr[0], arr[1], arr[2], oblecl)
    if (ct == "sphere"):
        arr = Coordinates.sphere(arr[0], arr[1], arr[2])
    return arr

"""
    Returns the number of the planet identifier used by the program
"""
def planetNum(p):
    if (p.lower() == "moon"):
        return 0
    elif (p.lower() == "mercury"):
        return 1
    elif (p.lower() == "venus"):
        return 2
    elif (p.lower() == "mars"):
        return 3
    elif (p.lower() == "jupiter"):
        return 4
    elif (p.lower() == "saturn"):
        return 5
    return -1


"""
    Returns heliocentric ecliptic coordinates for a given planet
    of either rectangle or spherical form
"""
def getPlanPos(n, d = 0.0, ct = "rect"):
    N, i, w, a, e, M = OrbitalElements.getPlan(n,d)

    E = eccentricAnomaly(M, e)
    r, v = distAnomaly(a, E, e)

    # Ecliptic Coords
    xeclip, yeclip, zeclip = Coordinates.EclipticCoords(r, N, v, w, i)

    # Convert to longitude, latitude, distance
    long, lat, r = Coordinates.sphere(xeclip, yeclip, zeclip)

    # Perturbations
    long, lat, r = Perturbations.perturPlanets(n, d, M, N, w, long, lat, r)

    if (ct == "sphere"):
        return long, lat, r
    
    xeclip, yeclip, zeclip = Coordinates.rect(long, lat, r)
    return xeclip, yeclip, zeclip


"""
    Given a date value and a planet, returns the geocentric coordinates of that planet.
"""
def getPlanPosGeo(n, d, ct = "rect", csys = "equat"):
    sun_arr = sunPos(d)
    plan_arr = getPlanPos(n, d)
    plan_geo = togeo(sun_arr, plan_arr)
    return correctCoords(plan_geo, sun_arr[3], ct, csys)


"""
    Converts to geocentric position using the sun's rectangular coordinates
"""
def togeo(sunCoords, planCoords):
    xgeoc = sunCoords[0] + planCoords[0]
    ygeoc = sunCoords[1] + planCoords[1]
    zgeoc = sunCoords[2] + planCoords[2]
    return xgeoc, ygeoc, zgeoc


"""
    Calculates the eccentric anomaly defined by Kepler's equation.
    The eccentric anomaly is an angle for calculating the position in
    an elliptic orbit
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
    Calculates the planet's distance from the Sun and the angle from
    the perihelion to the planet (as seen from the sun).
    Ther perihelion is the point of distance in a planet's orbit that
    is closest to the sun
"""
def distAnomaly(a, E, e):
    # Compute rectangular coords in orbital plane
    x = a * (Trig.cosd(E) - e)
    y = a * sqrt(1 - e * e) * Trig.sind(E)
    
    # Convert to Distance and true anomaly
    r = sqrt(x*x + y*y) 
    v = Trig.rev(Trig.atan2d(y, x))
    return r, v