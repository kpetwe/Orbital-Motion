from math import *

"""
    Calculates planetary positions of planets with an error of at most 1-2 arc minutes.
    All formulae used in this file are courtesy of https://www.stjarnhimlen.se/comp/tutorial.html#4.
    Contains methods for calculating between coordinate systems (elliptical and equatorial)
    and for standard heliocentric and geocentric models.
"""
class OrbitalMotion:
    """
        Computes the time value "d" using a date in the year.
    """
    # d = 0.0 at 2000 Jan 0.0 or 1999 Dec 31.0 0:00 UT
    # D = date, y = year in 4 digits, M = month (1-12)
    def date(D, m, y, UT = 0.0):
        d = 367 * y - 7 * (y + (m + 9) // 12) // 4 - 3 * ((y + (m - 9) // 7 ) // 100 + 1) // 4 + 275 * m // 9 + D - 730515
        return d + UT / 24.0
    
    """
        Longitudinal Precission: since all calculations are relative to the 
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
        L = rev(w + M)
        # Auxiliary Angle (Eccentric Anomaly)
        E = M + (180/pi) * e * sind(M) * (1 + e * cosd(M))
        
        r, v = OrbitalMotion.distAnomaly(a, E, e)
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

        E = OrbitalMotion.eccentricAnomaly(M, e)
        r, v = OrbitalMotion.distAnomaly(a, E, e)

        # Ecliptic Coords
        xeclip, yeclip, zeclip = Coordinates.EclipticCoords(r, N, v, w, i)

        # Convert to longitude, latitude, distance
        long, lat, r = Coordinates.SphericalCoords(xeclip, yeclip, zeclip)

        # Perturbations
        long, lat, r = Perturbations.perturPlanets(n, d, M, N, w, long, lat, r)
    
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
        plan_geo = OrbitalMotion.togeo(sun_arr, plan_arr)
        plan_geo_eq = Coordinates.ecEq(plan_geo[0], plan_geo[1], plan_geo[2], sun_arr[11])
        return plan_geo_eq


    """
        TODO: ADD description for this method
    """
    def eccentricAnomaly(M, e):
        # Eccentric Anomaly Calculation
        E0 = M + (180/pi) * e * sind(M) * (1 + e * cosd(M))

        # Since Eccentricity of the Moon's orbit is larger than of the Earth's,
        # Our initial approximation will have a greater error than 0.005 degrees
        E1 = E0 - (E0 - (180/pi) * e * sind(E0) - M) / (1 - e * cosd(E0))
        while (abs(E0 - E1) > 0.0005):
            E0 = E1
            E1 = E0 - (E0 - (180/pi) * e * sind(E0) - M) / (1 - e * cosd(E0))
        return E1
    

    """
        TODO: Add description for this method
    """
    def distAnomaly(a, E, e):
        # Compute rectangular coords in orbital plane
        x = a * (cosd(E) - e)
        y = a * sqrt(1 - e * e) * sind(E)
        
        # Convert to Distance and true anomaly
        r = sqrt(x*x + y*y) 
        v = rev(atan2d(y, x))
        return r, v
    


"""
    Series of methods that allow us to convert between various coordinate systems when
    calculating / presenting orbit data.
    Conversions mentioned include:
        Rectangular - Spherical coordinates
        Equatorial - Ecliptic coordinates
    
    This class also allows us to derive the actual Ecliptic and Spherical coordinates
    from the orbital elements.
"""
class Coordinates:
    """
        Derivations of the Ecliptical and Spherical coordinates using the orbital elements.
    """
    def EclipticCoords(r, N, v, w, i):
        xeclip = r * ( cosd(N) * cosd(v+w) - sind(N) * sind(v+w) * cosd(i) )
        yeclip = r * ( sind(N) * cosd(v+w) + cosd(N) * sind(v+w) * cosd(i) )
        zeclip = r * sind(v+w) * sind(i)
        return xeclip, yeclip, zeclip

    def SphericalCoords(xeclip, yeclip, zeclip):
        long =  atan2d(yeclip, xeclip)
        lat  =  atan2d(zeclip, sqrt(xeclip * xeclip + yeclip * yeclip))
        r    =  sqrt(xeclip * xeclip + yeclip * yeclip + zeclip * zeclip)
        return long, lat, r
    

    """
        Coordinate conversions between rectangular coordinates and spherical cooinates.
        
        Conversions can also apply for ecliptic and horizontal coordinates. 
        Replace RA, Decl with:
            Long, Lat (Ecliptic)
            Azimuth, Altitude (Horizontal)
    """
    # Spherical -> Rectangular Coordinates
    def rect(RA, Decl, r = 1):
        x = r * cosd(RA) * cosd(Decl)
        y = r * sind(RA) * cosd(Decl)
        z = r * sind(Decl)
        return x, y, z

    # Rectangular -> Spherical Coordinates
    def sphere(x, y, z):
        r = sqrt(x * x + y * y + z * z)
        RA = atan2d(y, x)
        Decl = atan2d(z, sqrt(x * x + y * y))
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
        yeclip = yequat * cosd(-oblecl) - zequat * sind(-oblecl)
        zeclip = yequat * sind(-oblecl) + zequat * cosd(-oblecl)
        return xeclip, yeclip, zeclip
        

    # Ecliptic -> Equatorial Conversion
    def ecEq(xeclip, yeclip, zeclip, oblecl = 23.4):
        xequat = xeclip
        yequat = yeclip * cosd(oblecl) - zeclip * sind(oblecl)
        zequat = yeclip * sind(oblecl) + zeclip * cosd(oblecl)
        return xequat, yequat, zequat
    


"""
    A series of methods which return the orbital elements of 
    the Sun, Moon, and planets Mercury - Saturn, given a date value d.

    Values returned are:

    N = longitude of the ascending node
    i = inclination to the ecliptic (plane of the Earth's orbit)
    w = argument of perihelion
    a = semi-major axis, or mean distance from Sun
    e = eccentricity (0=circle, 0-1=ellipse, 1=parabola)
    M = mean anomaly (0 at perihelion; increases uniformly with time)

"""
class OrbitalElements:

    # Assumes input will be between 0-5 for n
    def getPlan(n,d):
        if n == 0:
            return OrbitalElements.moonElm(d)
        if n == 1:
            return OrbitalElements.mercuryElm(d)
        elif n == 2:
            return OrbitalElements.venusElm(d)
        elif n == 3:
            return OrbitalElements.marsElm(d)
        elif n == 4:
            return OrbitalElements.jupiterElm(d)
        return OrbitalElements.saturnElm(d)

    # We treat N = 0, i = 0 in the sun position calculation
    def sunElm(d):
        # Longitude of Perihelion (deg)
        w = rev(282.9404 + 4.70935E-5 * d)
        # Mean distance, a.u.
        a = 1.0
        # Eccentricity
        e = 0.016709 - 1.151E-9 * d  
        # Mean anomoly (deg)
        M = rev(356.0470 + 0.9856002585 * d)
        return w, a, e, M

    def moonElm(d):
        # Longitude of Ascensding Node
        N = rev(125.1228 - 0.0529538083  * d)
        # Inclination    
        i = 5.1454
        # Argument of Perigee
        w = rev(318.0634 + 0.1643573223  * d)   
        a = 60.2666                                
        e = 0.054900                                
        M = rev(115.3654 + 13.0649929509 * d) 
        return N, i, w, a, e, M 
    
    def mercuryElm(d):
        N =  48.3313 + 3.24587E-5   * d    
        i =   7.0047 + 5.00E-8      * d
        # Argument of Perehelion    
        w =  29.1241 + 1.01444E-5   * d
        # Semi-major axis    
        a = 0.387098                              
        e = 0.205635 + 5.59E-10     * d   
        M = 168.6562 + 4.0923344368 * d
        return N, i, w, a, e, M 

    def venusElm(d):
        N =  76.6799 + 2.46590E-5   * d
        i =   3.3946 + 2.75E-8      * d
        w =  54.8910 + 1.38374E-5   * d
        a = 0.723330
        e = 0.006773     - 1.302E-9  * d
        M =  48.0052 + 1.6021302244 * d
        return N, i, w, a, e, M

    def marsElm(d):
        N =  49.5574 + 2.11081E-5 * d
        i = 1.8497 - 1.78E-8 * d
        w = 286.5016 + 2.92961E-5 * d
        a = 1.523688  
        e = 0.093405 + 2.516E-9 * d
        M =  18.6021 + 0.5240207766 * d
        return N, i, w, a, e, M

    def jupiterElm(d):
        N = 100.4542 + 2.76854E-5 * d
        i = 1.3030 - 1.557E-7 * d
        w = 273.8777 + 1.64505E-5 * d
        a = 5.20256
        e = 0.048498 + 4.469E-9 * d
        M =  19.8950 + 0.0830853001 * d
        return N, i, w, a, e, M

    def saturnElm(d):
        N = 113.6634 + 2.38980E-5 * d
        i = 2.4886 - 1.081E-7 * d
        w = 339.3939 + 2.97661E-5 * d
        a = 9.55475  
        e = 0.055546 - 9.499E-9 * d
        M = 316.9670 + 0.0334442282 * d 
        return N, i, w, a, e, M
    

"""
    A series of functions for calculating the perturbations of the moon,
    jupiter, and saturn, as they significantly affect the accuracy of the
    orbital calculation.
"""
class Perturbations:
    """
        TODO: ADD comment
    """
    def perturPlanets(n, d, M, N, w, long, lat, r):
        # Perturbations
        if (n == 4):
            sat_arr = OrbitalElements.getPlan(5,d)
            long += Perturbations.jupPertur(M, sat_arr[5])
        elif (n == 5):
            jup_arr = OrbitalElements.getPlan(4, d)
            diff_arr = Perturbations.satPertur(jup_arr[5], M)
            long += diff_arr[0]
            lat += diff_arr[1]
        elif (n == 0):
            sun_arr = OrbitalMotion.sunPos(d)
            diff_arr = Perturbations.moonPerturSetup(sun_arr[9], N, w, sun_arr[10], M)
            long += diff_arr[0]
            lat += diff_arr[1]
            r += diff_arr[2]
        return long, lat, r
        
    """
        Moon perturbation calculations
    """
    def moonPertur(Ms, Mm, D, F):
        dif_lon = 0.0
        dif_lon -= 1.274 * sind(Mm - 2*D)
        dif_lon += 0.658 * sind(2*D)
        dif_lon -= 0.186 * sind(Ms)
        dif_lon -= 0.059 * sind(2*Mm - 2*D)
        dif_lon -= 0.057 * sind(Mm - 2*D + Ms)
        dif_lon += 0.053 * sind(Mm + 2*D)
        dif_lon += 0.046 * sind(2*D - Ms)
        dif_lon += 0.041 * sind(Mm - Ms)
        dif_lon -= 0.035 * sind(D)
        dif_lon -= 0.031 * sind(Mm + Ms)
        dif_lon -= 0.015 * sind(2*F - 2*D)
        dif_lon += 0.011 * sind(Mm - 4*D)

        dif_lat = 0.0
        dif_lat -= 0.173 * sind(F - 2*D)
        dif_lat -= 0.055 * sind(Mm - F - 2*D)
        dif_lat -= 0.046 * sind(Mm + F - 2*D)
        dif_lat += 0.033 * sind(F + 2*D)
        dif_lat += 0.017* sind(2*Mm + F)

        dif_r = -0.58 * cosd(Mm - 2*D) -0.46 * cosd(2*D)
        return dif_lon, dif_lat, dif_r
    
    # Takes orbital calclation values and calculates values needed
    # for the moon perturbation
        # Ls  = Sun's mean longitude
        # Ms = Sun's mean anomaly 
        # All other values are for the moon
    def moonPerturSetup(Ls, N, w, Ms, Mm):
        # Moon's Mean Longitude 
        Lm = rev(N + w + Mm)
        # Moon's mean Elongation
        D = rev(Lm - Ls)
        # Moon's Argument of Latitude
        F = rev(Lm - N)
        # TODO: Rewrite so that this returns moonPertur(Ms, Mm, D, F) instead
        return Perturbations.moonPertur(Ms, Mm, D, F)
    
    """
        Jupiter and Saturn perturbation calculations

        Mj = mean anomaly of Jupiter
        Ms = mean anomaly of Saturn
    """
    def jupPertur(Mj, Ms):
        diff_long = 0.0
        diff_long -= 0.332 * sind(2 * Mj - 5 * Ms - 67.6)
        diff_long -= 0.056 * sind(2 * Mj - 2 * Ms + 21)
        diff_long += 0.042 * sind(3 * Mj - 5 * Ms + 21)
        diff_long -= 0.036 * sind(Mj - 2 * Ms)
        diff_long += 0.022 * cosd(Mj - Ms)
        diff_long += 0.023 * sind(2 * Mj - 3 * Ms + 52)
        diff_long -= 0.016 * sind(Mj - 5 * Ms - 69)
        return diff_long
    
    def satPertur(Mj, Ms):
        diff_long = 0.0
        diff_long += 0.812 * sind(2 * Mj - 5 * Ms - 67.6)
        diff_long -= 0.229 * cosd(2 * Mj - 4 * Ms - 2)
        diff_long += 0.119 * sind(Mj - 2*Ms - 3)
        diff_long += 0.046 * sind(2 * Mj - 6 * Ms - 69)
        diff_long += 0.014 * sind(Mj - 3 * Ms + 32)

        diff_lat = 0.0
        diff_lat -= 0.020 * cosd(2 * Mj - 4 * Ms - 2)
        diff_lat += 0.018 * sind(2 * Mj - 6 * Ms - 49)
        return diff_long, diff_lat




"""
    A series of miscellaneous functions for 
    calculating trigonmetric values using degrees instead of radians. 
    Simplifies planetary calculations.
"""

# Methods/Consts covered by Python:
    # pi = 3.141592653589793
    # radians(deg) -> degrees to radians
    # degrees(rad) -> radians to degrees
    # np.cbrt(x) -> cube root function

# Conversions between radians and degrees for various functions
sind = lambda x: sin(radians(x))
cosd = lambda x: cos(radians(x))
tand = lambda x: tan(radians(x))
asind = lambda x: degrees(asin(x))
acosd = lambda x: degrees(acos(x))
atand = lambda x: degrees(atan(x))
atan2d = lambda y,x: degrees(atan2(y,x))

# Reduces angle to be between 0 and 360 degrees
def rev(x):
    return x - floor(x/360.0) * 360.0