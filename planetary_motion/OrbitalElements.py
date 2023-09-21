from . import Trig

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


# Assumes input will be between 0-5 for n
def getPlan(n,d):
    if n == 0:
        return moonElm(d)
    if n == 1:
        return mercuryElm(d)
    elif n == 2:
        return venusElm(d)
    elif n == 3:
        return marsElm(d)
    elif n == 4:
        return jupiterElm(d)
    return saturnElm(d)

# We treat N = 0, i = 0 in the sun position calculation
def sunElm(d):
    # Longitude of Perihelion (deg)
    w = Trig.rev(282.9404 + 4.70935E-5 * d)
    # Mean distance, a.u.
    a = 1.0
    # Eccentricity
    e = 0.016709 - 1.151E-9 * d  
    # Mean anomoly (deg)
    M = Trig.rev(356.0470 + 0.9856002585 * d)
    return w, a, e, M

def moonElm(d):
    # Longitude of Ascensding Node
    N = Trig.rev(125.1228 - 0.0529538083  * d)
    # Inclination    
    i = 5.1454
    # Argument of Perigee
    w = Trig.rev(318.0634 + 0.1643573223  * d)   
    a = 60.2666                                
    e = 0.054900                                
    M = Trig.rev(115.3654 + 13.0649929509 * d) 
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

