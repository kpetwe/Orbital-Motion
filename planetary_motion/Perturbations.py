from . import OrbitalElements
from . import OrbitalMotion
from . import Trig

"""
    A series of functions for calculating the perturbations of the moon,
    jupiter, and saturn, as they significantly affect the accuracy of the
    orbital calculation.
"""
"""
    Modifies final calculation of longitude, latitude, and distance
    by accounting for perturbations between Saturn and Jupiter or the 
    Moon and the Sun due to gravity
"""
def perturPlanets(n, d, M, N, w, long, lat, r):
    # Perturbations
    if (n == 4):
        sat_arr = OrbitalElements.getPlan(5, d)
        # "sat_arr[5]"
        long += jupPertur(M, sat_arr[5])
    elif (n == 5):
        jup_arr = OrbitalElements.getPlan(4, d)
        diff_arr = satPertur(jup_arr[5], M)
        long += diff_arr[0]
        lat += diff_arr[1]
    elif (n == 0):
        sun_arr = OrbitalMotion.sunPos(d, csys = "pertur")
        diff_arr = moonPerturSetup(sun_arr[0], N, w, sun_arr[1], M)
        long += diff_arr[0]
        lat += diff_arr[1]
        r += diff_arr[2]
    return Trig.rev(long), Trig.rev(lat), r
    
"""
    Moon perturbation calculations
"""
def moonPertur(Ms, Mm, D, F):
    dif_lon = 0.0
    dif_lon -= 1.274 * Trig.sind(Mm - 2*D)
    dif_lon += 0.658 * Trig.sind(2*D)
    dif_lon -= 0.186 * Trig.sind(Ms)
    dif_lon -= 0.059 * Trig.sind(2*Mm - 2*D)
    dif_lon -= 0.057 * Trig.sind(Mm - 2*D + Ms)
    dif_lon += 0.053 * Trig.sind(Mm + 2*D)
    dif_lon += 0.046 * Trig.sind(2*D - Ms)
    dif_lon += 0.041 * Trig.sind(Mm - Ms)
    dif_lon -= 0.035 * Trig.sind(D)
    dif_lon -= 0.031 * Trig.sind(Mm + Ms)
    dif_lon -= 0.015 * Trig.sind(2*F - 2*D)
    dif_lon += 0.011 * Trig.sind(Mm - 4*D)

    dif_lat = 0.0
    dif_lat -= 0.173 * Trig.sind(F - 2*D)
    dif_lat -= 0.055 * Trig.sind(Mm - F - 2*D)
    dif_lat -= 0.046 * Trig.sind(Mm + F - 2*D)
    dif_lat += 0.033 * Trig.sind(F + 2*D)
    dif_lat += 0.017* Trig.sind(2*Mm + F)

    dif_r = -0.58 * Trig.cosd(Mm - 2*D) -0.46 * Trig.cosd(2*D)
    return dif_lon, dif_lat, dif_r

# Takes orbital calclation values and calculates values needed
# for the moon perturbation
    # Ls  = Sun's mean longitude
    # Ms = Sun's mean anomaly 
    # All other values are for the moon
def moonPerturSetup(Ls, N, w, Ms, Mm):
    # Moon's Mean Longitude 
    Lm = Trig.rev(N + w + Mm)
    # Moon's mean Elongation
    D = Trig.rev(Lm - Ls)
    # Moon's Argument of Latitude
    F = Trig.rev(Lm - N)
    # TODO: Rewrite so that this returns moonPertur(Ms, Mm, D, F) instead
    return moonPertur(Ms, Mm, D, F)

"""
    Jupiter and Saturn perturbation calculations

    Mj = mean anomaly of Jupiter
    Ms = mean anomaly of Saturn
"""
def jupPertur(Mj, Ms):
    diff_long = 0.0
    diff_long -= 0.332 * Trig.sind(2 * Mj - 5 * Ms - 67.6)
    diff_long -= 0.056 * Trig.sind(2 * Mj - 2 * Ms + 21)
    diff_long += 0.042 * Trig.sind(3 * Mj - 5 * Ms + 21)
    diff_long -= 0.036 * Trig.sind(Mj - 2 * Ms)
    diff_long += 0.022 * Trig.cosd(Mj - Ms)
    diff_long += 0.023 * Trig.sind(2 * Mj - 3 * Ms + 52)
    diff_long -= 0.016 * Trig.sind(Mj - 5 * Ms - 69)
    return diff_long

def satPertur(Mj, Ms):
    diff_long = 0.0
    diff_long += 0.812 * Trig.sind(2 * Mj - 5 * Ms - 67.6)
    diff_long -= 0.229 * Trig.cosd(2 * Mj - 4 * Ms - 2)
    diff_long += 0.119 * Trig.sind(Mj - 2*Ms - 3)
    diff_long += 0.046 * Trig.sind(2 * Mj - 6 * Ms - 69)
    diff_long += 0.014 * Trig.sind(Mj - 3 * Ms + 32)

    diff_lat = 0.0
    diff_lat -= 0.020 * Trig.cosd(2 * Mj - 4 * Ms - 2)
    diff_lat += 0.018 * Trig.sind(2 * Mj - 6 * Ms - 49)
    return diff_long, diff_lat